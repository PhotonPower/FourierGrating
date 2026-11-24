// FourierGrating - CAD slicer with Nurbforge (Plan section "M14: CAD-Slicer"): STEP -> FG
// coordinates (micrometres, axes) -> checks -> bands between critical heights -> homogeneous /
// Contour2D / SectionedLayer. The only module that includes Nurbforge headers.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <mutex>
#include <nurbforge/analysis/shell_section.hpp>
#include <nurbforge/core/error.hpp>
#include <nurbforge/geometry/transform.hpp>
#include <nurbforge/io/step.hpp>
#include <nurbforge/topology/boolean.hpp>
#include <sstream>

#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/cyl/angular.hpp"
#include "fg/geometry/line_section.hpp"
#include "fg/geometry/section.hpp"

namespace fg::cad {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string fmt(const char* f, double a, double b = 0, double c = 0) {
  char buf[256];
  std::snprintf(buf, sizeof buf, f, a, b, c);
  return buf;
}

// ------------------------------------------------------------------ model in FG coordinates

/// A body: shell in FG coordinates (micrometres, x, y = stacking, z), its material, its name.
struct Body {
  nf::Shell shell;
  std::string name, material;
  std::array<double, 6> box{};  ///< xmin, ymin, zmin, xmax, ymax, zmax
};

/// Bounding box of the trimmed body: the extreme critical heights along x, y, z (the extrema of a
/// coordinate lie on vertices, edge extrema, face extrema inside the trim or planar faces normal to
/// the axis). The control points of the untrimmed carrier surfaces would overestimate it, e.g. a
/// lid trimmed out of a larger rectangle, and shift the default y range and cell origin.
std::array<double, 6> bounds(const nf::Shell& s) {
  std::array<double, 6> b{};
  const nf::Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (std::size_t i = 0; i < 3; ++i) {
    const std::vector<double> h = nf::critical_heights(s, axes[i]);
    FG_REQUIRE(!h.empty(), "cad: body without vertices, edges or faces");
    b[i] = h.front();
    b[i + 3] = h.back();
  }
  return b;
}

nf::Vec3 axis(const std::string& a) {
  FG_REQUIRE(a.size() == 2 && (a[0] == '+' || a[0] == '-'),
             "cad.axes: axes are \"+X\", \"-Y\", ... (got '" + a + "')");
  const double s = a[0] == '+' ? 1.0 : -1.0;
  switch (a[1]) {
    case 'X':
    case 'x':
      return {s, 0, 0};
    case 'Y':
    case 'y':
      return {0, s, 0};
    case 'Z':
    case 'z':
      return {0, 0, s};
    default:
      throw Error("cad.axes: unknown axis '" + a + "'");
  }
}

nf::Shell transformed(const nf::Shell& s, const nf::Transform& T) {
  std::vector<nf::TrimmedSurface> faces;
  for (const auto& f : s.faces())
    faces.emplace_back(f.surface().transformed(T), f.outer(), f.holes());
  return nf::Shell::sew(std::move(faces), 1e-6);
}

std::string resolveBackground(const CadOptions& o, const Stack& st) {
  if (o.background == "superstrate") return st.superstrate;
  if (o.background == "substrate") return st.substrate;
  return o.background;
}

/// Bodies in FG coordinates, cell-relative (x - x0, z - z0); findings for reading, orientation,
/// units, assignment. Returns false if the model cannot be used further.
struct Model {
  std::vector<Body> bodies;
  std::array<double, 2> origin{0, 0};
  std::array<double, 2> yRange{0, 0};
  std::string background;
  // 1D gratings (Plan section "M23"): the cutting plane and the profiles of the bodies in (x, y)
  bool oneD = false;
  double zCut = kNaN;
  // solids of revolution (Plan "M20: Schichttypen"): axis at x = z = 0 after the shift, largest radius
  bool axisym = false;
  int rot = 0;  ///< P-fold rotational symmetry (M21): axis at x = z = 0 as for axisym
  double Rs = 0;
  std::shared_ptr<const std::vector<geometry::ProfileBody>> profiles;
};

bool loadModel(const CadOptions& o, const Stack& st, Model& m, Report& r) {
  // 2D: sections in the (x, z) cell; 1D (dz = 0): z-invariant bodies cut once at z_cut (Plan "M23")
  m.axisym = o.axisymmetric;
  m.rot = o.rotation;
  const bool axial = m.axisym || m.rot > 0;  // about an axis, no period
  FG_REQUIRE(axial || st.dx_um > 0, "cad: period_um.x must be positive");
  m.oneD = !axial && st.dz_um <= 0;
  auto error = [&](const std::string& code, const std::string& msg,
                   std::array<double, 3> w = {kNaN, kNaN, kNaN}) {
    r.errors.push_back({code, msg, w});
    r.ok = false;
  };
  // ---- read (micrometres; "unit" corrects a wrong declaration)
  std::vector<nf::NamedShell> named;
  double declared = 1e-3;
  try {
    nf::StepReadOptions ro;
    ro.unit = nf::StepLengthUnit::Micrometre;
    ro.tolerance = 1e-6;
    ro.file_unit_metres = &declared;
    named =
        o.stepText.empty() ? nf::read_step_named(o.file, ro) : nf::from_step_named(o.stepText, ro);
  } catch (const nf::GeometryError& e) {
    error("STEP_UNSUPPORTED", e.what());
    return false;
  } catch (const std::runtime_error& e) {
    error("STEP_UNSUPPORTED", e.what());
    return false;
  }
  if (named.empty()) {
    error("STEP_UNSUPPORTED", "no solid (MANIFOLD_SOLID_BREP / BREP_WITH_VOIDS) in the STEP data");
    return false;
  }
  double scale = 1.0;  // file coordinates were converted from `declared` to um
  if (o.unit == "mm")
    scale = 1e-3 / declared;
  else if (o.unit == "um")
    scale = 1e-6 / declared;
  else
    FG_REQUIRE(o.unit == "auto", "cad.unit must be \"auto\", \"mm\" or \"um\"");
  // ---- axes: FG x, y from CAD axes, z = x cross y (right-handed; mirrored axes rejected)
  const nf::Vec3 ex = axis(o.axisX), ey = axis(o.axisY);
  FG_REQUIRE(std::abs(nf::dot(ex, ey)) < 0.5, "cad.axes: x and y must be different CAD axes");
  const nf::Vec3 ez = nf::cross(ex, ey);
  if (!o.axisZ.empty()) {
    const nf::Vec3 g = axis(o.axisZ);
    FG_REQUIRE(nf::dot(g, ez) > 0.5,
               "cad.axes: z = " + o.axisZ +
                   " gives a left-handed (mirrored) system; mirroring is not allowed "
                   "(it changes chiral structures); FG z must be x cross y");
  }
  // rows: FG axes in CAD coordinates (rotation), uniform scale
  const nf::Transform R({scale * ex.x, scale * ex.y, scale * ex.z, 0, scale * ey.x, scale * ey.y,
                         scale * ey.z, 0, scale * ez.x, scale * ez.y, scale * ez.z, 0});
  // ---- assignments
  std::vector<bool> used(o.bodies.size(), false);
  std::array<double, 6> all{1e300, 1e300, 1e300, -1e300, -1e300, -1e300};
  for (std::size_t k = 0; k < named.size(); ++k) {
    std::string material;
    for (std::size_t a = 0; a < o.bodies.size(); ++a) {
      const auto& b = o.bodies[a];
      if ((b.index >= 0 && static_cast<std::size_t>(b.index) == k) ||
          (!b.name.empty() && b.name == named[k].name)) {
        material = b.material;
        used[a] = true;
        break;
      }
    }
    nf::Shell sh = transformed(named[k].shell, R);
    const auto bb = bounds(sh);
    if (material.empty()) {
      error("BODY_UNASSIGNED",
            "body " + std::to_string(k) + " ('" + named[k].name +
                "') has no material (cad.bodies: name or index)",
            {0.5 * (bb[0] + bb[3]), 0.5 * (bb[1] + bb[4]), 0.5 * (bb[2] + bb[5])});
      continue;
    }
    FG_REQUIRE(st.materials.count(material), "cad.bodies: unknown material '" + material + "'");
    for (int i = 0; i < 3; ++i) {
      all[static_cast<std::size_t>(i)] =
          std::min(all[static_cast<std::size_t>(i)], bb[static_cast<std::size_t>(i)]);
      all[static_cast<std::size_t>(i + 3)] =
          std::max(all[static_cast<std::size_t>(i + 3)], bb[static_cast<std::size_t>(i + 3)]);
    }
    m.bodies.push_back({std::move(sh), named[k].name, material, bb});
  }
  for (std::size_t a = 0; a < o.bodies.size(); ++a)
    if (!used[a])
      error("BODY_UNASSIGNED", "cad.bodies entry " + std::to_string(a) + " ('" + o.bodies[a].name +
                                   "') matches no body of the STEP file");
  if (m.bodies.empty()) return false;
  // ---- cell origin, y range, translation into the cell
  // axisymmetric: the origin is the axis (default: centre of the bounding box); rotation: a first guess, the
  // centroid of a section follows (centreOnAxis)
  m.origin = o.cellOrigin_um ? *o.cellOrigin_um
             : axial       ? std::array<double, 2>{0.5 * (all[0] + all[3]), 0.5 * (all[2] + all[5])}
                           : std::array<double, 2>{all[0], all[2]};
  if (m.oneD) m.origin[1] = 0;  // z stays FG z (z_cut_um)
  const nf::Transform shift = nf::Transform::translation({-m.origin[0], 0, -m.origin[1]});
  for (auto& b : m.bodies) {
    b.shell = transformed(b.shell, shift);
    b.box = bounds(b.shell);
    for (int i : {0, 2, 3, 5}) m.Rs = std::max(m.Rs, std::abs(b.box[static_cast<std::size_t>(i)]));
  }
  m.yRange = o.yRange_um ? *o.yRange_um : std::array<double, 2>{all[1], all[4]};
  FG_REQUIRE(m.yRange[1] > m.yRange[0], "cad.y_range_um: empty range");
  m.background = resolveBackground(o, st);
  FG_REQUIRE(st.materials.count(m.background),
             "cad.background: unknown material '" + m.background + "'");
  // ---- units: dimensions far away from the period / wavelength (factor 1e3)
  // 1D: the extrusion length in z is arbitrary and not checked
  const double extent = m.oneD ? std::max(all[3] - all[0], all[4] - all[1])
                               : std::max({all[3] - all[0], all[4] - all[1], all[5] - all[2]}),
               period = axial ? (o.lambda_um > 0 ? o.lambda_um : extent) : std::max(st.dx_um, st.dz_um);
  const double ref = o.lambda_um > 0 ? std::max(period, o.lambda_um) : period;
  if (extent > 100 * ref || extent < 1e-3 * std::min(period, ref))
    error("UNIT_SUSPECT", fmt("bodies are %.6g um large for a period of %.6g um (factor 1000 off? "
                              "cad.unit \"mm\" / \"um\", or "
                              "write the STEP file with unit micrometre)",
                              extent, period));
  return true;
}

// ------------------------------------------------------------------ checks

/// Bezier pieces (x, y) of a profile curve (section with a plane z = const) in FG 3D coordinates.
geometry::Contour2D profileLoop(const nf::SectionLoop& loop) {
  geometry::Contour2D c;
  for (const auto& sc : loop.curves)
    for (const auto& b : sc.curve3d.decompose_bezier()) {
      FG_REQUIRE(b.degree() <= 9, "cad: profile curve of degree > 9");
      geometry::BezierSegment s;
      for (std::size_t j = 0; j < b.num_control_points(); ++j)
        s.P.push_back({b.control_points()[j].x, b.control_points()[j].y});
      if (b.is_rational()) s.w = b.weights();
      c.segments.push_back(std::move(s));
    }
  return c;
}

/// 1D gratings (Plan section "M23: z-Invarianz und Schnittebene"): every body extruded along z, the
/// cutting plane z_cut inside all of them, the profiles in (x, y).
void checkProfiles(const CadOptions& o, Model& m, Report& r) {
  auto error = [&](const std::string& code, const std::string& msg,
                   std::array<double, 3> w = {kNaN, kNaN, kNaN}) {
    r.errors.push_back({code, msg, w});
    r.ok = false;
  };
  double zlo = -1e300, zhi = 1e300;
  for (const auto& b : m.bodies) {
    zlo = std::max(zlo, b.box[2]);
    zhi = std::min(zhi, b.box[5]);
    if (!nf::is_extrusion_band(b.shell, {0, 0, 1}, b.box[2], b.box[5]))
      error("NOT_Z_INVARIANT",
            "body '" + b.name +
                "' is not extruded along z: a 1D grating needs z-invariant bodies (a line extruded "
                "along FG z, by default -CAD Y); use a 2D grating for this body",
            {0.5 * (b.box[0] + b.box[3]) + m.origin[0], 0.5 * (b.box[1] + b.box[4]),
             0.5 * (b.box[2] + b.box[5])});
  }
  if (!r.ok) return;
  if (zlo >= zhi) {
    error("OUT_OF_RANGE", fmt("the bodies have no common z range (1D: one plane z = z_cut_um must cut all "
                              "of them; overlap %.6g ... %.6g um)",
                              zlo, zhi));
    return;
  }
  m.zCut = o.zCut_um ? *o.zCut_um : 0.5 * (zlo + zhi);
  for (const auto& b : m.bodies)
    if (!(m.zCut > b.box[2] + o.criticalTol_um && m.zCut < b.box[5] - o.criticalTol_um))
      error("OUT_OF_RANGE", "z_cut_um = " + fmt("%.6g um does not cut body '", m.zCut) + b.name +
                                fmt("' (z = %.6g ... %.6g um)", b.box[2], b.box[5]),
            {0.5 * (b.box[0] + b.box[3]) + m.origin[0], 0.5 * (b.box[1] + b.box[4]), m.zCut});
  if (!r.ok) return;
  auto profiles = std::make_shared<std::vector<geometry::ProfileBody>>();
  try {
    for (const auto& b : m.bodies) {
      nf::SectionOptions so;
      so.s_axis = {1, 0, 0};
      const nf::PlanarSection s = nf::section(b.shell, {{0, 0, m.zCut}, {0, 0, 1}}, so);
      geometry::ProfileBody p;
      p.material = b.material;
      p.name = b.name;
      for (const auto& loop : s.loops) p.loops.push_back(profileLoop(loop));
      profiles->push_back(std::move(p));
    }
  } catch (const std::exception& e) {
    error("SHELL_OPEN", std::string("profile section at z_cut_um failed: ") + e.what());
    return;
  }
  m.profiles = profiles;
}

void checkModel(const CadOptions& o, Model& m, Report& r) {
  auto error = [&](const std::string& code, const std::string& msg,
                   std::array<double, 3> w = {kNaN, kNaN, kNaN}) {
    r.errors.push_back({code, msg, w});
    r.ok = false;
  };
  for (auto& b : m.bodies) {
    const auto fe = b.shell.free_edges();
    if (!fe.empty()) {
      const auto& c = b.shell.edges()[fe.front()].curve;
      const nf::Vec3 p = c.evaluate(0.5 * (c.domain_start() + c.domain_end()));
      error("SHELL_OPEN",
            "body '" + b.name + "' is not closed (" + std::to_string(fe.size()) + " free edges)",
            {p.x + m.origin[0], p.y, p.z + m.origin[1]});
      continue;
    }
    if (!b.shell.is_consistently_oriented() || b.shell.volume(1e-6) < 0) {
      r.warnings.push_back(
          {"SHELL_ORIENTATION",
           "body '" + b.name + "' was not consistently oriented outwards (corrected)"});
      b.shell = b.shell.oriented();
    }
    if (b.box[1] < m.yRange[0] - o.criticalTol_um || b.box[4] > m.yRange[1] + o.criticalTol_um)
      error("OUT_OF_RANGE",
            "body '" + b.name + "' reaches " +
                fmt("y = %.6g ... %.6g um, outside y_range_um = [%.6g", b.box[1], b.box[4],
                    m.yRange[0]) +
                fmt(", %.6g]", m.yRange[1]),
            {0.5 * (b.box[0] + b.box[3]) + m.origin[0], b.box[1],
             0.5 * (b.box[2] + b.box[5]) + m.origin[1]});
  }
  if (!r.ok) return;
  if (m.oneD) checkProfiles(o, m, r);
}

// ------------------------------------------------------------------ sections

/// Bezier segments (x, z) of a section curve given in FG 3D coordinates.
std::vector<geometry::BezierSegment> bezier(const nf::NurbsCurve& c3) {
  std::vector<geometry::BezierSegment> out;
  for (const auto& b : c3.decompose_bezier()) {
    FG_REQUIRE(b.degree() <= 9, "cad: section curve of degree > 9");
    geometry::BezierSegment s;
    for (std::size_t j = 0; j < b.num_control_points(); ++j)
      s.P.push_back({b.control_points()[j].x, b.control_points()[j].z});
    if (b.is_rational()) s.w = b.weights();
    out.push_back(std::move(s));
  }
  return out;
}

/// Boundary samples of one section curve (chord deviation <= tol) with the 3D outward normals.
void samples(const nf::SectionCurve& c, double tol, BoundaryChain& chain) {
  std::vector<double> ts{c.curve3d.domain_start(), c.curve3d.domain_end()};
  for (int pass = 0; pass < 30; ++pass) {
    std::vector<double> next{ts.front()};
    bool refined = false;
    for (std::size_t k = 1; k < ts.size(); ++k) {
      const double a = ts[k - 1], b = ts[k], m = 0.5 * (a + b);
      const nf::Vec3 pa = c.curve3d.evaluate(a), pb = c.curve3d.evaluate(b),
                     pm = c.curve3d.evaluate(m);
      const double tx = pb.x - pa.x, tz = pb.z - pa.z, l = std::hypot(tx, tz);
      const double dev = l > 0 ? std::abs((pm.x - pa.x) * tz - (pm.z - pa.z) * tx) / l
                               : std::hypot(pm.x - pa.x, pm.z - pa.z);
      if (dev > tol || ts.size() < 9) {
        next.push_back(m);
        refined = true;
      }
      next.push_back(b);
    }
    ts = std::move(next);
    if (!refined) break;
  }
  for (std::size_t k = 0; k + 1 < ts.size(); ++k) {  // the end is the start of the next curve
    const nf::Vec3 p = c.curve3d.evaluate(ts[k]), n = c.normal(ts[k]);
    chain.samples.push_back({p.x, p.z, Vec3D(n.x, n.y, n.z), false});
  }
}

/// Section of one body at height y: regions (outer loops with their holes) and boundary samples.
void sectionOf(const Body& b, double y, double sampleTol, CrossSection& cs) {
  {
    if (y <= b.box[1] || y >= b.box[4]) return;
    nf::SectionOptions so;
    so.s_axis = {1, 0, 0};
    const nf::PlanarSection s = nf::section(b.shell, {{0, y, 0}, {0, 1, 0}}, so);
    std::vector<geometry::Region> outers;
    std::vector<geometry::Contour2D> holes;
    for (const auto& loop : s.loops) {
      geometry::Contour2D c;
      BoundaryChain chain;
      for (const auto& sc : loop.curves) {
        for (auto& seg : bezier(sc.curve3d)) c.segments.push_back(std::move(seg));
        samples(sc, sampleTol, chain);
      }
      cs.chains.push_back(std::move(chain));
      // the solid is on the left in (s, t) = (x, -z): area > 0 is an outer loop
      if (loop.area > 0) {
        geometry::Region r;
        r.outer = std::move(c);
        r.material = b.material;
        outers.push_back(std::move(r));
      } else {
        holes.push_back(std::move(c));
      }
    }
    // holes to the smallest outer loop containing them
    const double ftol = 1e-6 * std::max(b.box[3] - b.box[0], b.box[5] - b.box[2]);
    for (auto& r : outers) r.normalize(ftol);
    for (auto& h : holes) {
      const auto p = h.segments.front().P.front();
      int best = -1;
      double area = 1e300;
      for (std::size_t k = 0; k < outers.size(); ++k) {
        // a point slightly inside the hole would be better; the hole's start lies on its boundary,
        // which is inside the outer loop unless the loops touch
        if (outers[k].contains(p[0], p[1]) && std::abs(outers[k].area()) < area) {
          best = static_cast<int>(k);
          area = std::abs(outers[k].area());
        }
      }
      FG_REQUIRE(best >= 0, "cad: section hole without an enclosing outer loop");
      outers[static_cast<std::size_t>(best)].holes.push_back(std::move(h));
    }
    for (auto& r : outers) {
      r.normalize(ftol);  // again, with the holes
      cs.regions.push_back(std::move(r));
    }
  }
}

/// Section of all bodies at height y.
CrossSection sectionAt(const std::vector<Body>& bodies, double y, double sampleTol) {
  FG_PROF_SCOPE("cad.section");
  CrossSection cs;
  for (const Body& b : bodies) sectionOf(b, y, sampleTol, cs);
  return cs;
}

/// SectionLayer provider: sections of the bodies inside one band (y relative to the band bottom).
class NurbforgeSectionProvider final : public SectionProvider {
 public:
  NurbforgeSectionProvider(std::shared_ptr<const std::vector<Body>> b, double y0, double y1,
                           double tol, double cell)
      : bodies_(std::move(b)),
        y0_(y0),
        h_(y1 - y0),
        tol_(tol),
        margin_(1e-8 * std::max(cell, y1 - y0)) {}
  /// Never on a band boundary (a critical height, e.g. a lid): nf::section rejects planes within
  /// its tolerance of a face. Only endpoint evaluations (size probes, Radau nodes) are affected.
  CrossSection section(double y) const override {
    return sectionAt(*bodies_, y0_ + std::clamp(y, margin_, h_ - margin_), tol_);
  }
  double thickness_um() const override { return h_; }

 private:
  std::shared_ptr<const std::vector<Body>> bodies_;
  double y0_, h_, tol_, margin_;
};

/// Area fractions at mid height (bodies' regions, background rest).
std::vector<std::pair<std::string, double>> fractions(const CrossSection& cs, const std::string& bg,
                                                      double cell) {
  std::map<std::string, double> a;
  double sum = 0;
  for (const auto& r : cs.regions) {
    a[r.material] += std::abs(r.area());
    sum += std::abs(r.area());
  }
  a[bg] += cell - sum;
  std::vector<std::pair<std::string, double>> out;
  for (const auto& [m, v] : a) out.push_back({m, v / cell});
  return out;
}

/// Bands between the critical heights; merges thin bands (warning).
std::vector<std::array<double, 2>> bands(const Model& m, const CadOptions& o, Report& r) {
  std::vector<double> h{m.yRange[0], m.yRange[1]};
  for (const auto& b : m.bodies)
    for (double y : nf::critical_heights(b.shell, {0, 1, 0}, {o.criticalTol_um}))
      if (y > m.yRange[0] && y < m.yRange[1]) h.push_back(y);
  std::sort(h.begin(), h.end());
  std::vector<double> u;
  for (double y : h)
    if (u.empty() || y - u.back() > o.criticalTol_um) u.push_back(y);
  u.back() = m.yRange[1];
  std::vector<std::array<double, 2>> out;
  for (std::size_t k = 0; k + 1 < u.size(); ++k) out.push_back({u[k], u[k + 1]});
  for (std::size_t k = 0; k < out.size() && out.size() > 1;) {
    if (out[k][1] - out[k][0] >= o.minBand_um) {
      ++k;
      continue;
    }
    r.warnings.push_back(
        {"FEATURE_SMALL",
         fmt("band y = %.9g ... %.9g um thinner than min_band_um, merged", out[k][0], out[k][1])});
    if (k + 1 < out.size()) {
      out[k + 1][0] = out[k][0];
    } else {
      out[k - 1][1] = out[k][1];
    }
    out.erase(out.begin() + static_cast<std::ptrdiff_t>(k));
  }
  return out;
}

/// Centroid (x, z) of the regions of a section: first moments by Green's theorem, int x dA = 1/2 oint x^2 dz and
/// int z dA = -1/2 oint z^2 dx (outer loops counter-clockwise, holes clockwise; 48 Gauss nodes per half segment as
/// Contour2D::signedArea).
std::array<double, 2> centroid(const CrossSection& cs) {
  const auto& g = quadrature::gaussLegendre(48);
  double a = 0, mx = 0, mz = 0;
  for (const auto& reg : cs.regions) {
    std::vector<const geometry::Contour2D*> loops{&reg.outer};
    for (const auto& h : reg.holes) loops.push_back(&h);
    for (const auto* c : loops)
      for (const auto& seg : c->segments) {
        geometry::BezierSegment halves[2];
        seg.split(0.5, halves[0], halves[1]);
        for (const auto& h : halves)
          for (std::size_t l = 0; l < g.x.size(); ++l) {
            geometry::Point2 p, dp;
            h.eval(0.5 * (g.x[l] + 1), p, dp);
            const double w = 0.25 * g.w[l];
            a += w * (p[0] * dp[1] - p[1] * dp[0]);
            mx += w * p[0] * p[0] * dp[1];
            mz -= w * p[1] * p[1] * dp[0];
          }
      }
  }
  FG_REQUIRE(a != 0, "cad: section without area");
  return {mx / a, mz / a};
}

/// P-fold rotational symmetry (M21) without cad.cell_origin_um: the axis through the centroid of the lowest section
/// with boundaries (band mid heights; every section of a P-fold symmetric arrangement has its centroid on the axis,
/// P >= 2), the bodies shifted onto it.
void centreOnAxis(const CadOptions& o, Model& m) {
  if (m.rot == 0 || o.cellOrigin_um) return;
  Report dummy;
  for (const auto& band : bands(m, o, dummy)) {
    const CrossSection cs = sectionAt(m.bodies, 0.5 * (band[0] + band[1]), 1e-6 * m.Rs);
    if (cs.regions.empty()) continue;
    const auto c = centroid(cs);
    m.origin = {m.origin[0] + c[0], m.origin[1] + c[1]};
    const nf::Transform shift = nf::Transform::translation({-c[0], 0, -c[1]});
    m.Rs = 0;
    for (auto& b : m.bodies) {
      b.shell = transformed(b.shell, shift);
      b.box = bounds(b.shell);
      for (int i : {0, 2, 3, 5}) m.Rs = std::max(m.Rs, std::abs(b.box[static_cast<std::size_t>(i)]));
    }
    return;
  }
}

// ------------------------------------------------------------------ overlaps

/// Flattened section of one body (FG x, z) for the overlap test.
struct Flat {
  std::vector<geometry::Region> regions;
  std::vector<std::vector<geometry::Point2>> loops;
  std::array<double, 4> box{1e300, 1e300, -1e300, -1e300};  ///< xmin, zmin, xmax, zmax
};

Flat flatSection(const Body& b, double y, double tol) {
  CrossSection cs;
  sectionOf(b, y, tol, cs);
  Flat f;
  for (auto& r : cs.regions) {
    for (const auto* c : {&r.outer}) f.loops.push_back(c->flatten(tol));
    for (const auto& h : r.holes) f.loops.push_back(h.flatten(tol));
    f.regions.push_back(std::move(r));
  }
  for (const auto& l : f.loops)
    for (const auto& p : l) {
      f.box[0] = std::min(f.box[0], p[0]);
      f.box[1] = std::min(f.box[1], p[1]);
      f.box[2] = std::max(f.box[2], p[0]);
      f.box[3] = std::max(f.box[3], p[1]);
    }
  return f;
}

/// Length scale of the overlap probes: the smaller section extent.
double probeTol(const Flat& a, const Flat& b) {
  return std::min({a.box[2] - a.box[0], a.box[3] - a.box[1], b.box[2] - b.box[0], b.box[3] - b.box[1]});
}

/// Do the sections a and b + (sx, sz) overlap (proper boundary crossings or a probe of one inside the other)?
bool overlap(const Flat& a, const Flat& b, double sx, double sz) {
  if (a.loops.empty() || b.loops.empty()) return false;
  if (b.box[0] + sx >= a.box[2] || a.box[0] >= b.box[2] + sx || b.box[1] + sz >= a.box[3] ||
      a.box[1] >= b.box[3] + sz)
    return false;
  // rounding of touching bodies (e.g. faces of Boolean differences) is no crossing (v2.24)
  const double scale = std::max({a.box[2] - a.box[0], a.box[3] - a.box[1], b.box[2] - b.box[0], b.box[3] - b.box[1]});
  for (const auto& la : a.loops)
    for (const auto& lb : b.loops)
      for (std::size_t i = 0; i < la.size(); ++i) {
        const auto& p = la[i];
        const auto& q = la[(i + 1) % la.size()];
        for (std::size_t j = 0; j < lb.size(); ++j) {
          const geometry::Point2 r{lb[j][0] + sx, lb[j][1] + sz},
              t{lb[(j + 1) % lb.size()][0] + sx, lb[(j + 1) % lb.size()][1] + sz};
          if (geometry::properCrossing(p, q, r, t, 1e-9 * scale)) return true;
        }
      }
  // containment and collinear boundaries (e.g. boxes of the same z range): probes just inside every edge (v2.24;
  // the former test of a boundary point is ambiguous for touching bodies)
  const double eps = 1e-3 * probeTol(a, b);
  for (const auto& rb : b.regions)
    for (const auto& q : rb.interiorProbes(eps))
      for (const auto& ra : a.regions)
        if (ra.contains(q[0] + sx, q[1] + sz)) return true;
  for (const auto& ra : a.regions)
    for (const auto& q : ra.interiorProbes(eps))
      for (const auto& rb : b.regions)
        if (rb.contains(q[0] - sx, q[1] - sz)) return true;
  return false;
}

/// Do the boundaries of a and b + (sx, sz) share a piece longer than tol (collinear segments
/// within tol)? A body touching its own periodic image this way spans the period (e.g. a
/// z-invariant line extruded over exactly dz): its periodic union is a strip around the cell,
/// which closed contours cannot describe, and the shared edge is no material interface.
bool touches(const Flat& a, const Flat& b, double sx, double sz, double tol) {
  if (a.loops.empty() || b.loops.empty()) return false;
  if (b.box[0] + sx > a.box[2] + tol || a.box[0] > b.box[2] + sx + tol ||
      b.box[1] + sz > a.box[3] + tol || a.box[1] > b.box[3] + sz + tol)
    return false;
  for (const auto& la : a.loops)
    for (auto lb : b.loops) {
      for (auto& p : lb) p = {p[0] + sx, p[1] + sz};
      if (geometry::sharedBoundary(la, lb, tol)) return true;
    }
  return false;
}

/// BODY_OVERLAP / PERIODIC_OVERLAP from sections at three heights per band (the topology of the
/// sections does not change inside a band); "priority" and periodic_union resolve them with
/// Booleans.
void checkOverlaps(const CadOptions& o, const Stack& st, Model& m, Report& r);

geometry::LineSectionOptions lineOptions(const CadOptions& o) {
  geometry::LineSectionOptions lo;
  lo.periodicUnion = o.periodicUnion;
  lo.priority = o.overlap == "priority";
  return lo;
}

/// 1D (Plan section "M23: Linienschnitt"): the line sections at three heights per band; the same
/// function builds the layers, so check and solver cannot disagree.
void checkOverlaps1D(const CadOptions& o, const Stack& st, Model& m, Report& r) {
  const geometry::LineSectionOptions lo = lineOptions(o);
  std::vector<int> selfDone;
  std::vector<std::array<int, 2>> pairDone;
  Report dummy;
  for (const auto& band : bands(m, o, dummy))
    for (double f : {0.2, 0.5, 0.8}) {
      const double y = band[0] + f * (band[1] - band[0]);
      geometry::LineSectionFindings fd;
      try {
        (void)geometry::lineSection(*m.profiles, y, st.dx_um, m.background, lo, &fd);
      } catch (const std::exception& e) {
        r.errors.push_back({"SHELL_OPEN", std::string("line section failed: ") + e.what(),
                            {kNaN, y, m.zCut}});
        r.ok = false;
        return;
      }
      for (int i : fd.selfOverlap) {
        if (o.periodicUnion || std::find(selfDone.begin(), selfDone.end(), i) != selfDone.end()) continue;
        selfDone.push_back(i);
        const auto& b = m.bodies[static_cast<std::size_t>(i)];
        r.errors.push_back({"PERIODIC_OVERLAP",
                            "body '" + b.name + "' overlaps its periodic image (period " +
                                fmt("%.6g um; smaller body or cad.periodic_union)", st.dx_um),
                            {0.5 * (b.box[0] + b.box[3]) + m.origin[0], y, m.zCut}});
        r.ok = false;
      }
      for (const auto& ij : fd.bodyOverlap) {
        if (std::find(pairDone.begin(), pairDone.end(), ij) != pairDone.end()) continue;
        pairDone.push_back(ij);
        const auto& bi = m.bodies[static_cast<std::size_t>(ij[0])];
        const auto& bj = m.bodies[static_cast<std::size_t>(ij[1])];
        if (lo.priority) {
          r.warnings.push_back({"BODY_OVERLAP", "body '" + bi.name + "' overlaps '" + bj.name +
                                                    "' (priority: the overlap belongs to '" + bj.name + "')"});
        } else {
          r.errors.push_back({"BODY_OVERLAP",
                              "bodies '" + bi.name + "' and '" + bj.name +
                                  "' overlap (cad.overlap = \"priority\" resolves it)",
                              {0.5 * (bi.box[0] + bi.box[3]) + m.origin[0], y, m.zCut}});
          r.ok = false;
        }
      }
    }
}

/// 1D bands -> homogeneous / lamellar1d (extrusion in y, exact) / Sectioned1DLayer (Plan section "M23:
/// Bandklassifikation in 1D"); length fractions in the report.
std::vector<LayerEntry> layersOf1D(const Model& m, const CadOptions& o, const Stack& st, Report& r) {
  const double dx = st.dx_um;
  const int G = std::max(o.normalField.grid, 64);
  const geometry::LineSectionOptions lo = lineOptions(o);
  std::vector<std::string> materials{m.background};
  for (const auto& b : m.bodies)
    if (std::find(materials.begin(), materials.end(), b.material) == materials.end())
      materials.push_back(b.material);
  std::vector<LayerEntry> out;
  for (const auto& band : bands(m, o, r)) {
    const double y0 = band[0], y1 = band[1], ym = 0.5 * (y0 + y1);
    const LineSection mid = geometry::lineSection(*m.profiles, ym, dx, m.background, lo);
    std::map<std::string, double> len;
    for (std::size_t k = 0; k < mid.materials.size(); ++k) {
      const double l = mid.breaks[k + 1] - mid.breaks[k];
      len[mid.materials[k]] += l / dx;
      if (mid.materials.size() > 1 && l < 2 * dx / G)
        r.warnings.push_back(
            {"FEATURE_SMALL",
             fmt("interval of %.3g um at y = %.6g um is smaller than two normal-field raster spacings", l, ym),
             {0.5 * (mid.breaks[k] + mid.breaks[k + 1]) + m.origin[0], ym, m.zCut}});
    }
    Band rb{y0, y1, "", "", {}};
    for (const auto& [mat, f] : len) rb.areaFractions.push_back({mat, f});
    bool extrusion = true;
    for (const auto& b : m.bodies) {
      if (y1 <= b.box[1] || y0 >= b.box[4]) continue;
      extrusion = extrusion && nf::is_extrusion_band(b.shell, {0, 1, 0}, y0, y1);
    }
    LayerEntry le;
    le.factorization = o.factorization;
    if (mid.materials.size() == 1) {
      rb.kind = "homogeneous";
      rb.material = mid.materials.front();
      le.spec = HomogeneousLayer{y1 - y0, rb.material};
    } else if (extrusion) {
      rb.kind = "lamellar1d";  // exact (Li rule); the FFF variants coincide with it
      le.factorization = o.factorization == Factorization::Laurent ? Factorization::Laurent : Factorization::Auto;
      le.spec = Lamellar1DLayer{y1 - y0, mid.breaks, mid.materials};
    } else {
      rb.kind = "sectioned1d";
      Sectioned1DLayer s;
      s.provider = std::make_shared<ProfileLineProvider>(m.profiles, y0, y1 - y0, dx, m.background, lo);
      s.background = m.background;
      s.materials = materials;
      s.slices = o.rcwaSlices;
      s.chebyshevNodes = o.yInterpolation == "chebyshev" ? o.chebyshevNodes : 0;
      le.spec = s;
    }
    r.bands.push_back(rb);
    out.push_back(std::move(le));
  }
  return out;
}

/// Bands -> layers (and the report entries).
// ------------------------------------------------------------------ solids of revolution (M20)

/// Rings of a section of a solid of revolution (Plan "M20: Schichttypen", cad): every loop must be a circle
/// about the axis (x = z = 0) within 1e-9 R_s, every region a disc or an annulus; radii and materials from the
/// axis outwards, the background outside. False with a NOT_AXISYMMETRIC finding otherwise.
bool ringsOfSection(const CrossSection& cs, double Rs, std::array<double, 2> origin, const std::string& background,
                    double y, RingsLayer& out, Report& r) {
  const double tol = 1e-9 * Rs;
  auto fail = [&](const std::string& msg, std::array<double, 2> p) {
    r.errors.push_back({"NOT_AXISYMMETRIC", msg, {p[0] + origin[0], y, p[1] + origin[1]}});
    r.ok = false;
    return false;
  };
  auto circle = [&](const geometry::Contour2D& c, double& radius) {
    const auto pts = c.flatten(1e-3 * Rs);
    double lo = 1e300, hi = 0;
    for (const auto& p : pts) {
      lo = std::min(lo, std::hypot(p[0], p[1]));
      hi = std::max(hi, std::hypot(p[0], p[1]));
    }
    radius = 0.5 * (lo + hi);
    if (hi - lo <= 2 * tol) return true;
    geometry::Point2 worst = pts.front();
    for (const auto& p : pts)
      if (std::abs(std::hypot(p[0], p[1]) - radius) > std::abs(std::hypot(worst[0], worst[1]) - radius)) worst = p;
    return fail(fmt("section at y = %.9g um is not a circle about the axis: the radius varies over %.6g ... %.6g um "
                    "(tolerance 1e-9 R_s; the axis is cad.cell_origin_um, default the centre of the bounding box)",
                    y, lo, hi),
                worst);
  };
  struct Interval {
    double a, b;
    std::string material;
  };
  std::vector<Interval> iv;
  for (const auto& reg : cs.regions) {
    double ro = 0, ri = 0;
    if (!circle(reg.outer, ro)) return false;
    if (reg.holes.size() > 1)
      return fail(fmt("section at y = %.9g um: a region with %.0f holes cannot be rings about the axis", y,
                      static_cast<double>(reg.holes.size())),
                  reg.outer.flatten(1e-3 * Rs).front());
    if (reg.holes.size() == 1 && !circle(reg.holes.front(), ri)) return false;
    iv.push_back({ri, ro, reg.material});
  }
  std::sort(iv.begin(), iv.end(), [](const Interval& p, const Interval& q) { return p.a < q.a; });
  out.radii_um.clear();
  out.materials.clear();
  auto push = [&](double rEnd, const std::string& mat) {  // material mat out to rEnd
    if (!out.materials.empty() && out.materials.back() == mat) {
      out.radii_um.back() = rEnd;
    } else {
      out.radii_um.push_back(rEnd);
      out.materials.push_back(mat);
    }
  };
  double cur = 0;
  for (const auto& v : iv) {
    if (v.a < cur - tol) return fail(fmt("section at y = %.9g um: rings overlap at r = %.9g um", y, v.a), {v.a, 0});
    if (v.a > cur + tol) push(v.a, background);
    push(v.b, v.material);
    cur = v.b;
  }
  if (!out.materials.empty() && out.materials.back() == background) {
    out.radii_um.pop_back();
    out.materials.pop_back();
  }
  out.materials.push_back(background);
  return true;
}

/// Rings of the bodies of one band at any height (y relative to the band bottom, cached per height): the
/// radius functions of a rings_profile band. The rings may not change their materials inside a band.
class AxisymSections {
 public:
  AxisymSections(std::shared_ptr<const std::vector<Body>> b, const Model& m, double y0, double y1, double tol,
                 std::vector<std::string> materials)
      : bodies_(std::move(b)),
        Rs_(m.Rs),
        origin_(m.origin),
        background_(m.background),
        y0_(y0),
        h_(y1 - y0),
        tol_(tol),
        margin_(1e-8 * std::max(m.Rs, y1 - y0)),
        materials_(std::move(materials)) {}
  double radius(double y, std::size_t i) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (y != lastY_) {
      const double yy = y0_ + std::clamp(y, margin_, h_ - margin_);
      Report r;
      RingsLayer rl;
      if (!ringsOfSection(sectionAt(*bodies_, yy, tol_), Rs_, origin_, background_, yy, rl, r))
        throw Error("cad: " + r.errors.front().message);
      FG_REQUIRE(rl.materials == materials_, fmt("cad (axisymmetric): the rings change inside the band at y = %.9g um", yy));
      last_ = std::move(rl.radii_um);
      lastY_ = y;
    }
    return last_.at(i);
  }

 private:
  std::shared_ptr<const std::vector<Body>> bodies_;
  double Rs_;
  std::array<double, 2> origin_;
  std::string background_;
  double y0_, h_, tol_, margin_;
  std::vector<std::string> materials_;
  mutable std::mutex mutex_;
  mutable double lastY_ = kNaN;
  mutable std::vector<double> last_;
};

/// Bands of a solid of revolution: homogeneous, rings (extrusion bands) or rings_profile (radii from the
/// sections at the slice centres, slicing.rcwa_slices slices).
std::vector<LayerEntry> layersOfAxisym(const Model& m, const CadOptions& o, Report& r) {
  const double sampleTol = 1e-6 * m.Rs;
  const auto shared = std::make_shared<const std::vector<Body>>(m.bodies);
  std::vector<LayerEntry> out;
  for (const auto& band : bands(m, o, r)) {
    const double y0 = band[0], y1 = band[1], ym = 0.5 * (y0 + y1);
    const CrossSection mid = sectionAt(m.bodies, ym, sampleTol);
    Band rb{y0, y1, "", "", {}};
    LayerEntry le;
    le.factorization = o.factorization;
    bool any = false, extrusion = true;
    for (const auto& b : m.bodies) {
      if (y1 <= b.box[1] || y0 >= b.box[4]) continue;
      any = true;
      extrusion = extrusion && nf::is_extrusion_band(b.shell, {0, 1, 0}, y0, y1);
    }
    RingsLayer mr;
    if (mid.regions.empty()) {
      std::string mat = m.background;  // no boundary in the band: the material on the axis
      for (const auto& b : m.bodies)
        if (nf::point_in_shell(b.shell, {0, ym, 0})) mat = b.material;
      rb.kind = "homogeneous";
      rb.material = mat;
      le.spec = HomogeneousLayer{y1 - y0, mat};
    } else if (!ringsOfSection(mid, m.Rs, m.origin, m.background, ym, mr, r)) {
      continue;
    } else if (!any || extrusion) {
      rb.kind = "rings";
      mr.thickness_um = y1 - y0;
      le.spec = mr;
    } else {
      bool ok = true;  // the check covers the band (topology is constant between critical heights)
      for (const double f : {0.2, 0.8}) {
        const double yy = y0 + f * (y1 - y0);
        RingsLayer t;
        ok = ok && ringsOfSection(sectionAt(m.bodies, yy, sampleTol), m.Rs, m.origin, m.background, yy, t, r);
      }
      if (!ok) continue;
      rb.kind = "rings_profile";
      RingsProfileLayer rp;
      rp.thickness_um = y1 - y0;
      rp.materials = mr.materials;
      rp.slices = o.rcwaSlices;
      const auto src = std::make_shared<AxisymSections>(shared, m, y0, y1, sampleTol, mr.materials);
      for (std::size_t i = 0; i < mr.radii_um.size(); ++i) {
        RingsProfileLayer::Radius R;
        R.shape = RingsProfileLayer::Radius::Shape::Function;
        R.f = [src, i](double y) { return src->radius(y, i); };
        rp.radii.push_back(std::move(R));
      }
      le.spec = std::move(rp);
    }
    r.bands.push_back(rb);
    out.push_back(std::move(le));
  }
  return out;
}

// ------------------------------------------------------------------ P-fold rotational symmetry (M21)

/// Bands of a P-fold rotationally symmetric arrangement (Plan sec:cylP "CAD"): homogeneous, contour2d (y-invariant
/// bands: the section at mid height) or contour2d_staircase (slicing.rcwa_slices contour2d layers, the sections at
/// the slice centres). Every section must go over into itself under a rotation by 2 pi / P within 1e-9 of its
/// largest radius, as the cylindrical layers check it (E_ROT_SYMMETRY otherwise); r.rotationOrder is the largest
/// order n <= max(P, 12) of all sections at the band mid heights.
std::vector<LayerEntry> layersOfRotation(const Model& m, const CadOptions& o, Report& r) {
  const double sampleTol = 1e-6 * m.Rs;
  const int P = m.rot, maxOrder = std::max(P, 12);
  std::vector<bool> order(static_cast<std::size_t>(maxOrder + 1), true);
  auto tolOf = [](const CrossSection& cs) {
    double rmax = 0;
    for (const auto& reg : cs.regions) rmax = std::max(rmax, cyl::criticalRadii(reg).back());
    return 1e-9 * rmax;
  };
  auto symmetric = [&](const CrossSection& cs, double y) {
    const double tol = tolOf(cs);
    if (cyl::hasRotationSymmetry(cs.regions, P, tol)) return true;
    r.errors.push_back(
        {"E_ROT_SYMMETRY",
         "section at y = " + fmt("%.9g um is not %.0f-fold rotationally symmetric about the axis (largest order found: ", y,
                                 P) +
             std::to_string(cyl::rotationOrder(cs.regions, maxOrder, tol)) +
             "; tolerance 1e-9 of the largest radius; the axis is cad.cell_origin_um, default the centroid of a "
             "section)",
         {m.origin[0], y, m.origin[1]}});
    r.ok = false;
    return false;
  };
  auto contour = [&](const CrossSection& cs, double h) {
    Contour2DLayer c;
    c.thickness_um = h;
    c.background = m.background;
    c.regions = cs.regions;
    c.normalField = o.normalField;
    LayerEntry le;
    le.factorization = o.factorization;
    le.spec = std::move(c);
    return le;
  };
  std::vector<LayerEntry> out;
  for (const auto& band : bands(m, o, r)) {
    const double y0 = band[0], y1 = band[1], ym = 0.5 * (y0 + y1);
    const CrossSection mid = sectionAt(m.bodies, ym, sampleTol);
    Band rb{y0, y1, "", "", {}};
    bool any = false, extrusion = true;
    for (const auto& b : m.bodies) {
      if (y1 <= b.box[1] || y0 >= b.box[4]) continue;
      any = true;
      extrusion = extrusion && nf::is_extrusion_band(b.shell, {0, 1, 0}, y0, y1);
    }
    if (mid.regions.empty()) {
      std::string mat = m.background;  // no boundary in the band: the material on the axis
      for (const auto& b : m.bodies)
        if (nf::point_in_shell(b.shell, {0, ym, 0})) mat = b.material;
      rb.kind = "homogeneous";
      rb.material = mat;
      LayerEntry le;
      le.factorization = o.factorization;
      le.spec = HomogeneousLayer{y1 - y0, mat};
      r.bands.push_back(rb);
      out.push_back(std::move(le));
      continue;
    }
    const double tol = tolOf(mid);
    for (int n = 2; n <= maxOrder; ++n)
      if (order[static_cast<std::size_t>(n)]) order[static_cast<std::size_t>(n)] = cyl::hasRotationSymmetry(mid.regions, n, tol);
    if (!symmetric(mid, ym)) continue;
    if (!any || extrusion) {
      rb.kind = "contour2d";
      out.push_back(contour(mid, y1 - y0));
    } else {
      rb.kind = "contour2d_staircase";
      const int K = std::max(1, o.rcwaSlices);
      bool ok = true;
      std::vector<LayerEntry> steps;
      for (int k = 0; k < K && ok; ++k) {
        const double yy = y0 + (k + 0.5) * (y1 - y0) / K;
        const CrossSection cs = sectionAt(m.bodies, yy, sampleTol);
        ok = symmetric(cs, yy);
        steps.push_back(contour(cs, (y1 - y0) / K));
      }
      if (!ok) continue;
      for (auto& e : steps) out.push_back(std::move(e));
    }
    r.bands.push_back(rb);
  }
  r.rotationOrder = 1;
  for (int n = maxOrder; n > 1; --n)
    if (order[static_cast<std::size_t>(n)]) {
      r.rotationOrder = n;
      break;
    }
  return out;
}

std::vector<LayerEntry> layersOf(const Model& m, const CadOptions& o, const Stack& st, Report& r) {
  if (m.oneD) return layersOf1D(m, o, st, r);
  if (m.axisym) return layersOfAxisym(m, o, r);
  if (m.rot > 0) return layersOfRotation(m, o, r);
  const int G = std::max(o.normalField.grid, 64);
  const double h = std::min(st.dx_um, st.dz_um) / G, sampleTol = 1e-3 * h;
  const auto shared = std::make_shared<const std::vector<Body>>(m.bodies);
  std::vector<LayerEntry> out;
  for (const auto& band : bands(m, o, r)) {
    const double y0 = band[0], y1 = band[1], ym = 0.5 * (y0 + y1);
    const CrossSection mid = sectionAt(m.bodies, ym, sampleTol);
    Band rb{y0, y1, "", "", fractions(mid, m.background, st.dx_um * st.dz_um)};
    for (const auto& reg : mid.regions) {
      const auto b = reg.outer.bounds();
      if (std::min(b[2] - b[0], b[3] - b[1]) < 2 * std::max(st.dx_um, st.dz_um) / G)
        r.warnings.push_back(
            {"FEATURE_SMALL",
             fmt("feature of %.3g um at y = %.6g um is smaller than two normal-field raster "
                 "spacings",
                 std::min(b[2] - b[0], b[3] - b[1]), ym),
             {0.5 * (b[0] + b[2]) + m.origin[0], ym, 0.5 * (b[1] + b[3]) + m.origin[1]}});
    }
    LayerEntry le;
    le.factorization = o.factorization;
    bool any = false, extrusion = true;
    for (const auto& b : m.bodies) {
      if (y1 <= b.box[1] || y0 >= b.box[4]) continue;
      any = true;
      extrusion = extrusion && nf::is_extrusion_band(b.shell, {0, 1, 0}, y0, y1);
    }
    if (mid.regions.empty()) {
      // no boundary in the band: the material at the cell centre
      std::string mat = m.background;
      for (const auto& b : m.bodies)
        if (nf::point_in_shell(b.shell, {0.5 * st.dx_um, ym, 0.5 * st.dz_um})) mat = b.material;
      rb.kind = "homogeneous";
      rb.material = mat;
      le.spec = HomogeneousLayer{y1 - y0, mat};
    } else if (!any || extrusion) {
      rb.kind = "extrusion";
      Contour2DLayer c;
      c.thickness_um = y1 - y0;
      c.background = m.background;
      c.regions = mid.regions;
      c.normalField = o.normalField;
      le.spec = c;
    } else {
      rb.kind = "sectioned";
      SectionedLayer s;
      s.provider = std::make_shared<NurbforgeSectionProvider>(shared, y0, y1, sampleTol,
                                                              std::max(st.dx_um, st.dz_um));
      s.background = m.background;
      s.slices = o.rcwaSlices;
      s.normalField = o.normalField;
      s.chebyshevNodes = o.yInterpolation == "chebyshev" ? o.chebyshevNodes : 0;
      le.spec = s;
    }
    r.bands.push_back(rb);
    out.push_back(std::move(le));
  }
  return out;
}

void checkOverlaps(const CadOptions& o, const Stack& st, Model& m, Report& r) {
  if (m.oneD) {
    checkOverlaps1D(o, st, m, r);
    return;
  }
  // axisymmetric and rotational symmetry: no periodic images
  const bool axial = m.axisym || m.rot > 0;
  const double tol = 1e-4 * (axial ? m.Rs : std::min(st.dx_um, st.dz_um));
  const int span = axial ? 0 : 1;
  std::vector<std::array<int, 3>> reported;  // (i, j, periodic)
  // priority changes one pair per round, then everything is checked again; one round more than pairs
  // verifies the last change (v2.24; v2.21 stopped after four rounds without a final check)
  const std::size_t n = m.bodies.size(), rounds = n * (n - 1) / 2 + 2;
  bool changed = true;
  for (std::size_t round = 0; round < rounds && changed; ++round) {
    changed = false;
    Report dummy;
    for (const auto& band : bands(m, o, dummy)) {
      for (double f : {0.2, 0.5, 0.8}) {
        const double y = band[0] + f * (band[1] - band[0]);
        std::vector<Flat> fl;
        for (const auto& b : m.bodies) fl.push_back(flatSection(b, y, tol));
        for (std::size_t i = 0; i < m.bodies.size(); ++i)
          for (std::size_t j = 0; j <= i; ++j)
            for (int di = -span; di <= span; ++di)
              for (int dk = -span; dk <= span; ++dk) {
                const bool periodic = di != 0 || dk != 0;
                if (i == j && !periodic) continue;
                const std::array<int, 3> key{static_cast<int>(i), static_cast<int>(j),
                                             periodic ? 1 : 0};
                if (i == j && touches(fl[i], fl[j], di * st.dx_um, dk * st.dz_um, tol)) {
                  // not resolvable by periodic_union (the union would grow around the cell)
                  if (std::find(reported.begin(), reported.end(), key) != reported.end()) continue;
                  reported.push_back(key);
                  const auto& bx = fl[i].box;
                  r.errors.push_back(
                      {"PERIODIC_OVERLAP",
                       "body '" + m.bodies[i].name + "' touches its periodic image (shift " +
                           fmt("%.6g, %.6g um", di * st.dx_um, dk * st.dz_um) +
                           "): it spans the period, e.g. a z-invariant line extruded over the "
                           "whole period; such sections wrap around the cell and are not closed "
                           "contours (describe z-invariant gratings as 1D: lamellar1d / profile1d)",
                       {0.5 * (bx[0] + bx[2]) + m.origin[0], y,
                        0.5 * (bx[1] + bx[3]) + m.origin[1]}});
                  r.ok = false;
                  continue;
                }
                if (!overlap(fl[i], fl[j], di * st.dx_um, dk * st.dz_um)) continue;
                if (o.overlap == "priority" && i != j && !periodic) {
                  try {
                    nf::Shell d = nf::boolean(m.bodies[i].shell, m.bodies[j].shell, nf::BooleanOp::Difference);
                    if (!d.free_edges().empty()) throw Error("the result is not closed");
                    m.bodies[i].shell = std::move(d);
                  } catch (const std::exception& e) {  // e.g. an open result shell
                    r.errors.push_back({"BODY_OVERLAP",
                                        "bodies '" + m.bodies[i].name + "' and '" + m.bodies[j].name +
                                            "' overlap and the Boolean difference of cad.overlap = \"priority\" "
                                            "failed (" + e.what() + "); remove the overlap in the model"});
                    r.ok = false;
                    return;
                  }
                  m.bodies[i].box = bounds(m.bodies[i].shell);
                  r.warnings.push_back(
                      {"BODY_OVERLAP",
                       "body '" + m.bodies[i].name + "' overlaps '" + m.bodies[j].name +
                           "' (priority: the overlap belongs to '" + m.bodies[j].name + "')"});
                  changed = true;
                  break;
                }
                if (std::find(reported.begin(), reported.end(), key) != reported.end()) continue;
                reported.push_back(key);
                const auto& bx = fl[i].box;
                r.errors.push_back(
                    {periodic ? "PERIODIC_OVERLAP" : "BODY_OVERLAP",
                     periodic && i == j
                         ? "body '" + m.bodies[i].name + "' overlaps its periodic image (period " +
                               fmt("%.6g x %.6g um", st.dx_um, st.dz_um) +
                               "; make the body smaller; cad.periodic_union cannot resolve this in 2D: the "
                               "united body would wrap around the cell)"
                         : "bodies '" + m.bodies[i].name + "' and '" + m.bodies[j].name + "'" +
                               (periodic ? " (periodic image)" : "") +
                               " overlap (cad.overlap = \"priority\" resolves it)",
                     {0.5 * (bx[0] + bx[2]) + m.origin[0], y,
                      0.5 * (bx[1] + bx[3]) + m.origin[1]}});
                r.ok = false;
              }
        if (changed) break;
      }
      if (changed) break;
    }
  }
  if (changed) {
    r.errors.push_back({"BODY_OVERLAP", "cad.overlap = \"priority\": overlaps not resolved after " +
                                            std::to_string(rounds) + " Boolean differences"});
    r.ok = false;
  }
}

}  // namespace

bool available() {
  return true;
}

Report check(const CadOptions& o, const Stack& st) {
  Report r;
  Model m;
  r.axisymmetric = o.axisymmetric;
  r.rotation = o.rotation;
  r.oneD = !o.axisymmetric && o.rotation == 0 && st.dz_um <= 0;
  if (!loadModel(o, st, m, r)) return r;
  checkModel(o, m, r);
  if (r.ok) centreOnAxis(o, m);
  r.cellOrigin_um = m.origin;
  r.zCut_um = m.zCut;
  if (r.ok) checkOverlaps(o, st, m, r);
  if (!r.ok) return r;
  try {
    (void)layersOf(m, o, st, r);
  } catch (const std::exception& e) {
    r.errors.push_back({"SHELL_OPEN", std::string("slicing failed: ") + e.what()});
    r.ok = false;
  }
  return r;
}

std::vector<LayerEntry> slice(const CadOptions& o, const Stack& st, Report* report) {
  Report r;
  Model m;
  std::vector<LayerEntry> layers;
  r.axisymmetric = o.axisymmetric;
  r.rotation = o.rotation;
  r.oneD = !o.axisymmetric && o.rotation == 0 && st.dz_um <= 0;
  if (loadModel(o, st, m, r)) {
    checkModel(o, m, r);
    if (r.ok) centreOnAxis(o, m);
    r.cellOrigin_um = m.origin;
    r.zCut_um = m.zCut;
    if (r.ok) checkOverlaps(o, st, m, r);
    if (r.ok) layers = layersOf(m, o, st, r);
  }
  if (report) *report = r;
  if (!r.ok) {
    std::string msg = "cad: geometry check failed:";
    for (const auto& e : r.errors) msg += "\n  " + e.code + ": " + e.message;
    throw Error(msg);
  }
  return layers;
}

/// 1D: the profiles in the (x, y) plane over [0, dx] with the images at +- dx and the band boundaries.
std::string previewSvg1D(const CadOptions& o, const Stack& st, Model& m, Report& r) {
  checkModel(o, m, r);
  if (!m.profiles)
    throw Error("cad preview: " + (r.errors.empty() ? std::string() : r.errors.front().message));
  static const char* palette[] = {"#3b6fb6", "#d1495b", "#edae49", "#66a182",
                                  "#8d6cab", "#00798c", "#c97b84"};
  std::map<std::string, std::string> colour;
  for (const auto& p : *m.profiles)
    if (!colour.count(p.material)) colour[p.material] = palette[colour.size() % 7];
  const double dx = st.dx_um, ya = m.yRange[0], yb = m.yRange[1], H = yb - ya;
  const double W = 400, sc = W / std::max(dx, H);
  std::ostringstream svg;
  // y up: the drawing uses (x, -y)
  svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << dx * sc << "\" height=\"" << H * sc
      << "\" viewBox=\"0 " << -yb << " " << dx << " " << H << "\">\n";
  svg << "<title>profile at z = " << m.zCut << " um</title>\n";
  svg << "<rect x=\"0\" y=\"" << -yb << "\" width=\"" << dx << "\" height=\"" << H
      << "\" fill=\"#f2f2f2\" stroke=\"#888\" stroke-width=\"" << 0.003 * dx << "\"/>\n";
  for (int di = -1; di <= 1; ++di)
    for (const auto& p : *m.profiles) {
      svg << "<path fill-rule=\"evenodd\" fill=\"" << colour[p.material] << "\" transform=\"translate("
          << di * dx << " 0)\" d=\"";
      for (const auto& c : p.loops) {
        const auto q = c.flatten(1e-3 * std::max(dx, H));
        for (std::size_t k = 0; k < q.size(); ++k) svg << (k ? " L" : "M") << q[k][0] << " " << -q[k][1];
        svg << " Z ";
      }
      svg << "\"><title>" << p.name << " (" << p.material << ")</title></path>\n";
    }
  Report dummy;
  for (const auto& band : bands(m, o, dummy))
    svg << "<line x1=\"0\" x2=\"" << dx << "\" y1=\"" << -band[1] << "\" y2=\"" << -band[1]
        << "\" stroke=\"#555\" stroke-dasharray=\"" << 0.01 * dx << "\" stroke-width=\"" << 0.002 * dx
        << "\"/>\n";
  svg << "</svg>\n";
  return svg.str();
}

std::string previewSvg(const CadOptions& o, const Stack& st, double y) {
  Report r;
  Model m;
  if (!loadModel(o, st, m, r))
    throw Error("cad preview: " + (r.errors.empty() ? std::string() : r.errors.front().message));
  if (m.oneD) return previewSvg1D(o, st, m, r);
  FG_REQUIRE(!m.axisym && m.rot == 0,
             "cad preview: not available in a cylindrical stack (M20, M21); see the report bands");
  FG_REQUIRE(std::isfinite(y), "cad preview: y_um is required for 2D gratings");
  const double h = std::min(st.dx_um, st.dz_um) / 512;
  const CrossSection cs = sectionAt(m.bodies, y, 1e-3 * h);
  // colours per material (stable order)
  static const char* palette[] = {"#3b6fb6", "#d1495b", "#edae49", "#66a182",
                                  "#8d6cab", "#00798c", "#c97b84"};
  std::map<std::string, std::string> colour;
  for (const auto& reg : cs.regions)
    if (!colour.count(reg.material)) colour[reg.material] = palette[colour.size() % 7];
  const double W = 400, sc = W / std::max(st.dx_um, st.dz_um);
  std::ostringstream svg;
  svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << st.dx_um * sc << "\" height=\""
      << st.dz_um * sc << "\" viewBox=\"0 0 " << st.dx_um << " " << st.dz_um << "\">\n";
  svg << "<title>section y = " << y << " um</title>\n";
  svg << "<rect x=\"0\" y=\"0\" width=\"" << st.dx_um << "\" height=\"" << st.dz_um
      << "\" fill=\"#f2f2f2\" stroke=\"#888\" stroke-width=\"" << 0.003 * st.dx_um << "\"/>\n";
  for (int di = -1; di <= 1; ++di)
    for (int dk = -1; dk <= 1; ++dk)
      for (const auto& reg : cs.regions) {
        svg << "<path fill-rule=\"evenodd\" fill=\"" << colour[reg.material]
            << "\" transform=\"translate(" << di * st.dx_um << " " << dk * st.dz_um << ")\" d=\"";
        auto contour = [&](const geometry::Contour2D& c) {
          const auto p = c.flatten(1e-3 * std::max(st.dx_um, st.dz_um));
          for (std::size_t k = 0; k < p.size(); ++k)
            svg << (k ? " L" : "M") << p[k][0] << " " << p[k][1];
          svg << " Z ";
        };
        contour(reg.outer);
        for (const auto& hl : reg.holes) contour(hl);
        svg << "\"><title>" << reg.material << "</title></path>\n";
      }
  svg << "</svg>\n";
  return svg.str();
}

}  // namespace fg::cad
