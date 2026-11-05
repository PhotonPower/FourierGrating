// FourierGrating - CAD slicer with Nurbforge (Plan section "M14: CAD-Slicer"): STEP -> FG
// coordinates (micrometres, axes) -> checks -> bands between critical heights -> homogeneous /
// Contour2D / SectionedLayer. The only module that includes Nurbforge headers.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <nurbforge/analysis/shell_section.hpp>
#include <nurbforge/core/error.hpp>
#include <nurbforge/geometry/transform.hpp>
#include <nurbforge/io/step.hpp>
#include <nurbforge/topology/boolean.hpp>
#include <sstream>

#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
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

std::array<double, 6> bounds(const nf::Shell& s) {
  std::array<double, 6> b{1e300, 1e300, 1e300, -1e300, -1e300, -1e300};
  for (const auto& f : s.faces())
    for (const auto& p : f.surface().control_points()) {
      b[0] = std::min(b[0], p.x);
      b[1] = std::min(b[1], p.y);
      b[2] = std::min(b[2], p.z);
      b[3] = std::max(b[3], p.x);
      b[4] = std::max(b[4], p.y);
      b[5] = std::max(b[5], p.z);
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
};

bool loadModel(const CadOptions& o, const Stack& st, Model& m, Report& r) {
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
  m.origin = o.cellOrigin_um ? *o.cellOrigin_um : std::array<double, 2>{all[0], all[2]};
  const nf::Transform shift = nf::Transform::translation({-m.origin[0], 0, -m.origin[1]});
  for (auto& b : m.bodies) {
    b.shell = transformed(b.shell, shift);
    b.box = bounds(b.shell);
  }
  m.yRange = o.yRange_um ? *o.yRange_um : std::array<double, 2>{all[1], all[4]};
  FG_REQUIRE(m.yRange[1] > m.yRange[0], "cad.y_range_um: empty range");
  m.background = resolveBackground(o, st);
  FG_REQUIRE(st.materials.count(m.background),
             "cad.background: unknown material '" + m.background + "'");
  // ---- units: dimensions far away from the period / wavelength (factor 1e3)
  const double period = std::max(st.dx_um, st.dz_um),
               extent = std::max({all[3] - all[0], all[4] - all[1], all[5] - all[2]});
  const double ref = o.lambda_um > 0 ? std::max(period, o.lambda_um) : period;
  if (extent > 100 * ref || extent < 1e-3 * std::min(period, ref))
    error("UNIT_SUSPECT", fmt("bodies are %.6g um large for a period of %.6g um (factor 1000 off? "
                              "cad.unit \"mm\" / \"um\", or "
                              "write the STEP file with unit micrometre)",
                              extent, period));
  return true;
}

// ------------------------------------------------------------------ checks

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

/// Do the sections a and b + (sx, sz) overlap (boundaries cross or one contains the other)?
bool overlap(const Flat& a, const Flat& b, double sx, double sz) {
  if (a.loops.empty() || b.loops.empty()) return false;
  if (b.box[0] + sx >= a.box[2] || a.box[0] >= b.box[2] + sx || b.box[1] + sz >= a.box[3] ||
      a.box[1] >= b.box[3] + sz)
    return false;
  const auto o = [](const geometry::Point2& p, const geometry::Point2& q,
                    const geometry::Point2& r) {
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
  };
  for (const auto& la : a.loops)
    for (const auto& lb : b.loops)
      for (std::size_t i = 0; i < la.size(); ++i) {
        const auto& p = la[i];
        const auto& q = la[(i + 1) % la.size()];
        for (std::size_t j = 0; j < lb.size(); ++j) {
          const geometry::Point2 r{lb[j][0] + sx, lb[j][1] + sz},
              t{lb[(j + 1) % lb.size()][0] + sx, lb[(j + 1) % lb.size()][1] + sz};
          if (o(p, q, r) * o(p, q, t) < 0 && o(r, t, p) * o(r, t, q) < 0) return true;
        }
      }
  for (const auto& ra : a.regions)
    if (ra.contains(b.loops[0][0][0] + sx, b.loops[0][0][1] + sz)) return true;
  for (const auto& rb : b.regions)
    if (rb.contains(a.loops[0][0][0] - sx, a.loops[0][0][1] - sz)) return true;
  return false;
}

/// BODY_OVERLAP / PERIODIC_OVERLAP from sections at three heights per band (the topology of the
/// sections does not change inside a band); "priority" and periodic_union resolve them with
/// Booleans.
void checkOverlaps(const CadOptions& o, const Stack& st, Model& m, Report& r);

/// Bands -> layers (and the report entries).
std::vector<LayerEntry> layersOf(const Model& m, const CadOptions& o, const Stack& st, Report& r) {
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
  const double tol = 1e-4 * std::min(st.dx_um, st.dz_um);
  std::vector<std::array<int, 3>> reported;  // (i, j, periodic)
  bool changed = true;
  for (int round = 0; round < 4 && changed;
       ++round) {  // priority / union change the bodies: check again
    changed = false;
    Report dummy;
    for (const auto& band : bands(m, o, dummy)) {
      for (double f : {0.2, 0.5, 0.8}) {
        const double y = band[0] + f * (band[1] - band[0]);
        std::vector<Flat> fl;
        for (const auto& b : m.bodies) fl.push_back(flatSection(b, y, tol));
        for (std::size_t i = 0; i < m.bodies.size(); ++i)
          for (std::size_t j = 0; j <= i; ++j)
            for (int di = -1; di <= 1; ++di)
              for (int dk = -1; dk <= 1; ++dk) {
                const bool periodic = di != 0 || dk != 0;
                if (i == j && !periodic) continue;
                if (!overlap(fl[i], fl[j], di * st.dx_um, dk * st.dz_um)) continue;
                const std::array<int, 3> key{static_cast<int>(i), static_cast<int>(j),
                                             periodic ? 1 : 0};
                if (o.overlap == "priority" && i != j && !periodic) {
                  m.bodies[i].shell =
                      nf::boolean(m.bodies[i].shell, m.bodies[j].shell, nf::BooleanOp::Difference);
                  m.bodies[i].box = bounds(m.bodies[i].shell);
                  r.warnings.push_back(
                      {"BODY_OVERLAP",
                       "body '" + m.bodies[i].name + "' overlaps '" + m.bodies[j].name +
                           "' (priority: the overlap belongs to '" + m.bodies[j].name + "')"});
                  changed = true;
                  break;
                }
                if (o.periodicUnion && i == j) {
                  const nf::Shell img =
                      transformed(m.bodies[i].shell,
                                  nf::Transform::translation({di * st.dx_um, 0, dk * st.dz_um}));
                  m.bodies[i].shell = nf::boolean(m.bodies[i].shell, img, nf::BooleanOp::Union);
                  m.bodies[i].box = bounds(m.bodies[i].shell);
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
                               "; smaller body or cad.periodic_union)"
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
}

}  // namespace

bool available() {
  return true;
}

Report check(const CadOptions& o, const Stack& st) {
  Report r;
  Model m;
  if (!loadModel(o, st, m, r)) return r;
  r.cellOrigin_um = m.origin;
  checkModel(o, m, r);
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
  if (loadModel(o, st, m, r)) {
    r.cellOrigin_um = m.origin;
    checkModel(o, m, r);
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

std::string previewSvg(const CadOptions& o, const Stack& st, double y) {
  Report r;
  Model m;
  if (!loadModel(o, st, m, r))
    throw Error("cad preview: " + (r.errors.empty() ? std::string() : r.errors.front().message));
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
