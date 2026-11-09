#include "fg/io/config.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/physics/anisotropic.hpp"

namespace fg::io {

using nlohmann::json;

namespace {

constexpr double DEG = PI / 180.0;

cplx toComplex(const json& v) {
  if (v.is_number()) return {v.get<double>(), 0.0};
  FG_REQUIRE(v.is_array() && v.size() == 2, "complex numbers are written as number or [re, im]");
  return {v[0].get<double>(), v[1].get<double>()};
}

Material parseMaterial(const std::string& name, const json& j) {
  if (j.contains("n")) return Material::fromIndex(name, toComplex(j["n"]));
  if (j.contains("eps")) {
    Material m = Material::fromIndex(name, 1.0);
    m.eps = toComplex(j["eps"]);
    m.epsTensor = m.eps * Mat3C::Identity();
    return m;
  }
  if (j.contains("eps_tensor")) {
    const json& t = j["eps_tensor"];
    FG_REQUIRE(t.is_array() && t.size() == 3, "eps_tensor must be 3x3");
    Mat3C e;
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) e(i, k) = toComplex(t[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]);
    return Material::fromTensor(name, e);
  }
  if (j.contains("eps_principal")) {
    const json& p = j["eps_principal"];
    FG_REQUIRE(p.is_array() && p.size() == 3, "eps_principal needs three values");
    std::vector<double> eu{0.0, 0.0, 0.0};
    if (j.contains("euler_deg")) eu = j["euler_deg"].get<std::vector<double>>();
    FG_REQUIRE(eu.size() == 3, "euler_deg needs three angles (z-x-z)");
    return Material::fromTensor(name, rotatedTensor(toComplex(p[0]), toComplex(p[1]), toComplex(p[2]),
                                                    eu[0] * DEG, eu[1] * DEG, eu[2] * DEG));
  }
  throw Error("material '" + name + "' needs 'n', 'eps', 'eps_tensor' or 'eps_principal'");
}

Factorization parseFactorization(const json& j) {
  const std::string f = j.value("factorization", std::string("auto"));
  if (f == "auto") return Factorization::Auto;
  if (f == "laurent") return Factorization::Laurent;
  if (f == "li_lamellar") return Factorization::LiLamellar;
  if (f == "FFF_iso" || f == "fff_iso") return Factorization::FFF_iso;
  if (f == "FFF_sym" || f == "fff_sym") return Factorization::FFF_sym;
  if (f == "FFF_rig" || f == "fff_rig") return Factorization::FFF_rig;
  if (f == "FFF_A" || f == "fff_a") return Factorization::FFF_A;
  if (f == "FFF_B" || f == "fff_b") return Factorization::FFF_B;
  throw Error("unknown factorization '" + f + "'");
}

NormalFieldOptions parseNormalField(const json& nf) {
  NormalFieldOptions o;
  const std::string m = nf.value("method", std::string("tensor"));
  FG_REQUIRE(m == "tensor" || m == "vector", "normal_field.method must be \"tensor\" or \"vector\"");
  o.method = m == "vector" ? NormalFieldOptions::Method::Vector : NormalFieldOptions::Method::Tensor;
  o.grid = nf.value("grid", o.grid);
  o.sigma = nf.value("sigma", o.sigma);
  o.tau = nf.value("tau", o.tau);
  return o;
}

geometry::Point2 point2(const json& v, const char* what) {
  FG_REQUIRE(v.is_array() && v.size() == 2, std::string(what) + ": points are [x, z]");
  return {v[0].get<double>(), v[1].get<double>()};
}

/// One closed contour of a contour2d region or hole (Plan section "M13: Neuer Schichttyp Contour2D").
geometry::Contour2D parseContour(const json& j) {
  if (j.contains("circle")) {
    const json& c = j["circle"];
    return geometry::Contour2D::circle(point2(c.at("center_um"), "circle.center_um"), c.at("radius_um").get<double>());
  }
  if (j.contains("ellipse")) {
    const json& e = j["ellipse"];
    const auto r = e.at("radii_um").get<std::vector<double>>();
    FG_REQUIRE(r.size() == 2, "ellipse.radii_um: [a, b]");
    return geometry::Contour2D::ellipse(point2(e.at("center_um"), "ellipse.center_um"), r[0], r[1],
                                        e.value("angle_deg", 0.0) * PI / 180.0);
  }
  if (j.contains("polygon_um")) {
    std::vector<geometry::Point2> v;
    for (const auto& p : j["polygon_um"]) v.push_back(point2(p, "polygon_um"));
    return geometry::Contour2D::polygon(v);
  }
  if (j.contains("nurbs")) {
    const json& n = j["nurbs"];
    std::vector<geometry::Point2> pts;
    for (const auto& p : n.at("points_um")) pts.push_back(point2(p, "nurbs.points_um"));
    return geometry::Contour2D::nurbs(n.at("degree").get<int>(), n.at("knots").get<std::vector<double>>(), pts,
                                      n.value("weights", std::vector<double>{}));
  }
  throw Error("contour2d: region/hole needs 'circle', 'ellipse', 'polygon_um' or 'nurbs'");
}

/// The cad block (Plan section "M14: Konfiguration"); needs periods, materials and the scan of cfg.
cad::CadOptions parseCadOptions(const json& j, const std::string& baseDir, const RunConfig& cfg) {
  cad::CadOptions o;
  o.file = j.value("file", std::string());
  if (!o.file.empty() && std::filesystem::path(o.file).is_relative() && !baseDir.empty())
    o.file = (std::filesystem::path(baseDir) / o.file).string();
  o.stepText = j.value("step", std::string());  // contents instead of a file (Python fg.cad.check)
  FG_REQUIRE(!o.file.empty() || !o.stepText.empty(), "cad: 'file' (STEP) missing");
  o.unit = j.value("unit", o.unit);
  if (j.contains("axes")) {
    o.axisY = j["axes"].value("y", o.axisY);
    o.axisX = j["axes"].value("x", o.axisX);
    o.axisZ = j["axes"].value("z", o.axisZ);
  }
  if (j.contains("y_range_um")) o.yRange_um = std::array<double, 2>{j["y_range_um"][0].get<double>(), j["y_range_um"][1].get<double>()};
  if (j.contains("cell_origin_um"))
    o.cellOrigin_um = std::array<double, 2>{j["cell_origin_um"][0].get<double>(), j["cell_origin_um"][1].get<double>()};
  o.background = j.value("background", o.background);
  for (const auto& b : j.value("bodies", json::array())) {
    cad::BodyAssignment a;
    a.name = b.value("name", std::string());
    a.index = b.value("index", -1);
    a.material = b.at("material").get<std::string>();
    FG_REQUIRE(!a.name.empty() || a.index >= 0, "cad.bodies: every entry needs 'name' or 'index'");
    o.bodies.push_back(std::move(a));
  }
  o.overlap = j.value("overlap", o.overlap);
  FG_REQUIRE(o.overlap == "error" || o.overlap == "priority", "cad.overlap must be \"error\" or \"priority\"");
  o.periodicUnion = j.value("periodic_union", o.periodicUnion);
  if (j.contains("slicing")) {
    const json& s = j["slicing"];
    o.criticalTol_um = s.value("critical_tol_um", o.criticalTol_um);
    o.minBand_um = s.value("min_band_um", o.minBand_um);
    o.rcwaSlices = s.value("rcwa_slices", o.rcwaSlices);
    o.yInterpolation = s.value("y_interpolation", o.yInterpolation);
    FG_REQUIRE(o.yInterpolation == "none" || o.yInterpolation == "chebyshev",
               "cad.slicing.y_interpolation must be \"none\" or \"chebyshev\"");
    o.chebyshevNodes = s.value("chebyshev_nodes", o.chebyshevNodes);
  }
  o.factorization = parseFactorization(j);
  if (j.contains("normal_field")) o.normalField = parseNormalField(j["normal_field"]);
  for (const auto& w : cfg.scan) o.lambda_um = o.lambda_um > 0 ? std::min(o.lambda_um, w.lambda_um) : w.lambda_um;
  return o;
}

LayerEntry parseLayer(const json& j) {
  const std::string type = j.at("type").get<std::string>();
  LayerEntry le;
  le.factorization = parseFactorization(j);
  if (type == "homogeneous") {
    le.spec = HomogeneousLayer{j.at("thickness_um").get<double>(), j.at("material").get<std::string>()};
  } else if (type == "lamellar1d") {
    le.spec = Lamellar1DLayer{j.at("thickness_um").get<double>(),
                              j.at("breakpoints_um").get<std::vector<double>>(),
                              j.at("materials").get<std::vector<std::string>>()};
  } else if (type == "profile1d") {
    Profile1DLayer p;
    p.depth_um = j.at("depth_um").get<double>();
    p.below = j.at("below").get<std::string>();
    p.above = j.at("above").get<std::string>();
    p.slices = j.value("slices", 0);
    const std::string shape = j.value("shape", std::string("sinus"));
    if (shape == "sinus") {
      p.shape = Profile1DLayer::Shape::Sinus;
      p.x0_um = j.value("x0_um", 0.0);
    } else if (shape == "triangle" || shape == "points") {
      p.shape = Profile1DLayer::Shape::PiecewiseLinear;
      if (shape == "triangle") {
        p.xs_um = {0.0, j.at("apex_x_um").get<double>()};
        p.ys_um = {0.0, p.depth_um};
      } else {
        p.xs_um = j.at("x_um").get<std::vector<double>>();
        p.ys_um = j.at("y_um").get<std::vector<double>>();
        FG_REQUIRE(p.xs_um.size() == p.ys_um.size() && p.xs_um.size() >= 2, "profile points: x_um/y_um");
      }
    } else {
      throw Error("unknown profile shape '" + shape + "'");
    }
    le.spec = p;
  } else if (type == "polygon2d") {
    Polygon2DLayer g;
    g.thickness_um = j.at("thickness_um").get<double>();
    g.background = j.at("background").get<std::string>();
    for (const auto& ij : j.value("inclusions", json::array())) {
      Polygon2DLayer::Inclusion inc;
      inc.material = ij.at("material").get<std::string>();
      if (ij.contains("vertices_um")) {
        for (const auto& v : ij["vertices_um"]) {
          FG_REQUIRE(v.is_array() && v.size() == 2, "polygon2d: vertices_um entries are [x, z]");
          inc.vertices.push_back({v[0].get<double>(), v[1].get<double>()});
        }
      } else if (ij.contains("rectangle_um")) {
        const auto c = ij["rectangle_um"].at("center").get<std::vector<double>>();
        const auto w = ij["rectangle_um"].at("size").get<std::vector<double>>();
        FG_REQUIRE(c.size() == 2 && w.size() == 2, "polygon2d: rectangle_um needs center [x, z] and size [wx, wz]");
        const double x0 = c[0] - 0.5 * w[0], x1 = c[0] + 0.5 * w[0], z0 = c[1] - 0.5 * w[1], z1 = c[1] + 0.5 * w[1];
        inc.vertices = {{x0, z0}, {x1, z0}, {x1, z1}, {x0, z1}};
      } else {
        throw Error("polygon2d: inclusion needs 'vertices_um' or 'rectangle_um'");
      }
      FG_REQUIRE(inc.vertices.size() >= 3, "polygon2d: inclusion needs at least 3 vertices");
      g.inclusions.push_back(std::move(inc));
    }
    if (j.contains("normal_field")) g.normalField = parseNormalField(j["normal_field"]);
    le.spec = g;
  } else if (type == "contour2d") {
    Contour2DLayer c;
    c.thickness_um = j.at("thickness_um").get<double>();
    c.background = j.at("background").get<std::string>();
    for (const auto& rj : j.value("regions", json::array())) {
      geometry::Region r;
      r.material = rj.at("material").get<std::string>();
      r.outer = parseContour(rj);
      for (const auto& hj : rj.value("holes", json::array())) r.holes.push_back(parseContour(hj));
      c.regions.push_back(std::move(r));
    }
    if (j.contains("normal_field")) c.normalField = parseNormalField(j["normal_field"]);
    le.spec = c;
  } else if (type == "profile2d") {
    Profile2DLayer p;
    p.depth_um = j.at("depth_um").get<double>();
    p.below = j.at("below").get<std::string>();
    p.above = j.at("above").get<std::string>();
    p.slices = j.value("slices", 0);
    p.zSamples = j.value("z_samples", 0);
    p.x0_um = j.value("x0_um", 0.0);
    p.z0_um = j.value("z0_um", 0.0);
    const std::string shape = j.value("shape", std::string("egg_crate"));
    if (shape == "egg_crate") p.shape = Profile2DLayer::Shape::EggCrate;
    else if (shape == "sinus_x") p.shape = Profile2DLayer::Shape::SinusX;
    else throw Error("unknown profile2d shape '" + shape + "' (egg_crate | sinus_x)");
    le.spec = p;
  } else if (type == "rings") {  // M20 (Plan section "M20: Schichttypen")
    le.spec = RingsLayer{j.at("thickness_um").get<double>(), j.at("radii_um").get<std::vector<double>>(),
                         j.at("materials").get<std::vector<std::string>>()};
  } else {
    throw Error("unknown layer type '" + type + "'");
  }
  return le;
}

std::vector<double> parseRange(const json& v) {
  if (v.is_number()) return {v.get<double>()};
  if (v.is_array()) return v.get<std::vector<double>>();
  const double a = v.at("from").get<double>(), b = v.at("to").get<double>();
  const int n = v.at("steps").get<int>();
  FG_REQUIRE(n >= 1, "range needs steps >= 1");
  std::vector<double> out;
  for (int k = 0; k < n; ++k) out.push_back(n == 1 ? a : a + (b - a) * k / (n - 1));
  return out;
}

}  // namespace

RunConfig parseConfig(const std::string& text, const std::string& baseDir) {
  const json j = json::parse(text);
  RunConfig cfg;
  // rotationally symmetric structures (M20): no period, no truncation, "dimension" optional
  const bool cylindrical = j.contains("cylindrical");
  const std::string dim = cylindrical ? j.value("dimension", std::string("1D")) : j.at("dimension").get<std::string>();
  FG_REQUIRE(dim == "1D" || dim == "2D", "dimension must be \"1D\" or \"2D\"");
  const std::string solver = j.value("solver", std::string("rcwa"));
  FG_REQUIRE(solver == "rcwa" || solver == "rk", "solver must be \"rcwa\" or \"rk\"");
  cfg.solver.kind = solver == "rk" ? SolverKind::RK : SolverKind::RCWA;
  const json trunc = cylindrical ? j.value("truncation", json::object()) : j.at("truncation");
  cfg.solver.Nx = trunc.value("Nx", 10);
  cfg.solver.Nz = trunc.value("Nz", 0);
  if (j.contains("rk")) {
    cfg.solver.rkTableau = j["rk"].value("tableau", cfg.solver.rkTableau);
    cfg.solver.rkTol = j["rk"].value("tol", cfg.solver.rkTol);
  }
  if (j.contains("rcwa")) cfg.solver.rcwaSlices = j["rcwa"].value("slices", cfg.solver.rcwaSlices);
  cfg.solver.backend = j.value("backend", cfg.solver.backend);
  FG_REQUIRE(cfg.solver.backend == "auto" || cfg.solver.backend == "cpu" || cfg.solver.backend == "cuda",
             "backend must be \"auto\", \"cpu\" or \"cuda\"");
  cfg.solver.backendThreshold = j.value("backend_threshold", cfg.solver.backendThreshold);
  cfg.solver.splitPolarizations = j.value("split_polarizations", cfg.solver.splitPolarizations);
  cfg.solver.fftSampling = j.value("fft_sampling", cfg.solver.fftSampling);
  // scattered-field formulation (M18; automatic for aperiodic structures under oblique incidence)
  cfg.solver.scatteredField = j.value("scattered_field", cfg.solver.scatteredField);
  // thread model (Plan section "M15", step 2): {"scan": n | "auto", "slices": n | "auto", "blas": n}
  if (j.contains("threads")) {
    const json& t = j["threads"];
    FG_REQUIRE(t.is_object(), "threads must be an object {\"scan\", \"slices\", \"blas\"}");
    auto count = [&](const char* key, int def, bool autoAllowed) {
      if (!t.contains(key)) return def;
      const json& v = t[key];
      if (autoAllowed && v.is_string() && v.get<std::string>() == "auto") return 0;
      FG_REQUIRE(v.is_number_integer() && v.get<int>() >= 1,
                 std::string("threads.") + key + (autoAllowed ? " must be a positive integer or \"auto\"" : " must be a positive integer"));
      return v.get<int>();
    };
    cfg.solver.threadsScan = count("scan", cfg.solver.threadsScan, true);
    cfg.solver.threadsSlices = count("slices", cfg.solver.threadsSlices, true);
    cfg.solver.threadsBlas = count("blas", cfg.solver.threadsBlas, false);
  }
  // working precision (Plan section "M16")
  cfg.solver.precision = j.value("precision", cfg.solver.precision);
  if (cfg.solver.precision != "auto") (void)linalg::parsePrecision(cfg.solver.precision);
  if (j.contains("precision_fallback")) cfg.solver.precisionFallback = j["precision_fallback"].get<bool>() ? 1 : 0;
  if (j.contains("growth_limit")) cfg.solver.growthLimit = j["growth_limit"].get<double>();
  if (j.contains("checkpoint_every")) {
    cfg.solver.checkpointEvery = j["checkpoint_every"].get<int>();
    FG_REQUIRE(cfg.solver.checkpointEvery >= 1, "checkpoint_every must be >= 1");
  }
  const std::string reduction = j.value("reduction", std::string("sequential"));
  FG_REQUIRE(reduction == "sequential" || reduction == "tree", "reduction must be \"sequential\" or \"tree\"");
  cfg.solver.treeReduction = reduction == "tree";

  if (!cylindrical) {
    cfg.stack.dx_um = j.at("period_um").at("x").get<double>();
    cfg.stack.dz_um = dim == "2D" ? j.at("period_um").at("z").get<double>() : 0.0;
  }
  FG_REQUIRE(dim == "1D" || cfg.stack.dz_um > 0, "2D gratings need period_um.z > 0");
  for (const auto& [name, mj] : j.at("materials").items()) cfg.stack.materials[name] = parseMaterial(name, mj);
  cfg.stack.superstrate = j.at("superstrate").get<std::string>();
  cfg.stack.substrate = j.at("substrate").get<std::string>();
  (void)cfg.stack.material(cfg.stack.superstrate);
  (void)cfg.stack.material(cfg.stack.substrate);
  const json& inc = j.at("incidence");
  if (cylindrical) {
    // Plan section "M20: Schichttypen, CAD und Konfiguration"
    const json& cj = j["cylindrical"];
    CylindricalSpec& c = cfg.stack.cylindrical;
    c.enabled = true;
    FG_REQUIRE(cj.value("symmetry", std::string("axial")) == "axial", "cylindrical.symmetry must be \"axial\" (M20)");
    if (cj.contains("harmonics") && !(cj["harmonics"].is_string() && cj["harmonics"].get<std::string>() == "auto"))
      c.harmonics = cj["harmonics"].get<std::vector<int>>();
    if (cj.contains("spectrum")) {
      const json& sp = cj["spectrum"];
      c.kappaMax = sp.value("kappa_max", c.kappaMax);
      c.dkFactor = sp.value("dk_factor", c.dkFactor);
      c.nodesPerPanel = sp.value("nodes_per_panel", c.nodesPerPanel);
      if (sp.contains("deform")) {
        const json& d = sp["deform"];
        if (d.is_boolean() && !d.get<bool>()) c.deform = CylindricalSpec::Deform::Off;
        else if (d.is_object()) c.deform = CylindricalSpec::Deform::On, c.deltaRel = d.value("delta", c.deltaRel);
        else FG_REQUIRE(d.is_string() && d.get<std::string>() == "auto",
                        "cylindrical.spectrum.deform must be \"auto\", false or {\"delta\": ...}");
      }
      FG_REQUIRE(c.kappaMax > 0 && c.dkFactor > 0 && c.nodesPerPanel >= 2, "cylindrical.spectrum: invalid values");
    }
    const std::string fac = cj.value("factorization", std::string("li_r"));
    FG_REQUIRE(fac == "li_r" || fac == "laurent", "cylindrical.factorization must be \"li_r\" or \"laurent\"");
    c.liR = fac == "li_r";
    FG_REQUIRE(inc.value("type", std::string()) == "gaussian_beam",
               "cylindrical: incidence.type must be \"gaussian_beam\" (plane waves follow in M20)");
    c.waist_um = inc.at("waist_um").get<double>();
    c.focus_y_um = inc.value("focus_y_um", 0.0);
    c.polarization = inc.value("polarization", c.polarization);
    FG_REQUIRE(c.waist_um > 0, "incidence.waist_um must be > 0");
    for (const char* a : {"phi_deg", "psi_deg", "delta_deg"})
      FG_REQUIRE(!inc.contains(a), std::string("gaussian_beam: ") + a + " is not used (beam on the axis)");
  }
  // scan: cartesian product lambda x phi x psi x delta (each a number, an array or a range)
  for (double lam : parseRange(inc.at("lambda_um")))
    for (double phi : parseRange(inc.value("phi_deg", json(0.0))))
      for (double psi : parseRange(inc.value("psi_deg", json(0.0))))
        for (double delta : parseRange(inc.value("delta_deg", json(0.0))))
          cfg.scan.push_back({lam, phi * DEG, psi * DEG, delta * DEG});

  // layers; a cad block expands into its bands (M14), relative STEP paths against baseDir
  for (const auto& lj : j.value("layers", json::array())) {
    if (lj.value("type", std::string()) == "cad") {
      const cad::CadOptions o = parseCadOptions(lj, baseDir, cfg);
      for (auto& le : cad::slice(o, cfg.stack)) cfg.stack.layers.push_back(std::move(le));
      continue;
    }
    cfg.stack.layers.push_back(parseLayer(lj));
  }
  // coordinate transformation (Plan section "M17"): {"x": {"type": "asr", "eta": 0.8, "breakpoints_um": "auto"}}
  if (j.contains("coordinate_transform")) {
    const json& ct = j["coordinate_transform"];
    FG_REQUIRE(ct.is_object(), "coordinate_transform must be an object {\"x\", \"z\"}");
    for (const auto& [key, aj] : ct.items()) {
      FG_REQUIRE(key == "x" || key == "z", "coordinate_transform: unknown axis '" + key + "' (x or z)");
      AxisTransformSpec& a = key == "x" ? cfg.stack.transform.x : cfg.stack.transform.z;
      const std::string type = aj.value("type", std::string("asr"));
      FG_REQUIRE(type == "asr" || type == "none", "coordinate_transform." + key + ".type must be \"asr\" or \"none\"");
      a.type = type == "asr" ? AxisTransformSpec::Type::ASR : AxisTransformSpec::Type::None;
      a.eta = aj.value("eta", 0.0);
      FG_REQUIRE(a.eta >= 0 && a.eta < 1, "coordinate_transform." + key + ".eta must satisfy 0 <= eta < 1");
      if (aj.contains("breakpoints_um") && !aj["breakpoints_um"].is_string())
        a.breakpoints_um = aj["breakpoints_um"].get<std::vector<double>>();
      else
        FG_REQUIRE(aj.value("breakpoints_um", std::string("auto")) == "auto",
                   "coordinate_transform." + key + ".breakpoints_um must be \"auto\" or an array");
    }
    FG_REQUIRE(dim == "2D" || cfg.stack.transform.z.type == AxisTransformSpec::Type::None,
               "coordinate_transform.z needs a 2D grating");
    (void)resolveTransform(cfg.stack);  // reports missing breakpoints early
  }

  if (j.contains("aperiodic")) {
    const json& ap = j["aperiodic"];
    FG_REQUIRE(ap.is_object(), "aperiodic must be an object {\"pml\", \"excitation\"}");
    AperiodicSpec& a = cfg.stack.aperiodic;
    a.enabled = true;
    if (ap.contains("pml"))
      for (const auto& [key, pj] : ap["pml"].items()) {
        FG_REQUIRE(key == "x" || key == "z", "aperiodic.pml: unknown axis '" + key + "' (x or z)");
        PmlAxisSpec& p = key == "x" ? a.x : a.z;
        p.width_um = pj.at("width_um").get<double>();
        p.beta = pj.value("beta", 1.0);
        FG_REQUIRE(p.width_um > 0 && p.beta >= 0, "aperiodic.pml." + key + ": width_um > 0 and beta >= 0 required");
      }
    FG_REQUIRE(a.x.active() || a.z.active(), "aperiodic: at least one PML axis (pml.x or pml.z)");
    FG_REQUIRE(dim == "2D" || !a.z.active(), "aperiodic.pml.z needs a 2D grating");
    const std::string form = ap.value("formulation", std::string("auto"));
    FG_REQUIRE(form == "auto" || form == "scattered",
               "aperiodic.formulation must be \"auto\" (total field at normal incidence) or \"scattered\"");
    if (form == "scattered") cfg.solver.scatteredField = true;
    if (ap.contains("excitation") && ap["excitation"].is_object()) {  // {"mode": 0} or {"n_eff": 1.8}
      const json& ex = ap["excitation"];
      a.modeExcitation = true;
      a.mode = ex.value("mode", 0);
      a.nEff = ex.value("n_eff", 0.0);
      FG_REQUIRE(a.mode >= 0 && a.nEff >= 0, "aperiodic.excitation: mode >= 0, n_eff > 0");
    } else {
      FG_REQUIRE(ap.value("excitation", std::string("plane_wave")) == "plane_wave",
                 "aperiodic.excitation must be \"plane_wave\", {\"mode\": k} or {\"n_eff\": n}");
    }
    // structured half spaces (Plan "Einfall (b)"): y-invariant layer descriptions without thickness
    for (const char* key : {"superstrate_structure", "substrate_structure"})
      if (j.contains(key)) {
        json lj = j[key];
        if (!lj.contains("thickness_um")) lj["thickness_um"] = 1.0;
        (std::string(key) == "superstrate_structure" ? cfg.stack.superstrateStructure : cfg.stack.substrateStructure) =
            parseLayer(lj);
      }
    checkPmlZones(cfg.stack);  // reports structures in the zones early
  }

  if (j.contains("nearfield") && j["nearfield"].value("enabled", false)) {
    cfg.nearField = true;
    const json& nf = j["nearfield"];
    cfg.field.horizontal_y_um = nf.value("horizontal_y_um", 0.0);
    cfg.field.vertical_z_um = nf.value("vertical_z_um", 0.0);
    if (nf.contains("y_range_um")) {
      cfg.field.y_min_um = nf["y_range_um"][0].get<double>();
      cfg.field.y_max_um = nf["y_range_um"][1].get<double>();
    }
    cfg.field.nx = nf.value("nx", cfg.field.nx);
    cfg.field.ny = nf.value("ny", cfg.field.ny);
    cfg.field.nz = nf.value("nz", cfg.field.nz);
    cfg.field.periods = nf.value("periods", cfg.field.periods);
    cfg.field.reconstruction =
        nf.value("reconstruction", std::string("fff")) == "direct" ? Reconstruction::Direct : Reconstruction::FFF;
    if (nf.contains("quantities")) cfg.field.quantities = nf["quantities"].get<std::vector<std::string>>();
    const std::string field = nf.value("field", std::string("total"));
    FG_REQUIRE(field == "total" || field == "scattered", "nearfield.field must be \"total\" or \"scattered\"");
    cfg.field.scattered = field == "scattered";
  }
  if (j.contains("output")) {
    cfg.outputDir = j["output"].value("dir", cfg.outputDir);
    cfg.efficiencyFile = j["output"].value("efficiencies", cfg.efficiencyFile);
  }
  return cfg;
}

RunConfig loadConfig(const std::string& path) {
  std::ifstream in(path);
  FG_REQUIRE(in.good(), "cannot open configuration file " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return parseConfig(ss.str(), std::filesystem::path(path).parent_path().string());
}

namespace {
void writeHeader(std::ofstream& out, bool sweep) {
  if (sweep) out << "sweep_parameter,sweep_value,";
  out << "lambda_um,phi_deg,psi_deg,delta_deg,side,n,m,efficiency,eff_s,eff_p,amp_s_re,amp_s_im,amp_p_re,"
         "amp_p_im\n";
  out << std::setprecision(12);
}
void writeRows(std::ofstream& out, const std::string& prefix, const std::vector<Wavelength>& scan,
               const std::vector<Result>& results) {
  for (std::size_t k = 0; k < results.size(); ++k) {
    const auto& w = scan[k];
    auto row = [&](const char* side, const OrderEfficiency& o) {
      out << prefix << w.lambda_um << ',' << w.phi / DEG << ',' << w.psi / DEG << ',' << w.delta / DEG << ','
          << side << ',' << o.n << ',' << o.m << ',' << o.efficiency << ',' << o.eff_s << ',' << o.eff_p << ','
          << o.amp_s.real() << ',' << o.amp_s.imag() << ',' << o.amp_p.real() << ',' << o.amp_p.imag() << '\n';
    };
    for (const auto& o : results[k].reflected) row("R", o);
    for (const auto& o : results[k].transmitted) row("T", o);
    const ScatteringResult& sc = results[k].scattering;  // aperiodic (M18): widths in the efficiency column
    if (sc.present)
      for (const auto& [side, v] : {std::pair<const char*, double>{"S_up", sc.up}, {"S_down", sc.down}, {"S_abs", sc.absorbed}})
        row(side, OrderEfficiency{0, 0, v, 0, 0});
    const ModalResult& mr = results[k].modal;  // modal excitation (M18): guided-mode powers, n = mode number
    if (mr.present) {
      for (const auto& m : mr.reflected) row("MR", OrderEfficiency{m.index, 0, m.power, m.amplitude, 0, 0, 0});
      for (const auto& m : mr.transmitted) row("MT", OrderEfficiency{m.index, 0, m.power, m.amplitude, 0, 0, 0});
      row("Mrad", OrderEfficiency{0, 0, mr.radiated, 0, 0});
    }
    const CylindricalResult& cr = results[k].cylindrical;  // M20: R, T, A of the beam; per harmonic n = m
    if (cr.present) {
      for (const auto& [side, v] : {std::pair<const char*, double>{"C_R", cr.reflected}, {"C_T", cr.transmitted},
                                    {"C_A", cr.absorbed}})
        row(side, OrderEfficiency{0, 0, v, 0, 0});
      for (std::size_t i = 0; i < cr.harmonics.size(); ++i) {
        row("C_Rm", OrderEfficiency{cr.harmonics[i], 0, cr.reflectedPerHarmonic[i], 0, 0});
        row("C_Tm", OrderEfficiency{cr.harmonics[i], 0, cr.transmittedPerHarmonic[i], 0, 0});
      }
    }
  }
}
std::ofstream openCsv(const std::string& file) {
  const auto dir = std::filesystem::path(file).parent_path();
  if (!dir.empty()) std::filesystem::create_directories(dir);
  std::ofstream out(file);
  FG_REQUIRE(out.good(), "cannot write " + file);
  return out;
}
}  // namespace

void writeEfficiencies(const std::string& file, const std::vector<Wavelength>& scan,
                       const std::vector<Result>& results) {
  std::ofstream out = openCsv(file);
  writeHeader(out, false);
  writeRows(out, "", scan, results);
}

void writeEfficiencies(const std::string& file, const std::vector<RunConfig>& runs,
                       const std::vector<std::vector<Result>>& results) {
  FG_REQUIRE(runs.size() == results.size(), "writeEfficiencies: size mismatch");
  const bool sweep = !runs.empty() && !runs.front().sweepParameter.empty();
  std::ofstream out = openCsv(file);
  writeHeader(out, sweep);
  for (std::size_t r = 0; r < runs.size(); ++r) {
    std::ostringstream pre;
    pre << std::setprecision(12);
    if (sweep) pre << runs[r].sweepParameter << ',' << runs[r].sweepValue << ',';
    writeRows(out, pre.str(), runs[r].scan, results[r]);
  }
}

std::vector<RunConfig> parseSweep(const std::string& text, const std::string& baseDir) {
  const json j = json::parse(text);
  if (!j.contains("sweep")) return {parseConfig(text, baseDir)};
  const json& sw = j["sweep"];
  const std::string param = sw.at("parameter").get<std::string>();
  const json::json_pointer ptr(param);
  FG_REQUIRE(j.contains(ptr), "sweep: parameter '" + param + "' not found in the configuration");
  FG_REQUIRE(j.at(ptr).is_number(), "sweep: parameter '" + param + "' is not a number");
  std::vector<RunConfig> out;
  for (double v : parseRange(sw.at("values"))) {
    json jk = j;
    jk.erase("sweep");
    jk[ptr] = v;
    RunConfig c = parseConfig(jk.dump(), baseDir);
    c.sweepParameter = param;
    c.sweepValue = v;
    out.push_back(std::move(c));
  }
  FG_REQUIRE(!out.empty(), "sweep: no values");
  return out;
}

std::vector<RunConfig> loadSweep(const std::string& path) {
  std::ifstream in(path);
  FG_REQUIRE(in.good(), "cannot open configuration file " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return parseSweep(ss.str(), std::filesystem::path(path).parent_path().string());
}

namespace {
/// The configuration without its cad blocks (context: periods, materials, scan) and the cad blocks.
std::pair<RunConfig, std::vector<cad::CadOptions>> cadBlocks(const std::string& text, const std::string& baseDir,
                                                             const std::string& stepText) {
  json j = json::parse(text);
  std::vector<json> blocks;
  json rest = json::array();
  for (const auto& lj : j.value("layers", json::array())) {
    if (lj.value("type", std::string()) == "cad") blocks.push_back(lj);
    else rest.push_back(lj);
  }
  j["layers"] = rest;
  RunConfig cfg = parseConfig(j.dump(), baseDir);
  std::vector<cad::CadOptions> out;
  for (auto b : blocks) {
    if (!stepText.empty()) b["step"] = stepText;
    out.push_back(parseCadOptions(b, baseDir, cfg));
  }
  FG_REQUIRE(!out.empty(), "configuration has no layer of type \"cad\"");
  return {std::move(cfg), std::move(out)};
}
}  // namespace

std::string cadReport(const std::string& text, const std::string& baseDir, const std::string& stepText) {
  const auto [cfg, blocks] = cadBlocks(text, baseDir, stepText);
  json out = json::array();
  for (const auto& o : blocks) out.push_back(json::parse(cad::check(o, cfg.stack).json()));
  return out.size() == 1 ? out[0].dump(2) : out.dump(2);
}

std::string cadPreview(const std::string& text, const std::string& baseDir, const std::string& stepText, double y_um) {
  const auto [cfg, blocks] = cadBlocks(text, baseDir, stepText);
  return cad::previewSvg(blocks.front(), cfg.stack, y_um);
}

}  // namespace fg::io
