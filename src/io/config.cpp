#include "fg/io/config.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
#include "fg/geometry/dispersion.hpp"
#include "fg/geometry/nurbs_profile.hpp"
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

/// Table rows from a CSV file (Plan sec:disp_config): separators comma, semicolon or white space; lines
/// with '#' and a header line without numbers are skipped.
std::vector<std::vector<double>> readTableFile(const std::string& file, const std::string& where) {
  std::ifstream in(file);
  FG_REQUIRE(in.good(), where + ": cannot open table file " + file);
  std::vector<std::vector<double>> rows;
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (const auto c = line.find('#'); c != std::string::npos) line.erase(c);
    for (char& ch : line)
      if (ch == ',' || ch == ';' || ch == '\t' || ch == '\r') ch = ' ';
    std::istringstream ls(line);
    std::vector<double> row;
    std::string tok;
    bool numeric = true;
    while (ls >> tok) {
      std::size_t pos = 0;
      double v = 0;
      try {
        v = std::stod(tok, &pos);
      } catch (const std::exception&) {
        pos = 0;
      }
      if (pos != tok.size()) {
        numeric = false;
        break;
      }
      row.push_back(v);
    }
    if (row.empty() && numeric) continue;  // empty line
    if (!numeric) {
      FG_REQUIRE(header && rows.empty(), where + ": non-numeric line in " + file + ": '" + line + "'");
      header = false;
      continue;
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

/// A dispersion model (Plan sec:disp_models); `table` is the quantity of a table under this key.
std::shared_ptr<const Dispersion> parseDispersion(const json& v, const std::string& ptr, Dispersion::Quantity table,
                                                  const std::string& baseDir) {
  FG_REQUIRE(v.is_object() && v.size() == 1,
             ptr + ": a model is {\"sellmeier\" | \"cauchy\" | \"drude_lorentz\" | \"table\": {...}}");
  const std::string key = v.begin().key();
  const json& m = v.at(key);
  const std::string where = ptr + "/" + key;
  FG_REQUIRE(m.is_object(), where + " must be an object");
  auto known = [&](std::initializer_list<const char*> keys) {
    for (const auto& [k, x] : m.items()) {
      (void)x;
      bool ok = false;
      for (const char* a : keys) ok = ok || k == a;
      FG_REQUIRE(ok, where + ": unknown entry '" + k + "'");
    }
  };
  using Q = Dispersion::Quantity;
  if (key == "sellmeier") {
    known({"A", "B", "C_um2", "range_um"});
    Dispersion::Sellmeier s;
    s.A = m.value("A", 1.0);
    s.B = m.at("B").get<std::vector<double>>();
    s.C_um2 = m.at("C_um2").get<std::vector<double>>();
    if (m.contains("range_um")) {
      const auto r = m["range_um"].get<std::vector<double>>();
      FG_REQUIRE(r.size() == 2 && 0 <= r[0] && r[0] < r[1], where + ": range_um must be [lo, hi] with 0 <= lo < hi");
      s.rangeLo_um = r[0];
      s.rangeHi_um = r[1];
    }
    return std::make_shared<Dispersion>(ptr, Q::Eps, s);
  }
  if (key == "cauchy") {
    known({"A", "B_um2", "C_um4", "k"});
    Dispersion::Cauchy c;
    c.A = m.at("A").get<double>();
    c.B_um2 = m.value("B_um2", 0.0);
    c.C_um4 = m.value("C_um4", 0.0);
    c.k = m.value("k", 0.0);
    FG_REQUIRE(c.k >= 0, where + ": k must be >= 0 (losses, time dependence exp(-i omega t))");
    return std::make_shared<Dispersion>(ptr, Q::N, c);
  }
  if (key == "drude_lorentz") {
    known({"eps_inf", "omega_p_eV", "f0", "gamma0_eV", "oscillators"});
    Dispersion::DrudeLorentz d;
    d.epsInf = m.value("eps_inf", 1.0);
    d.omegaP_eV = m.at("omega_p_eV").get<double>();
    d.f0 = m.value("f0", 0.0);
    d.gamma0_eV = m.value("gamma0_eV", 0.0);
    FG_REQUIRE(d.gamma0_eV >= 0, where + ": gamma0_eV must be >= 0");
    for (const auto& o : m.value("oscillators", json::array())) {
      d.oscillators.push_back({o.at("f").get<double>(), o.at("omega_eV").get<double>(), o.at("gamma_eV").get<double>()});
      FG_REQUIRE(d.oscillators.back().gamma_eV >= 0, where + ": gamma_eV must be >= 0");
    }
    return std::make_shared<Dispersion>(ptr, Q::Eps, d);
  }
  FG_REQUIRE(key == "table", ptr + ": unknown model '" + key + "' (sellmeier | cauchy | drude_lorentz | table)");
  known({"file", "data", "unit", "interpolation", "extrapolate"});
  FG_REQUIRE(m.contains("file") != m.contains("data"), where + ": exactly one of 'file' and 'data'");
  Dispersion::Table t;
  std::vector<std::vector<double>> rows;
  if (m.contains("data")) {
    t.inlineData = true;
    rows = m["data"].get<std::vector<std::vector<double>>>();
  } else {
    std::string file = m["file"].get<std::string>();
    if (std::filesystem::path(file).is_relative() && !baseDir.empty())
      file = (std::filesystem::path(baseDir) / file).string();
    rows = readTableFile(file, where);
  }
  const std::string unit = m.value("unit", std::string("um"));
  FG_REQUIRE(unit == "um" || unit == "nm" || unit == "eV", where + ": unit must be \"um\", \"nm\" or \"eV\"");
  const std::string interp = m.value("interpolation", std::string("linear"));
  FG_REQUIRE(interp == "linear" || interp == "cubic", where + ": interpolation must be \"linear\" or \"cubic\"");
  t.cubic = interp == "cubic";
  const std::string ex = m.value("extrapolate", std::string("error"));
  FG_REQUIRE(ex == "error" || ex == "constant", where + ": extrapolate must be \"error\" or \"constant\"");
  t.extrapolateConstant = ex == "constant";
  FG_REQUIRE(!rows.empty(), where + ": no rows");
  const std::size_t cols = rows.front().size();
  FG_REQUIRE(cols == 2 || cols == 3, where + ": rows need 2 or 3 columns (lambda, a[, b])");
  std::vector<std::size_t> order(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    FG_REQUIRE(rows[i].size() == cols, where + ": row " + std::to_string(i) + " has a different number of columns");
    FG_REQUIRE(rows[i][0] > 0, where + ": row " + std::to_string(i) + ": the first column must be positive");
    order[i] = i;
  }
  auto lam = [&](std::size_t i) {
    const double x = rows[i][0];
    return unit == "um" ? x : unit == "nm" ? x * 1e-3 : HC_EV_UM / x;
  };
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return lam(a) < lam(b); });
  for (std::size_t i : order) {
    t.lambda_um.push_back(lam(i));
    t.a.push_back(rows[i][1]);
    t.b.push_back(cols == 3 ? rows[i][2] : 0.0);
    t.row.push_back(static_cast<int>(i));
  }
  return std::make_shared<Dispersion>(ptr, table, std::move(t));
}

Material parseMaterial(const std::string& name, const json& j, const std::string& baseDir) {
  using Q = Dispersion::Quantity;
  const std::string base = "/materials/" + name;
  if (j.contains("n") || j.contains("eps")) {
    const bool isN = j.contains("n");
    const json& v = isN ? j["n"] : j["eps"];
    if (v.is_object()) {  // M22: dispersive value
      Material m = Material::fromIndex(name, 1.0);
      m.model = parseDispersion(v, base + (isN ? "/n" : "/eps"), isN ? Q::N : Q::Eps, baseDir);
      m.eps = cplx(std::numeric_limits<double>::quiet_NaN(), 0.0);
      m.epsTensor = m.eps * Mat3C::Identity();
      return m;
    }
    if (isN) return Material::fromIndex(name, toComplex(v));
    Material m = Material::fromIndex(name, 1.0);
    m.eps = toComplex(v);
    m.epsTensor = m.eps * Mat3C::Identity();
    return m;
  }
  if (j.contains("eps_tensor")) {
    const json& t = j["eps_tensor"];
    FG_REQUIRE(t.is_array() && t.size() == 3, "eps_tensor must be 3x3");
    Mat3C e;
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) {
        const json& v = t[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
        if (v.is_object())
          throw NotImplemented(base + "/eps_tensor/" + std::to_string(i) + "/" + std::to_string(k) +
                               ": dispersive eps_tensor entries (M22: use eps_principal with euler_deg)");
        e(i, k) = toComplex(v);
      }
    return Material::fromTensor(name, e);
  }
  if (j.contains("eps_principal")) {
    const json& p = j["eps_principal"];
    FG_REQUIRE(p.is_array() && p.size() == 3, "eps_principal needs three values");
    std::vector<double> eu{0.0, 0.0, 0.0};
    if (j.contains("euler_deg")) eu = j["euler_deg"].get<std::vector<double>>();
    FG_REQUIRE(eu.size() == 3, "euler_deg needs three angles (z-x-z)");
    std::array<std::shared_ptr<const Dispersion>, 3> models;
    std::array<cplx, 3> fixed{};
    for (std::size_t i = 0; i < 3; ++i) {
      if (p[i].is_object())
        models[i] = parseDispersion(p[i], base + "/eps_principal/" + std::to_string(i), Q::Eps, baseDir);
      else
        fixed[i] = toComplex(p[i]);
    }
    const std::array<double, 3> euler{eu[0] * DEG, eu[1] * DEG, eu[2] * DEG};
    if (models[0] || models[1] || models[2]) {
      Material m = Material::fromTensor(name, Mat3C::Constant(cplx(std::numeric_limits<double>::quiet_NaN(), 0.0)));
      m.principalModels = models;
      m.principalAxes = true;
      m.principal = fixed;
      m.euler = euler;
      return m;
    }
    Material m = Material::fromTensor(name, rotatedTensor(fixed[0], fixed[1], fixed[2], euler[0], euler[1], euler[2]));
    m.principalAxes = true;
    m.principal = fixed;
    m.euler = euler;
    return m;
  }
  throw Error("material '" + name + "' needs 'n', 'eps', 'eps_tensor' or 'eps_principal'");
}

/// Early check of dispersive materials at every scan wavelength (Plan sec:disp_config, "Fruehe Pruefung"):
/// valid range, Sellmeier poles and the lossless superstrate, before anything is computed.
void checkDispersion(const RunConfig& cfg) {
  for (const auto& [name, m] : cfg.stack.materials) {
    if (!m.dispersive()) continue;
    for (const auto& w : cfg.scan) {
      const Material r = m.at(w.lambda_um);  // throws with the pointer
      if (name == cfg.stack.superstrate) {
        const cplx e = r.isotropic ? r.eps : r.epsTensor(0, 0);
        if (!(r.isotropic && e.imag() == 0.0 && e.real() > 0.0)) {
          std::ostringstream msg;
          msg << (m.model ? m.model->pointer() : "/materials/" + name) << ": the superstrate must be lossless at every "
              << "scan point (lambda " << w.lambda_um << " um: eps = " << e.real() << (e.imag() < 0 ? " - " : " + ")
              << std::abs(e.imag()) << "i)";
          throw Error(msg.str());
        }
      }
    }
  }
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
  // solids of revolution (M20): only in a cylindrical stack, and there only these
  o.axisymmetric = j.value("axisymmetric", false);
  FG_REQUIRE(o.axisymmetric == cfg.stack.cylindrical.enabled,
             o.axisymmetric ? "cad.axisymmetric needs a \"cylindrical\" stack (M20)"
                            : "cad in a \"cylindrical\" stack needs \"axisymmetric\": true (a solid of revolution)");
  const bool oneD = !o.axisymmetric && cfg.stack.dz_um <= 0;  // Plan section "M23: Konfiguration"
  if (j.contains("cell_origin_um")) {
    const json& c = j["cell_origin_um"];
    if (oneD) {
      FG_REQUIRE(c.is_number() || (c.is_array() && c.size() == 1),
                 "cad.cell_origin_um: a 1D grating has only x0 (a number or [x0])");
      o.cellOrigin_um = std::array<double, 2>{c.is_number() ? c.get<double>() : c[0].get<double>(), 0.0};
    } else {
      FG_REQUIRE(c.is_array() && c.size() == 2, "cad.cell_origin_um: [x0, z0]");
      o.cellOrigin_um = std::array<double, 2>{c[0].get<double>(), c[1].get<double>()};
    }
  }
  if (j.contains("z_cut_um")) {
    FG_REQUIRE(oneD, "cad.z_cut_um is only used for 1D gratings");
    o.zCut_um = j["z_cut_um"].get<double>();
  }
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

/// profile1d with "shape": "nurbs" (Plan section "Erweiterung ab Version 2.37: NURBS-Profil in 1D-Gittern (M24)"):
/// one layer, a Sectioned1DLayer with the NURBS line provider and its band boundaries. Errors carry the pointer
/// of the layer and a code [NURBS_*].
LayerEntry parseNurbsProfile(const json& j, const std::string& ptr, double dx_um, bool is1D, bool cylindrical,
                             bool aperiodic) {
  auto fail = [&](const std::string& code, const std::string& msg) {
    throw Error(ptr + ": [" + code + "] profile1d nurbs: " + msg);
  };
  if (cylindrical) fail("NURBS_STACK", "not in a cylindrical stack");
  if (!is1D) fail("NURBS_STACK", "needs a 1D grating (\"dimension\": \"1D\")");
  if (aperiodic) fail("NURBS_STACK", "not with a PML (aperiodic)");
  if (!j.contains("nurbs") || !j["nurbs"].is_object()) fail("NURBS_INPUT", "needs a \"nurbs\" object");
  const json& nj = j["nurbs"];
  geometry::NurbsCurve c;
  try {
    c.degree = nj.at("degree").get<int>();
    c.knots = nj.at("knots").get<std::vector<double>>();
    for (const auto& pj : nj.at("points_um")) {
      if (!pj.is_array() || pj.size() != 2) fail("NURBS_INPUT", "points_um entries are [x, y]");
      c.points.push_back({pj[0].get<double>(), pj[1].get<double>()});
    }
    if (nj.contains("weights")) c.weights = nj["weights"].get<std::vector<double>>();
  } catch (const json::exception& e) {
    fail("NURBS_INPUT", std::string("nurbs.degree, knots, points_um (and weights): ") + e.what());
  }
  const std::string below = j.at("below").get<std::string>(), above = j.at("above").get<std::string>();
  std::shared_ptr<const geometry::NurbsProfile> prof;
  try {
    prof = std::make_shared<const geometry::NurbsProfile>(c, dx_um, below, above);
  } catch (const Error& e) {
    throw Error(ptr + ": " + e.what());
  }
  // thickness: depth_um, or the top of the curve; the curve must lie in [0, depth] (Plan "M24: Konfiguration")
  const double depth = j.contains("depth_um") ? j["depth_um"].get<double>() : prof->yMax();
  std::ostringstream range;
  range << "the curve spans y = " << prof->yMin() << " ... " << prof->yMax() << " um";
  if (!(depth > 0)) fail("NURBS_RANGE", "depth_um must be positive (" + range.str() + ")");
  const double tol = 1e-12 * std::max(depth, dx_um);
  if (prof->yMin() < -tol || prof->yMax() > depth + tol)
    fail("NURBS_RANGE", range.str() + ", outside the layer [0, depth_um = " + std::to_string(depth) + "]");
  auto provider = std::make_shared<const NurbsProfileLineProvider>(prof, depth);
  Sectioned1DLayer l;
  l.yBreaks_um = provider->bandBreaks();
  // materials alternate in every band (self-intersection that the flattened check could miss)
  std::vector<double> edges{0.0};
  edges.insert(edges.end(), l.yBreaks_um.begin(), l.yBreaks_um.end());
  edges.push_back(depth);
  for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
    try {
      (void)provider->section(0.5 * (edges[b] + edges[b + 1]));
    } catch (const Error& e) {
      throw Error(ptr + ": " + e.what());
    }
  }
  l.provider = provider;
  l.background = above;
  l.materials = below == above ? std::vector<std::string>{below} : std::vector<std::string>{below, above};
  l.slices = j.value("slices", 0);
  l.chebyshevNodes = j.value("chebyshev_nodes", 0);
  if (l.slices < 0 || l.chebyshevNodes < 0) fail("NURBS_INPUT", "slices and chebyshev_nodes must be >= 0");
  LayerEntry le;
  le.factorization = parseFactorization(j);
  le.spec = l;
  return le;
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
      throw Error("unknown profile shape '" + shape + "' (sinus | triangle | points | nurbs)");
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
  } else if (type == "rings_profile") {  // M20: radii as functions of y
    RingsProfileLayer rp;
    rp.thickness_um = j.at("thickness_um").get<double>();
    rp.materials = j.at("materials").get<std::vector<std::string>>();
    rp.slices = j.value("slices", 0);
    using Shape = RingsProfileLayer::Radius::Shape;
    for (const json& r : j.at("radii")) {
      RingsProfileLayer::Radius R;
      if (r.is_number()) {
        R.r0_um = r.get<double>();
      } else if (r.contains("cone")) {
        const auto c = r["cone"].get<std::vector<double>>();
        FG_REQUIRE(c.size() == 2, "rings_profile: cone [r_bottom_um, r_top_um]");
        R.shape = Shape::Cone;
        R.r0_um = c[0];
        R.r1_um = c[1];
      } else if (r.contains("sphere")) {
        R.shape = Shape::Sphere;
        R.centerY_um = r["sphere"].at("center_y_um").get<double>();
        R.sphereR_um = r["sphere"].at("radius_um").get<double>();
      } else if (r.contains("table")) {
        R.shape = Shape::Table;
        R.ys_um = r["table"].at("y_um").get<std::vector<double>>();
        R.rs_um = r["table"].at("r_um").get<std::vector<double>>();
      } else {
        throw Error("rings_profile: a radius is a number, {\"cone\": [r0, r1]}, {\"sphere\": {...}} or {\"table\": {...}}");
      }
      rp.radii.push_back(std::move(R));
    }
    le.spec = std::move(rp);
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
  for (const auto& [name, mj] : j.at("materials").items()) cfg.stack.materials[name] = parseMaterial(name, mj, baseDir);
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
    // "axial" (M20) or {"rotation": P, "L": n} (M21, Plan sec:cylP "CAD und Konfiguration")
    if (cj.contains("symmetry")) {
      const json& sy = cj["symmetry"];
      if (sy.is_object()) {
        for (const auto& [key, v] : sy.items())
          FG_REQUIRE(key == "rotation" || key == "L", "cylindrical.symmetry: unknown key '" + key + "' (rotation, L)");
        FG_REQUIRE(sy.contains("rotation") && sy["rotation"].is_number_integer() && sy["rotation"].get<int>() >= 1,
                   "cylindrical.symmetry.rotation must be an integer P >= 1 (M21)");
        c.rotation = sy["rotation"].get<int>();
        if (sy.contains("L")) {
          FG_REQUIRE(sy["L"].is_number_integer() && sy["L"].get<int>() >= 0, "cylindrical.symmetry.L must be an integer >= 0");
          c.classL = sy["L"].get<int>();
        }
      } else {
        FG_REQUIRE(sy.is_string() && sy.get<std::string>() == "axial",
                   "cylindrical.symmetry must be \"axial\" (M20) or {\"rotation\": P} (M21)");
      }
    }
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
      FG_REQUIRE(c.deltaRel > 0, "cylindrical.spectrum.deform.delta must be > 0 (the path dips into Im k < 0)");
    }
    const std::string fac = cj.value("factorization", std::string("li_r"));
    FG_REQUIRE(fac == "li_r" || fac == "laurent", "cylindrical.factorization must be \"li_r\" or \"laurent\"");
    c.liR = fac == "li_r";
    c.opticalTheorem = cj.value("optical_theorem", false);  // Plan v2.25: control of sigma_ext (volume integrals)
    if (cj.contains("far_field")) {  // Plan "M20: Ausgaben": far field on (theta nodes) x n_phi azimuths
      c.farFieldPhi = cj["far_field"].value("n_phi", 72);
      FG_REQUIRE(c.farFieldPhi >= 1, "cylindrical.far_field.n_phi must be >= 1");
    }
    const std::string type = inc.value("type", std::string());
    FG_REQUIRE(type == "gaussian_beam" || type == "plane_wave",
               "cylindrical: incidence.type must be \"gaussian_beam\" or \"plane_wave\"");
    c.planeWave = type == "plane_wave";
    if (c.planeWave) {
      // plane wave: angles as in the Cartesian scan, scattered field (Plan sec:cylexc)
      for (const char* a : {"waist_um", "focus_y_um", "polarization"})
        FG_REQUIRE(!inc.contains(a), std::string("plane_wave: ") + a + " is not used (polarization: delta_deg)");
    } else {
      c.waist_um = inc.at("waist_um").get<double>();
      c.focus_y_um = inc.value("focus_y_um", 0.0);
      c.polarization = inc.value("polarization", c.polarization);
      FG_REQUIRE(c.waist_um > 0, "incidence.waist_um must be > 0");
      for (const char* a : {"phi_deg", "psi_deg", "delta_deg"})
        FG_REQUIRE(!inc.contains(a), std::string("gaussian_beam: ") + a + " is not used (beam on the axis)");
    }
  }
  // scan: cartesian product lambda x phi x psi x delta (each a number, an array or a range)
  for (double lam : parseRange(inc.at("lambda_um")))
    for (double phi : parseRange(inc.value("phi_deg", json(0.0))))
      for (double psi : parseRange(inc.value("psi_deg", json(0.0))))
        for (double delta : parseRange(inc.value("delta_deg", json(0.0))))
          cfg.scan.push_back({lam, phi * DEG, psi * DEG, delta * DEG});
  checkDispersion(cfg);  // M22: before anything is built or computed

  // layers; a cad block expands into its bands (M14), relative STEP paths against baseDir
  std::size_t layerIndex = 0;
  for (const auto& lj : j.value("layers", json::array())) {
    const std::string ptr = "/layers/" + std::to_string(layerIndex++);
    if (lj.value("type", std::string()) == "cad") {
      const cad::CadOptions o = parseCadOptions(lj, baseDir, cfg);
      for (auto& le : cad::slice(o, cfg.stack)) cfg.stack.layers.push_back(std::move(le));
      continue;
    }
    if (lj.value("type", std::string()) == "profile1d" && lj.value("shape", std::string()) == "nurbs") {  // M24
      cfg.stack.layers.push_back(
          parseNurbsProfile(lj, ptr, cfg.stack.dx_um, dim == "1D", cylindrical, j.contains("aperiodic")));
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
    cfg.field.extent_um = nf.value("extent_um", 0.0);  // cylindrical stacks (M20)
    FG_REQUIRE(cfg.field.extent_um >= 0, "nearfield.extent_um must be >= 0 (0: 1.5 R_s)");
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
    cfg.materialsFile = j["output"].value("materials", cfg.materialsFile);
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
      for (const auto& [side, v] : {std::pair<const char*, double>{"S_up", sc.up}, {"S_down", sc.down}, {"S_lat", sc.lateral},
                                    {"S_abs", sc.absorbed}})
        row(side, OrderEfficiency{0, 0, v, 0, 0});
    const ModalResult& mr = results[k].modal;  // modal excitation (M18): guided-mode powers, n = mode number
    if (mr.present) {
      for (const auto& m : mr.reflected) row("MR", OrderEfficiency{m.index, 0, m.power, m.amplitude, 0, 0, 0});
      for (const auto& m : mr.transmitted) row("MT", OrderEfficiency{m.index, 0, m.power, m.amplitude, 0, 0, 0});
      row("Mrad", OrderEfficiency{0, 0, mr.radiated, 0, 0});
    }
    const CylindricalResult& cr = results[k].cylindrical;  // M20: R, T, A of the beam; per harmonic n = m
    if (cr.present && cr.planeWave) {  // plane wave: cross sections (um^2) up, down, total; per harmonic
      for (const auto& [side, v] : {std::pair<const char*, double>{"C_Sup", cr.scatteredUp_um2},
                                    {"C_Sdown", cr.scatteredDown_um2}, {"C_Ssca", cr.scattering_um2},
                                    {"C_Sabs", cr.absorption_um2}, {"C_Sext", cr.extinction_um2},
                                    {"C_Sext_ot", cr.extinctionTheorem_um2}})
        if (!std::isnan(v)) row(side, OrderEfficiency{0, 0, v, 0, 0});  // C_Sext_ot only when evaluated
      for (std::size_t i = 0; i < cr.harmonics.size(); ++i) {
        row("C_Supm", OrderEfficiency{cr.harmonics[i], 0, cr.reflectedPerHarmonic[i], 0, 0});
        row("C_Sdownm", OrderEfficiency{cr.harmonics[i], 0, cr.transmittedPerHarmonic[i], 0, 0});
      }
    } else if (cr.present) {
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

std::vector<std::string> writeCylindricalSpectra(const std::string& dir, const std::vector<RunConfig>& runs,
                                                 const std::vector<std::vector<Result>>& results) {
  FG_REQUIRE(runs.size() == results.size(), "writeCylindricalSpectra: size mismatch");
  bool any = false;
  for (const auto& rs : results)
    for (const auto& r : rs) any = any || !r.cylindrical.farField.empty();
  if (!any) return {};
  const bool sweep = !runs.front().sweepParameter.empty();
  const auto fa = (std::filesystem::path(dir) / "angular_spectrum.csv").string(),
             ff = (std::filesystem::path(dir) / "far_field.csv").string();
  std::ofstream oa = openCsv(fa), of = openCsv(ff);
  const std::string head = std::string(sweep ? "sweep_parameter,sweep_value," : "") + "lambda_um,phi_deg,psi_deg,delta_deg,";
  oa << head << "side,m,theta_deg,weight,dP_s,dP_p\n" << std::setprecision(12);
  of << head << "side,theta_deg,phi_k_deg,dP_s,dP_p\n" << std::setprecision(12);
  for (std::size_t r = 0; r < runs.size(); ++r)
    for (std::size_t k = 0; k < results[r].size(); ++k) {
      std::ostringstream pre;
      pre << std::setprecision(12);
      if (sweep) pre << runs[r].sweepParameter << ',' << runs[r].sweepValue << ',';
      const Wavelength& w = runs[r].scan[k];
      pre << w.lambda_um << ',' << w.phi / DEG << ',' << w.psi / DEG << ',' << w.delta / DEG << ',';
      const CylindricalResult& c = results[r][k].cylindrical;
      for (const auto& a : c.angularSpectrum)
        for (std::size_t j = 0; j < a.theta_deg.size(); ++j)
          oa << pre.str() << (a.side > 0 ? "up" : "down") << ',' << a.m << ',' << a.theta_deg[j] << ',' << a.weight[j]
             << ',' << a.s[j] << ',' << a.p[j] << '\n';
      for (const auto& f : c.farField)
        for (std::size_t j = 0; j < f.theta_deg.size(); ++j)
          for (std::size_t l = 0; l < f.phi_deg.size(); ++l)
            of << pre.str() << (f.side > 0 ? "up" : "down") << ',' << f.theta_deg[j] << ',' << f.phi_deg[l] << ','
               << f.s[j * f.phi_deg.size() + l] << ',' << f.p[j * f.phi_deg.size() + l] << '\n';
    }
  return {fa, ff};
}

std::vector<MaterialValue> materialAt(const Stack& stack, const std::string& name, double lambda_um) {
  const Material& m = stack.material(name);
  const Material r = m.at(lambda_um);
  auto nOf = [](cplx e) {  // branch Im n >= 0
    const cplx n = std::sqrt(e);
    return n.imag() < 0 ? -n : n;
  };
  if (r.isotropic) {
    const bool isN = m.model && m.model->quantity() == Dispersion::Quantity::N;
    return {{name, 0, isN ? m.model->value(lambda_um) : nOf(r.eps), r.eps}};
  }
  FG_REQUIRE(r.principalAxes, "material '" + name + "' is given as eps_tensor (no principal values)");
  std::vector<MaterialValue> out;
  for (int i = 0; i < 3; ++i) {
    const cplx e = r.principal[static_cast<std::size_t>(i)];
    out.push_back({name, i + 1, nOf(e), e});
  }
  return out;
}

void writeMaterials(const std::string& file, const RunConfig& cfg) {
  std::ofstream out = openCsv(file);
  out << "lambda_um,material,component,n_re,n_im,eps_re,eps_im\n" << std::setprecision(12);
  std::vector<double> lams;
  for (const auto& w : cfg.scan)
    if (std::find(lams.begin(), lams.end(), w.lambda_um) == lams.end()) lams.push_back(w.lambda_um);
  for (double lam : lams)
    for (const auto& [name, m] : cfg.stack.materials) {
      if (!m.isotropic && !m.principalAxes) continue;  // full tensors have no principal values
      for (const auto& v : materialAt(cfg.stack, name, lam))
        out << lam << ',' << name << ',' << (v.component == 0 ? std::string("iso") : std::to_string(v.component)) << ','
            << v.n.real() << ',' << v.n.imag() << ',' << v.eps.real() << ',' << v.eps.imag() << '\n';
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
