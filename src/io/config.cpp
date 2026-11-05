#include "fg/io/config.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/core/error.hpp"
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
    if (j.contains("normal_field")) {
      const json& nf = j["normal_field"];
      const std::string m = nf.value("method", std::string("tensor"));
      FG_REQUIRE(m == "tensor" || m == "vector", "normal_field.method must be \"tensor\" or \"vector\"");
      g.normalField.method = m == "vector" ? NormalFieldOptions::Method::Vector : NormalFieldOptions::Method::Tensor;
      g.normalField.grid = nf.value("grid", g.normalField.grid);
      g.normalField.sigma = nf.value("sigma", g.normalField.sigma);
      g.normalField.tau = nf.value("tau", g.normalField.tau);
    }
    le.spec = g;
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

RunConfig parseConfig(const std::string& text) {
  const json j = json::parse(text);
  RunConfig cfg;
  const std::string dim = j.at("dimension").get<std::string>();
  FG_REQUIRE(dim == "1D" || dim == "2D", "dimension must be \"1D\" or \"2D\"");
  const std::string solver = j.value("solver", std::string("rcwa"));
  FG_REQUIRE(solver == "rcwa" || solver == "rk", "solver must be \"rcwa\" or \"rk\"");
  cfg.solver.kind = solver == "rk" ? SolverKind::RK : SolverKind::RCWA;
  cfg.solver.Nx = j.at("truncation").value("Nx", 10);
  cfg.solver.Nz = j.at("truncation").value("Nz", 0);
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

  cfg.stack.dx_um = j.at("period_um").at("x").get<double>();
  cfg.stack.dz_um = dim == "2D" ? j.at("period_um").at("z").get<double>() : 0.0;
  FG_REQUIRE(dim == "1D" || cfg.stack.dz_um > 0, "2D gratings need period_um.z > 0");
  for (const auto& [name, mj] : j.at("materials").items()) cfg.stack.materials[name] = parseMaterial(name, mj);
  cfg.stack.superstrate = j.at("superstrate").get<std::string>();
  cfg.stack.substrate = j.at("substrate").get<std::string>();
  (void)cfg.stack.material(cfg.stack.superstrate);
  (void)cfg.stack.material(cfg.stack.substrate);
  for (const auto& lj : j.value("layers", json::array())) cfg.stack.layers.push_back(parseLayer(lj));

  const json& inc = j.at("incidence");
  // scan: cartesian product lambda x phi x psi x delta (each a number, an array or a range)
  for (double lam : parseRange(inc.at("lambda_um")))
    for (double phi : parseRange(inc.value("phi_deg", json(0.0))))
      for (double psi : parseRange(inc.value("psi_deg", json(0.0))))
        for (double delta : parseRange(inc.value("delta_deg", json(0.0))))
          cfg.scan.push_back({lam, phi * DEG, psi * DEG, delta * DEG});

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
  return parseConfig(ss.str());
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

std::vector<RunConfig> parseSweep(const std::string& text) {
  const json j = json::parse(text);
  if (!j.contains("sweep")) return {parseConfig(text)};
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
    RunConfig c = parseConfig(jk.dump());
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
  return parseSweep(ss.str());
}

}  // namespace fg::io
