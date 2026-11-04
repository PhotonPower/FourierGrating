#include "fg/io/config.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/core/error.hpp"

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
  throw Error("material '" + name + "' needs 'n', 'eps' or 'eps_tensor'");
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
  } else if (type == "polygon2d" || type == "profile2d") {
    throw NotImplemented("layer type '" + type + "' in configuration (milestone M7)");
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

  cfg.stack.dx_um = j.at("period_um").at("x").get<double>();
  cfg.stack.dz_um = dim == "2D" ? j.at("period_um").at("z").get<double>() : 0.0;
  for (const auto& [name, mj] : j.at("materials").items()) cfg.stack.materials[name] = parseMaterial(name, mj);
  cfg.stack.superstrate = j.at("superstrate").get<std::string>();
  cfg.stack.substrate = j.at("substrate").get<std::string>();
  (void)cfg.stack.material(cfg.stack.superstrate);
  (void)cfg.stack.material(cfg.stack.substrate);
  for (const auto& lj : j.value("layers", json::array())) cfg.stack.layers.push_back(parseLayer(lj));

  const json& inc = j.at("incidence");
  for (double lam : parseRange(inc.at("lambda_um")))
    for (double phi : parseRange(inc.value("phi_deg", json(0.0))))
      cfg.scan.push_back({lam, phi * DEG, inc.value("psi_deg", 0.0) * DEG, inc.value("delta_deg", 0.0) * DEG});

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

void writeEfficiencies(const std::string& file, const std::vector<Wavelength>& scan,
                       const std::vector<Result>& results) {
  const auto dir = std::filesystem::path(file).parent_path();
  if (!dir.empty()) std::filesystem::create_directories(dir);
  std::ofstream out(file);
  FG_REQUIRE(out.good(), "cannot write " + file);
  out << "lambda_um,phi_deg,psi_deg,delta_deg,side,n,m,efficiency,amp_s_re,amp_s_im,amp_p_re,amp_p_im\n";
  out << std::setprecision(12);
  for (std::size_t k = 0; k < results.size(); ++k) {
    const auto& w = scan[k];
    auto row = [&](const char* side, const OrderEfficiency& o) {
      out << w.lambda_um << ',' << w.phi / DEG << ',' << w.psi / DEG << ',' << w.delta / DEG << ','
          << side << ',' << o.n << ',' << o.m << ',' << o.efficiency << ',' << o.amp_s.real() << ','
          << o.amp_s.imag() << ',' << o.amp_p.real() << ',' << o.amp_p.imag() << '\n';
    };
    for (const auto& o : results[k].reflected) row("R", o);
    for (const auto& o : results[k].transmitted) row("T", o);
  }
}

}  // namespace fg::io
