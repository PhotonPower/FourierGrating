// Validation M22 (Plan section "Dispersive Materialien: Tests und Abnahme"): V23a (bit-identical constant models),
// V23c (scan equivalence), V23d (physics), V23e (gradients), V23f (model cache). V23b and V23g are unit tests
// (tests/unit/test_dispersion.cpp), V23h is a Python test (tests/python/test_materials.py).
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <nlohmann/json.hpp>

#include "fg/adjoint/gradient.hpp"
#include "fg/cad/cad.hpp"
#include "fg/geometry/dispersion.hpp"
#include "fg/io/config.hpp"
#include "fg/solvers/simulation.hpp"

using json = nlohmann::json;
using fg::cplx;

namespace {

const json kSio2 = {{"sellmeier",
                     {{"A", 1.0},
                      {"B", {0.6961663, 0.4079426, 0.8974794}},
                      {"C_um2", {0.00467914826, 0.0135120631, 97.9340025}}}}};
// silver after Rakic et al. [D3] (Lorentz-Drude)
const json kAg = {{"drude_lorentz",
                   {{"eps_inf", 1.0},
                    {"omega_p_eV", 9.01},
                    {"f0", 0.845},
                    {"gamma0_eV", 0.048},
                    {"oscillators",
                     {{{"f", 0.065}, {"omega_eV", 0.816}, {"gamma_eV", 3.886}},
                      {{"f", 0.124}, {"omega_eV", 4.481}, {"gamma_eV", 0.452}},
                      {{"f", 0.011}, {"omega_eV", 8.185}, {"gamma_eV", 0.065}},
                      {{"f", 0.840}, {"omega_eV", 9.083}, {"gamma_eV", 0.916}},
                      {{"f", 5.646}, {"omega_eV", 20.29}, {"gamma_eV", 2.419}}}}}}};
// aluminium after Rakic et al. [D3]
const json kAl = {{"drude_lorentz",
                   {{"eps_inf", 1.0},
                    {"omega_p_eV", 14.98},
                    {"f0", 0.523},
                    {"gamma0_eV", 0.047},
                    {"oscillators",
                     {{{"f", 0.227}, {"omega_eV", 0.162}, {"gamma_eV", 0.333}},
                      {{"f", 0.050}, {"omega_eV", 1.544}, {"gamma_eV", 0.312}},
                      {{"f", 0.166}, {"omega_eV", 1.808}, {"gamma_eV", 1.351}},
                      {{"f", 0.030}, {"omega_eV", 3.473}, {"gamma_eV", 3.382}}}}}}};
const json kTio2 = {{"cauchy", {{"A", 2.40}, {"B_um2", 0.03}, {"C_um4", 0.002}, {"k", 0.004}}}};
// n, k table (inline, linear)
const json kAu = {{"table",
                   {{"data", {{0.40, 1.47, 1.95}, {0.50, 0.97, 1.87}, {0.60, 0.25, 2.98}, {0.70, 0.16, 3.95},
                              {0.85, 0.16, 5.17}, {1.00, 0.22, 6.35}}}}}};

/// 1D lamellar grating with four dispersive materials (Sellmeier substrate, Cauchy film, Drude-Lorentz bars),
/// conical incidence.
json grating1D() {
  return json{{"dimension", "1D"},
              {"truncation", {{"Nx", 8}}},
              {"period_um", {{"x", 0.9}}},
              {"incidence", {{"lambda_um", {0.45, 0.55, 0.65, 0.8}}, {"phi_deg", 15}, {"psi_deg", 10}, {"delta_deg", 30}}},
              {"materials", {{"air", {{"n", 1.0}}}, {"sio2", {{"n", kSio2}}}, {"tio2", {{"n", kTio2}}}, {"ag", {{"eps", kAg}}},
                             {"au", {{"n", kAu}}}}},
              {"superstrate", "air"},
              {"substrate", "sio2"},
              {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.05}, {"material", "tio2"}},
                                      json{{"type", "lamellar1d"},
                                           {"thickness_um", 0.08},
                                           {"breakpoints_um", {0.0, 0.3, 0.55, 0.9}},
                                           {"materials", {"ag", "air", "au"}}},
                                      json{{"type", "homogeneous"}, {"thickness_um", 0.1}, {"material", "air"}}})}};
}

/// The configuration at lambda with every model written as its number (Plan V23c: "feste Werte n(lambda_i)").
json fixedAt(const json& cfg, double lam) {
  const auto c = fg::io::parseConfig(cfg.dump());
  json f = cfg;
  for (const auto& [name, m] : c.stack.materials) {
    if (m.model) {
      const cplx q = m.model->value(lam);
      json& mj = f["materials"][name];
      mj.erase("n");
      mj.erase("eps");
      mj[m.model->quantity() == fg::Dispersion::Quantity::N ? "n" : "eps"] = {q.real(), q.imag()};
    }
    for (const auto& d : m.principalModels)
      if (d) {
        const cplx e = d->eps(lam);
        f[json::json_pointer(d->pointer())] = {e.real(), e.imag()};
      }
  }
  f["incidence"]["lambda_um"] = lam;
  return f;
}

void expectIdentical(const fg::Result& a, const fg::Result& b, const std::string& what) {
  ASSERT_EQ(a.reflected.size(), b.reflected.size()) << what;
  ASSERT_EQ(a.transmitted.size(), b.transmitted.size()) << what;
  for (std::size_t k = 0; k < a.reflected.size(); ++k) {
    EXPECT_EQ(a.reflected[k].efficiency, b.reflected[k].efficiency) << what;
    EXPECT_EQ(a.reflected[k].amp_s, b.reflected[k].amp_s) << what;
    EXPECT_EQ(a.reflected[k].amp_p, b.reflected[k].amp_p) << what;
  }
  for (std::size_t k = 0; k < a.transmitted.size(); ++k) {
    EXPECT_EQ(a.transmitted[k].efficiency, b.transmitted[k].efficiency) << what;
    EXPECT_EQ(a.transmitted[k].amp_s, b.transmitted[k].amp_s) << what;
  }
  EXPECT_EQ(a.absorption, b.absorption) << what;
  EXPECT_EQ(a.scattering.up, b.scattering.up) << what;
  EXPECT_EQ(a.scattering.down, b.scattering.down) << what;
  EXPECT_EQ(a.cylindrical.reflected, b.cylindrical.reflected) << what;
  EXPECT_EQ(a.cylindrical.transmitted, b.cylindrical.transmitted) << what;
}

fg::Result solveFirst(const json& cfg) {
  const auto c = fg::io::parseConfig(cfg.dump());
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
}

/// Scan equivalence: solveScan of the dispersive configuration = single computations with fixed values.
void scanEquivalence(const json& cfg, const std::string& what) {
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto rs = fg::solveScan(c.solver, c.stack, c.scan);
  for (std::size_t i = 0; i < c.scan.size(); ++i)
    expectIdentical(rs[i], solveFirst(fixedAt(cfg, c.scan[i].lambda_um)), what + " lambda " + std::to_string(c.scan[i].lambda_um));
}

fg::Objective objective(std::initializer_list<std::tuple<char, int, int, double>> w) {
  fg::Objective o;
  for (const auto& [s, n, m, v] : w) o.weights.push_back({fg::OrderRef{s, n, m}, v});
  return o;
}

/// Central difference of L at the JSON pointer (h relative to the value).
double fd(const json& cfg, const fg::Objective& obj, const std::string& ptr, double rel = 1e-6) {
  const json::json_pointer jp(ptr);
  const double v = cfg.at(jp).get<double>(), h = rel * std::max(1.0, std::abs(v));
  json a = cfg, b = cfg;
  a[jp] = v + h;
  b[jp] = v - h;
  return (fg::gradient(a.dump(), obj, {}).value - fg::gradient(b.dump(), obj, {}).value) / (2 * h);
}

void check(const json& cfg, const fg::Objective& obj, const std::vector<std::string>& params, double tol,
           const char* what) {
  const auto g = fg::gradient(cfg.dump(), obj, params);
  for (std::size_t k = 0; k < params.size(); ++k) {
    const double num = fd(cfg, obj, params[k]);
    std::printf("%s %-52s adjoint % .10e  fd % .10e  rel %.2e\n", what, params[k].c_str(), g.grad[k], num,
                std::abs(g.grad[k] - num) / std::max(1e-12, std::abs(num)));
    EXPECT_NEAR(g.grad[k], num, tol * std::max(1.0, std::abs(num))) << what << " " << params[k];
  }
}

}  // namespace

// ---------------------------------------------------------------- V23a: constant models = numbers (bit-identical)
TEST(V23a_Dispersion, ConstantModelsAreBitIdentical) {
  // 1D (lamellar, isotropic): number, Cauchy with A only, constant table (n under "n", eps under "eps")
  json lam{{"dimension", "1D"},
           {"truncation", {{"Nx", 10}}},
           {"period_um", {{"x", 1.0}}},
           {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 10}, {"delta_deg", 30}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", 2.0}}}, {"film", {{"eps", {2.2, 0.1}}}}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.1}, {"material", "film"}},
                                   json{{"type", "lamellar1d"},
                                        {"thickness_um", 0.4},
                                        {"breakpoints_um", {0.0, 0.45, 1.0}},
                                        {"materials", {"bar", "air"}}}})}};
  const auto ref = solveFirst(lam);
  json a = lam, b = lam;
  a["materials"]["bar"]["n"] = {{"cauchy", {{"A", 2.0}}}};
  a["materials"]["glass"]["n"] = {{"cauchy", {{"A", 1.5}, {"B_um2", 0.0}}}};
  b["materials"]["bar"]["n"] = {{"table", {{"data", {{0.5, 2.0}, {0.7, 2.0}}}}}};
  b["materials"]["film"]["eps"] = {{"table", {{"data", {{0.4, 2.2, 0.1}, {0.6, 2.2, 0.1}, {0.9, 2.2, 0.1}}}}}};
  expectIdentical(ref, solveFirst(a), "1D cauchy");
  expectIdentical(ref, solveFirst(b), "1D table");
  // 2D (polygon pillars)
  json pil{{"dimension", "2D"},
           {"truncation", {{"Nx", 3}, {"Nz", 3}}},
           {"period_um", {{"x", 0.8}, {"z", 0.7}}},
           {"incidence", {{"lambda_um", 0.6}, {"phi_deg", 10}, {"psi_deg", 25}, {"delta_deg", 40}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"pillar", {{"eps", {4.0, 0.2}}}}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({json{{"type", "polygon2d"},
                                        {"thickness_um", 0.2},
                                        {"background", "air"},
                                        {"factorization", "laurent"},
                                        {"inclusions", {{{"material", "pillar"},
                                                         {"rectangle_um", {{"center", {0.4, 0.35}}, {"size", {0.4, 0.3}}}}}}}}})}};
  json pt = pil;
  pt["materials"]["pillar"]["eps"] = {{"table", {{"data", {{0.5, 4.0, 0.2}, {0.7, 4.0, 0.2}}}, {"interpolation", "linear"}}}};
  expectIdentical(solveFirst(pil), solveFirst(pt), "2D table");
  // anisotropic (principal values with Euler angles)
  json ani = lam;
  ani["materials"]["lnb"] = {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}};
  ani["layers"][0]["material"] = "lnb";
  json at = ani;
  at["materials"]["lnb"]["eps_principal"][1] = {{"table", {{"data", {{0.5, 4.98}, {0.7, 4.98}}}}}};
  at["materials"]["lnb"]["eps_principal"][2] = {{"sellmeier", {{"A", 4.6}, {"B", json::array()}, {"C_um2", json::array()}}}};
  expectIdentical(solveFirst(ani), solveFirst(at), "anisotropic");
}

// ---------------------------------------------------------------- V23c: scan equivalence
TEST(V23c_Dispersion, ScanEqualsFixedValues) {
  const json cfg = grating1D();
  scanEquivalence(cfg, "rcwa");
  json rk = cfg;
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-6}};
  scanEquivalence(rk, "rk");
  // the diagnostics carry n and eps of the dispersive materials
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan[1]);
  ASSERT_EQ(r.diag.materials.size(), 4u);
  for (const auto& v : r.diag.materials) {
    const cplx e = fg::io::materialAt(c.stack, v.name, 0.55).front().eps;
    EXPECT_EQ(v.eps, e) << v.name;
    EXPECT_NEAR(std::abs(v.n * v.n - v.eps), 0.0, 1e-13 * std::abs(v.eps)) << v.name;
  }
}

TEST(V23c_Dispersion, SimulationSweepNearFieldSMatrix) {
  const json cfg = grating1D();
  // Simulation: any order of wavelengths, also jumping back
  fg::Simulation sim(cfg.dump());
  fg::Wavelength w = sim.config().scan.front();
  for (double lam : {0.55, 0.45, 0.8, 0.45, 0.55, 0.55, 0.65}) {
    w.lambda_um = lam;
    expectIdentical(sim.solve(w), solveFirst(fixedAt(cfg, lam)), "simulation " + std::to_string(lam));
    // S-matrix in the s/p basis of the orders
    const auto S = sim.orderSMatrix(w);
    const auto fc = fg::io::parseConfig(fixedAt(cfg, lam).dump());
    const auto Sf = fg::makeSolver(fc.solver)->orderSMatrix(fc.stack, w);
    EXPECT_EQ((S.R_top - Sf.R_top).cwiseAbs().maxCoeff(), 0.0);
    EXPECT_EQ((S.T_down - Sf.T_down).cwiseAbs().maxCoeff(), 0.0);
  }
  // sweep over the wavelength and over a model value
  json sw = cfg;
  sw["incidence"]["lambda_um"] = 0.6;
  sw["sweep"] = {{"parameter", "/incidence/lambda_um"}, {"values", {0.5, 0.7}}};
  for (const auto& run : fg::io::parseSweep(sw.dump())) {
    const auto rs = fg::solveScan(run.solver, run.stack, run.scan);
    expectIdentical(rs.front(), solveFirst(fixedAt(cfg, run.sweepValue)), "sweep");
  }
  sw["sweep"] = {{"parameter", "/materials/ag/eps/drude_lorentz/omega_p_eV"}, {"values", {8.9, 9.1}}};
  for (const auto& run : fg::io::parseSweep(sw.dump())) {
    json m = cfg;
    m["incidence"]["lambda_um"] = 0.6;
    m["materials"]["ag"]["eps"]["drude_lorentz"]["omega_p_eV"] = run.sweepValue;
    expectIdentical(fg::solveScan(run.solver, run.stack, run.scan).front(), solveFirst(fixedAt(m, 0.6)), "sweep model");
  }
  // near field
  const auto c = fg::io::parseConfig(cfg.dump());
  fg::NearFieldRequest req;
  req.nx = 24;
  req.ny = 12;
  req.y_min_um = -0.05;
  req.y_max_um = 0.3;
  req.horizontal_y_um = 0.1;
  w = c.scan[2];
  const auto nf = fg::makeSolver(c.solver)->nearField(c.stack, w, req);
  const auto fc = fg::io::parseConfig(fixedAt(cfg, w.lambda_um).dump());
  const auto nff = fg::makeSolver(fc.solver)->nearField(fc.stack, w, req);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(nf.vertical.E[i], nff.vertical.E[i]);
    EXPECT_EQ(nf.horizontal.H[i], nff.horizontal.H[i]);
  }
  EXPECT_EQ(nf.vertical.material, nff.vertical.material);
}

TEST(V23c_Dispersion, AnisotropicPmlCylinderCad) {
  // principal values with a Sellmeier entry, 2D, cross-section
  json ani = grating1D();
  const json cauchy = {{"cauchy", {{"A", 2.1}, {"B_um2", 0.02}}}};
  ani["materials"]["lnb"] = {{"eps_principal", json::array({kSio2, 4.98, cauchy})}, {"euler_deg", {20, 70, 10}}};
  ani["layers"][0]["material"] = "lnb";
  scanEquivalence(ani, "principal");
  // isolated ridge with PML (M18), oblique incidence: scattered-field formulation
  json pml = json::parse(R"({"dimension": "1D", "period_um": {"x": 3.0}, "truncation": {"Nx": 24},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"lambda_um": [0.55, 0.65], "phi_deg": 12},
    "layers": [{"type": "lamellar1d", "thickness_um": 0.3, "breakpoints_um": [0.0, 1.25, 1.75, 3.0], "materials": ["air", "bar", "air"]}],
    "aperiodic": {"pml": {"x": {"width_um": 0.45, "beta": 2.0}}, "excitation": "plane_wave"}})");
  pml["materials"]["bar"] = {{"n", kTio2}};
  pml["materials"]["glass"] = {{"n", kSio2}};
  scanEquivalence(pml, "pml");
  // rotationally symmetric pillar (M20), Gaussian beam
  json cyl = json::parse(R"({"cylindrical": {"harmonics": "auto", "spectrum": {"kappa_max": 3}, "factorization": "li_r"},
    "solver": "rcwa", "materials": {"air": {"n": 1.0}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": [0.9, 1.0], "waist_um": 0.5, "focus_y_um": 0.8, "polarization": "x"},
    "layers": [{"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["sin", "air"]}]})");
  cyl["materials"]["glass"] = {{"n", kSio2}};
  cyl["materials"]["sin"] = {{"n", {{"cauchy", {{"A", 1.98}, {"B_um2", 0.012}}}}}};
  scanEquivalence(cyl, "cylinder");
  fg::Simulation sim(cyl.dump());  // the cylindrical path resolves the materials per solve
  fg::Wavelength w = sim.config().scan.back();
  expectIdentical(sim.solve(w), solveFirst(fixedAt(cyl, w.lambda_um)), "cylinder simulation");
  EXPECT_THROW((void)sim.orderSMatrix(w), fg::Error);
  // CAD (M14): needs the slicer
  if (!fg::cad::available()) GTEST_SKIP() << "CAD slicer not built (FG_ENABLE_NURBFORGE=OFF)";
  json cad = json::parse(R"({"dimension": "2D", "period_um": {"x": 0.5, "z": 0.5}, "truncation": {"Nx": 2, "Nz": 2},
    "materials": {"air": {"n": 1.0}, "sio2": {"n": 1.46}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"lambda_um": [0.6, 0.7], "phi_deg": 10, "delta_deg": 45},
    "layers": [{"type": "cad", "file": ")" FG_SOURCE_DIR R"(/examples/cad/pillar_array.step", "unit": "auto",
                "axes": {"y": "+Z", "x": "+X"}, "cell_origin_um": [0.0, -0.5], "background": "superstrate",
                "bodies": [{"name": "Pillar", "material": "si"}, {"name": "Lid", "material": "sio2"}],
                "slicing": {"rcwa_slices": 6}, "normal_field": {"grid": 64}}]})");
  cad["materials"]["glass"] = {{"n", kSio2}};
  cad["materials"]["si"] = {{"n", kTio2}};
  scanEquivalence(cad, "cad");
}

// ---------------------------------------------------------------- V23d: physics
TEST(V23d_Dispersion, DrudeLorentzHalfSpaceIsFresnel) {
  json cfg{{"dimension", "1D"},
           {"truncation", {{"Nx", 2}}},
           {"period_um", {{"x", 0.5}}},
           {"incidence", {{"lambda_um", {{"from", 0.3}, {"to", 2.0}, {"steps", 18}}}, {"phi_deg", 0}, {"delta_deg", 90}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"al", {{"eps", kAl}}}}},
           {"superstrate", "air"},
           {"substrate", "al"},
           {"layers", json::array()}};
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto rs = fg::solveScan(c.solver, c.stack, c.scan);
  for (std::size_t i = 0; i < c.scan.size(); ++i) {
    const cplx n = std::sqrt(fg::io::materialAt(c.stack, "al", c.scan[i].lambda_um).front().eps);
    const double R = std::norm((1.0 - n) / (1.0 + n));
    double R0 = 0;
    for (const auto& o : rs[i].reflected)
      if (o.n == 0) R0 = o.efficiency;
    EXPECT_NEAR(R0, R, 1e-12) << c.scan[i].lambda_um;
  }
}

TEST(V23d_Dispersion, LosslessSellmeierGratingConservesEnergy) {
  json cfg{{"dimension", "1D"},
           {"truncation", {{"Nx", 12}}},
           {"period_um", {{"x", 0.8}}},
           {"incidence", {{"lambda_um", {{"from", 0.4}, {"to", 1.6}, {"steps", 7}}}, {"phi_deg", 25}, {"psi_deg", 15}, {"delta_deg", 30}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"sio2", {{"n", kSio2}}}}},
           {"superstrate", "air"},
           {"substrate", "sio2"},
           {"layers", json::array({json{{"type", "lamellar1d"},
                                        {"thickness_um", 0.35},
                                        {"breakpoints_um", {0.0, 0.4, 0.8}},
                                        {"materials", {"sio2", "air"}}}})}};
  for (const char* solver : {"rcwa", "rk"}) {
    cfg["solver"] = solver;
    const auto c = fg::io::parseConfig(cfg.dump());
    const auto rs = fg::solveScan(c.solver, c.stack, c.scan);
    for (std::size_t i = 0; i < rs.size(); ++i)
      EXPECT_NEAR(rs[i].diag.energyBalance, 1.0, 1e-10) << solver << " " << c.scan[i].lambda_um;
  }
}

TEST(V23d_Dispersion, RayleighAnomalyInADispersiveSubstrate) {
  // order +1 in the substrate becomes propagating at the zero of lambda - d (nu_1(lambda) - nu_M(lambda) sin phi)
  const double d = 0.5, phi = 10.0;
  json cfg{{"dimension", "1D"},
           {"truncation", {{"Nx", 4}}},
           {"period_um", {{"x", d}}},
           {"incidence", {{"lambda_um", 0.6}, {"phi_deg", phi}, {"delta_deg", 90}}},
           {"materials", {{"water", {{"n", {{"cauchy", {{"A", 1.32}, {"B_um2", 0.003}}}}}}}, {"sio2", {{"n", kSio2}}}}},
           {"superstrate", "water"},
           {"substrate", "sio2"},
           {"layers", json::array({json{{"type", "lamellar1d"},
                                        {"thickness_um", 0.05},
                                        {"breakpoints_um", {0.0, 0.25, 0.5}},
                                        {"materials", {"sio2", "water"}}}})}};
  const auto c = fg::io::parseConfig(cfg.dump());
  auto nu = [&](const char* m, double lam) { return std::sqrt(fg::io::materialAt(c.stack, m, lam).front().eps.real()); };
  auto f = [&](double lam) { return lam - d * (nu("sio2", lam) - nu("water", lam) * std::sin(phi * fg::PI / 180)); };
  double a = 0.4, b = 1.0;
  ASSERT_LT(f(a), 0.0);
  ASSERT_GT(f(b), 0.0);
  while (b - a > 1e-12) ((f(0.5 * (a + b)) < 0) ? a : b) = 0.5 * (a + b);
  const double lamR = 0.5 * (a + b);
  std::printf("Rayleigh anomaly of the order +1 in the substrate: lambda = %.12f um\n", lamR);
  auto hasOrder = [&](double lam) {
    fg::Wavelength w = c.scan.front();
    w.lambda_um = lam;
    for (const auto& o : fg::makeSolver(c.solver)->solve(c.stack, w).transmitted)
      if (o.n == 1) return true;
    return false;
  };
  EXPECT_TRUE(hasOrder(lamR - 1e-10));
  EXPECT_FALSE(hasOrder(lamR + 1e-10));
}

// ---------------------------------------------------------------- V23e: gradients
TEST(V23e_Dispersion, GradientsIsotropic) {
  json cfg = grating1D();
  cfg["incidence"]["lambda_um"] = 0.63;
  cfg["layers"][1]["factorization"] = "li_lamellar";
  // the superstrate dispersive too (lossless Cauchy)
  cfg["materials"]["cover"] = {{"n", {{"cauchy", {{"A", 1.0}, {"B_um2", 0.004}}}}}};
  cfg["superstrate"] = "cover";
  cfg["layers"][2]["material"] = "cover";
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}, {'R', 0, 0, -0.3}});
  const std::vector<std::string> params{"/incidence/lambda_um",
                                        "/materials/sio2/n/sellmeier/B/0",
                                        "/materials/sio2/n/sellmeier/C_um2/1",
                                        "/materials/tio2/n/cauchy/A",
                                        "/materials/tio2/n/cauchy/k",
                                        "/materials/ag/eps/drude_lorentz/oscillators/1/f",
                                        "/materials/ag/eps/drude_lorentz/omega_p_eV",
                                        "/materials/ag/eps/drude_lorentz/gamma0_eV",
                                        "/materials/au/n/table/data/2/1",
                                        "/materials/au/n/table/data/3/2",
                                        "/materials/cover/n/cauchy/B_um2",
                                        "/layers/1/thickness_um"};
  check(cfg, obj, params, 1e-6, "rcwa");
  json rk = cfg;
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-10}};
  rk["truncation"] = {{"Nx", 6}};
  check(rk, obj, {"/incidence/lambda_um", "/materials/sio2/n/sellmeier/B/0", "/materials/ag/eps/drude_lorentz/oscillators/1/f",
                  "/materials/au/n/table/data/2/1"},
        1e-6, "rk");
  // cubic table
  json cub = cfg;
  cub["materials"]["au"]["n"]["table"]["interpolation"] = "cubic";
  check(cub, obj, {"/incidence/lambda_um", "/materials/au/n/table/data/2/1", "/materials/au/n/table/data/3/1", "/materials/au/n/table/data/4/2"}, 1e-6,
        "cubic");
  // wavelength gradient = fixed-material part + chain rule (Plan eq:disp_grad_lambda): the fixed part alone differs
  const auto g = fg::gradient(cfg.dump(), obj, {"/incidence/lambda_um"});
  const auto gf = fg::gradient(fixedAt(cfg, 0.63).dump(), obj, {"/incidence/lambda_um"});
  EXPECT_GT(std::abs(g.grad[0] - gf.grad[0]), 1e-3 * std::abs(g.grad[0]));
  // not parameters: table wavelengths, the model object, an imaginary part of the superstrate
  EXPECT_THROW((void)fg::gradient(cfg.dump(), obj, {"/materials/au/n/table/data/2/0"}), fg::Error);
  EXPECT_THROW((void)fg::gradient(cfg.dump(), obj, {"/materials/cover/n/cauchy/k"}), fg::Error);
}

TEST(V23e_Dispersion, GradientsPrincipalAndSimulation) {
  json cfg = grating1D();
  cfg["incidence"]["lambda_um"] = 0.63;
  const json cauchy = {{"cauchy", {{"A", 2.1}, {"B_um2", 0.02}, {"k", 0.01}}}};
  cfg["materials"]["lnb"] = {{"eps_principal", json::array({kSio2, 4.98, cauchy})}, {"euler_deg", {20, 70, 10}}};
  cfg["layers"][0]["material"] = "lnb";
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}});
  const std::vector<std::string> params{"/incidence/lambda_um", "/materials/lnb/eps_principal/0/sellmeier/B/2",
                                        "/materials/lnb/eps_principal/2/cauchy/B_um2",
                                        "/materials/lnb/eps_principal/2/cauchy/k", "/materials/lnb/euler_deg/1"};
  check(cfg, obj, params, 1e-6, "principal");
  // Simulation: valueAndGrad with the cached models = gradient(), also after a change of wavelength
  fg::Simulation sim(cfg.dump());
  fg::Wavelength w = sim.config().scan.front();
  for (double lam : {0.63, 0.7, 0.63}) {
    w.lambda_um = lam;
    json c = cfg;
    c["incidence"]["lambda_um"] = lam;
    const auto gs = sim.valueAndGrad(w, obj, params);
    const auto gg = fg::gradient(c.dump(), obj, params);
    EXPECT_EQ(gs.value, gg.value);
    for (std::size_t k = 0; k < params.size(); ++k) EXPECT_EQ(gs.grad[k], gg.grad[k]) << params[k];
  }
}

// ---------------------------------------------------------------- V23f: model cache of Simulation
TEST(V23f_Dispersion, SimulationRebuildsOnlyDispersiveLayers) {
  const json cfg = grating1D();  // layers 0 (tio2) and 1 (ag/au) dispersive, layer 2 (air) not
  fg::Simulation sim(cfg.dump());
  EXPECT_EQ(sim.modelBuilds(), 3);
  fg::Wavelength w = sim.config().scan.front();
  (void)sim.solve(w);
  EXPECT_EQ(sim.modelBuilds(), 3);
  w.lambda_um = 0.6;
  (void)sim.solve(w);
  EXPECT_EQ(sim.modelBuilds(), 5);
  (void)sim.solve(w);
  EXPECT_EQ(sim.modelBuilds(), 5);
  w.lambda_um = 0.45;
  (void)sim.solve(w);
  EXPECT_EQ(sim.modelBuilds(), 7);
  // without dispersion (and with a dispersive superstrate/substrate only) a change of wavelength rebuilds nothing
  json fixed = fixedAt(cfg, 0.5);
  fixed["materials"]["sio2"] = {{"n", kSio2}};
  fg::Simulation s2(fixed.dump());
  w.lambda_um = 0.5;
  (void)s2.solve(w);
  w.lambda_um = 0.7;
  (void)s2.solve(w);
  EXPECT_EQ(s2.modelBuilds(), 3);
  json one{{"dimension", "1D"},
           {"truncation", {{"Nx", 4}}},
           {"period_um", {{"x", 1.0}}},
           {"incidence", {{"lambda_um", 0.6}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.1}, {"material", "glass"}}})}};
  fg::Simulation s3(one.dump());
  for (double lam : {0.5, 0.7, 0.9}) {
    w.lambda_um = lam;
    (void)s3.solve(w);
  }
  EXPECT_EQ(s3.modelBuilds(), 1);
  // diagnostics: extrapolation warning
  json ex = cfg;
  ex["materials"]["au"]["n"]["table"]["extrapolate"] = "constant";
  fg::Simulation s4(ex.dump());
  w.lambda_um = 1.2;
  const auto r = s4.solve(w);
  bool warned = false;
  for (const auto& m : r.diag.warnings) warned = warned || m.find("/materials/au/n/table: lambda 1.2 um outside") != std::string::npos;
  EXPECT_TRUE(warned);
}
