// Dispersive material models (Plan section "Dispersive Materialien (M22)"): V23b (models) and V23g (errors).
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>

#include "fg/core/error.hpp"
#include "fg/geometry/dispersion.hpp"
#include "fg/io/config.hpp"

using json = nlohmann::json;
using fg::cplx;
using fg::Dispersion;

namespace {

// fused silica after Malitson [D1]
Dispersion::Sellmeier malitson() {
  Dispersion::Sellmeier s;
  s.A = 1.0;
  s.B = {0.6961663, 0.4079426, 0.8974794};
  s.C_um2 = {0.0684043 * 0.0684043, 0.1162414 * 0.1162414, 9.896161 * 9.896161};
  return s;
}

// aluminium after Rakic et al. [D3] (Lorentz-Drude)
Dispersion::DrudeLorentz rakicAl() {
  Dispersion::DrudeLorentz d;
  d.epsInf = 1.0;
  d.omegaP_eV = 14.98;
  d.f0 = 0.523;
  d.gamma0_eV = 0.047;
  d.oscillators = {{0.227, 0.162, 0.333}, {0.050, 1.544, 0.312}, {0.166, 1.808, 1.351}, {0.030, 3.473, 3.382}};
  return d;
}

Dispersion::Table table(std::vector<double> x, std::vector<double> a, std::vector<double> b, bool cubic) {
  Dispersion::Table t;
  t.lambda_um = std::move(x);
  t.a = std::move(a);
  t.b = std::move(b);
  t.cubic = cubic;
  t.inlineData = true;
  return t;
}

/// Central difference of f at x (step h).
template <class F>
cplx central(F f, double x, double h) {
  return (f(x + h) - f(x - h)) / (2 * h);
}

std::string minimal(const json& materials, const json& incidence = {{"lambda_um", 0.6}},
                    const std::string& sup = "air") {
  json j{{"dimension", "1D"},
         {"truncation", {{"Nx", 3}}},
         {"period_um", {{"x", 1.0}}},
         {"incidence", incidence},
         {"materials", materials},
         {"superstrate", sup},
         {"substrate", "sub"},
         {"layers", json::array()}};
  return j.dump();
}

/// Message of the fg::Error thrown by f ("" if none).
template <class F>
std::string errorOf(F f) {
  try {
    f();
  } catch (const fg::Error& e) {
    return e.what();
  }
  return "";
}

}  // namespace

TEST(Dispersion, SellmeierMalitson) {
  const Dispersion d("/materials/sio2/n", Dispersion::Quantity::Eps, malitson());
  for (double lam : {0.25, 0.4, 0.5876, 1.0, 1.55, 3.0}) {
    const double l2 = lam * lam;
    const double ref = 1 + 0.6961663 * l2 / (l2 - 0.0684043 * 0.0684043) +
                       0.4079426 * l2 / (l2 - 0.1162414 * 0.1162414) + 0.8974794 * l2 / (l2 - 9.896161 * 9.896161);
    EXPECT_NEAR(d.eps(lam).real(), ref, 1e-14 * ref) << lam;
    EXPECT_EQ(d.eps(lam).imag(), 0.0);
  }
  EXPECT_NEAR(std::sqrt(d.eps(0.5876).real()), 1.4585, 1e-4);  // published value (Malitson, n_d)
}

TEST(Dispersion, DrudeLorentzLossy) {
  const Dispersion d("/materials/al/eps", Dispersion::Quantity::Eps, rakicAl());
  for (double lam = 0.2; lam <= 3.0; lam += 0.01) EXPECT_GT(d.eps(lam).imag(), 0.0) << lam;
  EXPECT_LT(d.eps(0.8).real(), -10.0);  // metallic
}

TEST(Dispersion, Cauchy) {
  Dispersion::Cauchy c;
  c.A = 1.5;
  c.B_um2 = 0.004;
  c.C_um4 = 1e-4;
  c.k = 0.01;
  const Dispersion d("/materials/x/n", Dispersion::Quantity::N, c);
  const double lam = 0.55;
  const cplx n(1.5 + 0.004 / (lam * lam) + 1e-4 / std::pow(lam, 4), 0.01);
  EXPECT_NEAR(std::abs(d.value(lam) - n), 0.0, 1e-15);
  EXPECT_NEAR(std::abs(d.eps(lam) - n * n), 0.0, 1e-14);
  // a constant model is exactly the number (Plan V23a)
  Dispersion::Cauchy a;
  a.A = 1.5;
  EXPECT_EQ(Dispersion("/m/n", Dispersion::Quantity::N, a).value(0.7), cplx(1.5, 0.0));
}

TEST(Dispersion, TableLinearAndCubic) {
  const std::vector<double> x{0.4, 0.47, 0.6, 0.66, 0.8, 1.05};
  auto poly = [](double l) { return 2.0 - 1.5 * l + 0.8 * l * l - 0.3 * l * l * l; };
  std::vector<double> a, b;
  for (double l : x) a.push_back(poly(l)), b.push_back(0.5 * poly(l));
  const Dispersion lin("/m/eps", Dispersion::Quantity::Eps, table(x, a, b, false));
  const Dispersion cub("/m/eps", Dispersion::Quantity::Eps, table(x, a, b, true));
  for (std::size_t i = 0; i < x.size(); ++i) {  // exact at the nodes
    EXPECT_EQ(lin.value(x[i]), cplx(a[i], b[i]));
    EXPECT_EQ(cub.value(x[i]), cplx(a[i], b[i]));
  }
  for (std::size_t i = 0; i + 1 < x.size(); ++i) {  // linear exact in the middle
    const double m = 0.5 * (x[i] + x[i + 1]);
    EXPECT_NEAR(lin.value(m).real(), 0.5 * (a[i] + a[i + 1]), 1e-15);
  }
  for (double l = 0.4; l <= 1.05; l += 0.013) {  // not-a-knot reproduces cubic polynomials
    EXPECT_NEAR(cub.value(l).real(), poly(l), 1e-13) << l;
    EXPECT_NEAR(cub.value(l).imag(), 0.5 * poly(l), 1e-13) << l;
  }
  // n, k under "n": eps = (n + i k)^2
  const Dispersion nk("/m/n", Dispersion::Quantity::N, table(x, a, b, false));
  EXPECT_EQ(nk.eps(x[2]), cplx(a[2], b[2]) * cplx(a[2], b[2]));
  // range
  EXPECT_THROW((void)lin.value(0.39), fg::Error);
  EXPECT_THROW((void)lin.value(1.06), fg::Error);
  auto t = table(x, a, b, false);
  t.extrapolateConstant = true;
  const Dispersion ex("/m/eps", Dispersion::Quantity::Eps, t);
  EXPECT_EQ(ex.value(0.3), cplx(a.front(), b.front()));
  EXPECT_EQ(ex.value(2.0), cplx(a.back(), b.back()));
  EXPECT_TRUE(ex.extrapolated(2.0));
  EXPECT_FALSE(ex.extrapolated(0.5));
}

TEST(Dispersion, WavelengthDerivatives) {
  std::vector<Dispersion> models{
      Dispersion("/a/eps", Dispersion::Quantity::Eps, malitson()),
      Dispersion("/b/eps", Dispersion::Quantity::Eps, rakicAl()),
      Dispersion("/c/n", Dispersion::Quantity::N, Dispersion::Cauchy{1.6, 0.01, 2e-4, 0.02}),
      Dispersion("/d/n", Dispersion::Quantity::N,
                 table({0.4, 0.5, 0.65, 0.8, 1.0}, {1.2, 0.9, 0.5, 0.3, 0.25}, {2.0, 2.6, 3.4, 4.6, 6.0}, false)),
      Dispersion("/e/n", Dispersion::Quantity::N,
                 table({0.4, 0.5, 0.65, 0.8, 1.0}, {1.2, 0.9, 0.5, 0.3, 0.25}, {2.0, 2.6, 3.4, 4.6, 6.0}, true)),
  };
  for (const auto& d : models)
    for (double lam : {0.43, 0.57, 0.71, 0.93}) {  // tables: away from the nodes
      const cplx num = central([&](double l) { return d.eps(l); }, lam, 1e-6);
      const cplx ana = d.dEpsDLambda(lam);
      EXPECT_LT(std::abs(ana - num), 1e-8 * std::max(1.0, std::abs(ana))) << d.pointer() << " " << lam;
    }
}

TEST(Dispersion, ParameterDerivatives) {
  // derivative of the natural value with respect to model numbers against central differences
  const double lam = 0.63, h = 1e-6;
  auto checkParam = [&](auto make, const std::vector<std::string>& tok, auto set, Dispersion::Quantity q) {
    const Dispersion d("/m/v", q, make(0.0));
    const cplx num = (Dispersion("/m/v", q, make(h)).value(lam) - Dispersion("/m/v", q, make(-h)).value(lam)) / (2 * h);
    const cplx ana = d.dValueDParam(tok, lam);
    EXPECT_LT(std::abs(ana - num), 1e-8 * std::max(1.0, std::abs(ana))) << tok[1];
    (void)set;
  };
  using Q = Dispersion::Quantity;
  checkParam([](double e) { auto s = malitson(); s.B[1] += e; return s; }, {"sellmeier", "B", "1"}, 0, Q::Eps);
  checkParam([](double e) { auto s = malitson(); s.C_um2[0] += e; return s; }, {"sellmeier", "C_um2", "0"}, 0, Q::Eps);
  checkParam([](double e) { auto s = malitson(); s.A += e; return s; }, {"sellmeier", "A"}, 0, Q::Eps);
  checkParam([](double e) { return Dispersion::Cauchy{1.6 + e, 0.01, 2e-4, 0.02}; }, {"cauchy", "A"}, 0, Q::N);
  checkParam([](double e) { return Dispersion::Cauchy{1.6, 0.01 + e, 2e-4, 0.02}; }, {"cauchy", "B_um2"}, 0, Q::N);
  checkParam([](double e) { return Dispersion::Cauchy{1.6, 0.01, 2e-4 + e, 0.02}; }, {"cauchy", "C_um4"}, 0, Q::N);
  checkParam([](double e) { return Dispersion::Cauchy{1.6, 0.01, 2e-4, 0.02 + e}; }, {"cauchy", "k"}, 0, Q::N);
  checkParam([](double e) { auto d = rakicAl(); d.epsInf += e; return d; }, {"drude_lorentz", "eps_inf"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.omegaP_eV += e; return d; }, {"drude_lorentz", "omega_p_eV"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.f0 += e; return d; }, {"drude_lorentz", "f0"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.gamma0_eV += e; return d; }, {"drude_lorentz", "gamma0_eV"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.oscillators[2].f += e; return d; },
             {"drude_lorentz", "oscillators", "2", "f"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.oscillators[1].omega_eV += e; return d; },
             {"drude_lorentz", "oscillators", "1", "omega_eV"}, 0, Q::Eps);
  checkParam([](double e) { auto d = rakicAl(); d.oscillators[0].gamma_eV += e; return d; },
             {"drude_lorentz", "oscillators", "0", "gamma_eV"}, 0, Q::Eps);
  for (bool cubic : {false, true}) {
    const std::vector<double> x{0.4, 0.5, 0.65, 0.8, 1.0};
    auto mk = [&](std::size_t col, double e) {
      std::vector<double> a{1.2, 0.9, 0.5, 0.3, 0.25}, b{2.0, 2.6, 3.4, 4.6, 6.0};
      (col == 1 ? a : b)[2] += e;
      return table(x, a, b, cubic);
    };
    checkParam([&](double e) { return mk(1, e); }, {"table", "data", "2", "1"}, 0, Q::N);
    checkParam([&](double e) { return mk(2, e); }, {"table", "data", "2", "2"}, 0, Q::N);
    const Dispersion d("/m/v", Q::N, mk(1, 0.0));
    EXPECT_THROW((void)d.dValueDParam({"table", "data", "2", "0"}, lam), fg::Error);  // wavelengths: no parameters
  }
}

TEST(Dispersion, UnitsAndFiles) {
  // the same table in um, nm and eV (eV rows in any order) from inline data and from a CSV file
  const std::vector<double> x{0.4, 0.5, 0.65, 0.8};
  json um = json::array(), nm = json::array(), ev = json::array();
  std::ostringstream csv;
  csv << "# test table\nwavelength_nm;n;k\n";
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double n = 1.5 + 0.1 * static_cast<double>(i), k = 0.01 * static_cast<double>(i);
    um.push_back({x[i], n, k});
    nm.push_back({x[i] * 1000, n, k});
    ev.push_back({fg::HC_EV_UM / x[i], n, k});
    csv << x[i] * 1000 << "; " << n << "; " << k << "\n";
  }
  const auto dir = std::filesystem::temp_directory_path() / "fg_test_dispersion";
  std::filesystem::create_directories(dir);
  std::ofstream(dir / "t.csv") << csv.str();
  const json mats{{"air", {{"n", 1.0}}},
                  {"sub", {{"n", 1.5}}},
                  {"um", {{"n", {{"table", {{"data", um}}}}}}},
                  {"nm", {{"n", {{"table", {{"data", nm}, {"unit", "nm"}}}}}}},
                  {"ev", {{"n", {{"table", {{"data", ev}, {"unit", "eV"}}}}}}},
                  {"file", {{"n", {{"table", {{"file", "t.csv"}, {"unit", "nm"}}}}}}}};
  const auto cfg = fg::io::parseConfig(minimal(mats), dir.string());
  for (double lam : {0.4, 0.46, 0.7, 0.8}) {
    const cplx ref = fg::io::materialAt(cfg.stack, "um", lam).front().eps;
    for (const char* m : {"nm", "ev", "file"})
      EXPECT_NEAR(std::abs(fg::io::materialAt(cfg.stack, m, lam).front().eps - ref), 0.0, 1e-14) << m << " " << lam;
  }
}

// V23g: errors before anything is computed, with the JSON pointer in the message
TEST(Dispersion, ErrorsBeforeTheComputation) {
  const json table{{"table", {{"data", {{0.3, 1.0, 0.5}, {1.9, 0.5, 4.0}}}}}};
  const json mats{{"air", {{"n", 1.0}}}, {"sub", {{"n", 1.5}}}, {"au", {{"n", table}}}};
  // lambda outside a table
  const json scan{{"lambda_um", {{"from", 0.25}, {"to", 0.5}, {"steps", 2}}}};
  EXPECT_NE(errorOf([&] { (void)fg::io::parseConfig(minimal(mats, scan)); })
                .find("/materials/au/n/table: lambda 0.25 um outside [0.3, 1.9] um"),
            std::string::npos);
  // Sellmeier pole
  json pole = mats;
  pole["sio2"] = {{"n", {{"sellmeier", {{"B", {1.0}}, {"C_um2", {0.36}}}}}}};
  EXPECT_NE(errorOf([&] { (void)fg::io::parseConfig(minimal(pole, {{"lambda_um", 0.6}})); })
                .find("/materials/sio2/n/sellmeier: lambda 0.6 um at the pole"),
            std::string::npos);
  // lossy superstrate at a scan point
  json sup = mats;
  sup["cover"] = {{"n", {{"cauchy", {{"A", 1.4}, {"k", 0.01}}}}}};
  EXPECT_NE(errorOf([&] { (void)fg::io::parseConfig(minimal(sup, {{"lambda_um", 0.6}}, "cover")); })
                .find("/materials/cover/n: the superstrate must be lossless"),
            std::string::npos);
  // file and data at once
  json both = mats;
  both["au"]["n"]["table"]["file"] = "x.csv";
  EXPECT_NE(errorOf([&] { (void)fg::io::parseConfig(minimal(both)); }).find("/materials/au/n/table: exactly one of"),
            std::string::npos);
  // model in eps_tensor
  json ten = mats;
  ten["co"] = {{"eps_tensor", {{table, 0, 0}, {0, 2, 0}, {0, 0, 2}}}};
  EXPECT_NE(errorOf([&] { (void)fg::io::parseConfig(minimal(ten)); }).find("/materials/co/eps_tensor/0/0"),
            std::string::npos);
  // a sweep of the wavelength is checked for every value
  json j = json::parse(minimal(mats));
  j["sweep"] = {{"parameter", "/incidence/lambda_um"}, {"values", {0.5, 2.0}}};
  EXPECT_NE(errorOf([&] { (void)fg::io::parseSweep(j.dump()); }).find("lambda 2 um outside"), std::string::npos);
}
