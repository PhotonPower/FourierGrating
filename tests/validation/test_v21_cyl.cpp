// V21 (M20, Plan section "Tests und Abnahme (M20)"): rotationally symmetric structures.
// V21b: guided modes of a step-index fiber are eigenvalues +-i beta of Plan eq:cyl_M.
#include <gtest/gtest.h>

#include "../support/std_bessel.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>

#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;
namespace cyl = fg::cyl;

namespace {

// Exact eigenvalue equation of the step-index fiber (Snyder/Love), normalized units; kind 0: TE,
// 1: TM (m = 0), 2: hybrid (m >= 1).
double fiberEquation(int m, int kind, double eps1, double eps2, double a, double b) {
  const double U = a * std::sqrt(eps1 - b * b), W = a * std::sqrt(b * b - eps2);
  auto J = [](int n, double x) {
    return fgtest::cyl_bessel_j(double(std::abs(n)), x) * ((n < 0 && n % 2) ? -1 : 1);
  };
  auto K = [](int n, double x) { return fgtest::cyl_bessel_k(double(std::abs(n)), x); };
  const double jj = 0.5 * (J(m - 1, U) - J(m + 1, U)) / (U * J(m, U));
  const double kk = -0.5 * (K(m - 1, W) + K(m + 1, W)) / (W * K(m, W));
  if (kind == 0) return jj + kk;
  if (kind == 1) return eps1 * jj + eps2 * kk;
  const double s = 1 / (U * U) + 1 / (W * W);
  return (jj + kk) * (eps1 * jj + eps2 * kk) - (m * b) * (m * b) * s * s;
}

/// Largest root in (n2, n1) (the fundamental mode of the given kind).
double fiberBeta(int m, int kind, double eps1, double eps2, double a) {
  auto f = [&](double b) { return fiberEquation(m, kind, eps1, eps2, a, b); };
  const double lo = std::sqrt(eps2) + 1e-9, hi = std::sqrt(eps1) - 1e-9;
  const int n = 20000;
  for (int i = n; i > 0; --i) {
    double b0 = lo + (hi - lo) * (i - 1) / n, b1 = lo + (hi - lo) * i / n;
    double f0 = f(b0), f1 = f(b1);
    if (!std::isfinite(f0) || !std::isfinite(f1) || f0 * f1 >= 0) continue;
    for (int it = 0; it < 200 && b1 - b0 > 1e-15; ++it) {
      const double bm = 0.5 * (b0 + b1), fm = f(bm);
      if (fm * f0 > 0)
        b0 = bm, f0 = fm;
      else
        b1 = bm, f1 = fm;
    }
    const double r = 0.5 * (b0 + b1);
    if (std::abs(f(r)) < 1e-6) return r;  // not a pole of J_m(U)
  }
  return std::nan("");
}

/// Error of the eigenvalue closest to beta_exact (Plan eq:cyl_M, reduced: M^2 = diag(PR, RP)).
double betaError(int m, double kmax, cyl::CylFactorization f, double betaExact) {
  const double eps1 = 4, eps2 = 1, a = 2;
  cyl::SpectrumOptions o;
  o.kmax = kmax;
  o.Rs = a;
  o.branchPoints = {std::sqrt(eps2), std::sqrt(eps1)};
  const auto grid = cyl::spectralGrid(o);
  const cyl::Rings rings{{a}, {eps1}, eps2};
  const MatC M = cyl::propagationMatrix(grid.k, cyl::cylQ(grid, m, rings, f));
  const Eigen::Index n2 = M.rows() / 2;
  const MatC PR = M.topRightCorner(n2, n2) * M.bottomLeftCorner(n2, n2);
  const VecC mu = PR.eigenvalues();
  double err = 1e300;
  for (Eigen::Index i = 0; i < mu.size(); ++i)
    err = std::min(err, std::abs(std::sqrt(-mu(i)) - betaExact));
  return err;
}

}  // namespace

TEST(V21b, FiberModesLiRAndLaurent) {
  FG_SKIP_WITHOUT_STD_BESSEL();
  // Plan: Li in r <= 2e-4 at k_max = 8 n_core for TM01, HE11, m = 2; Laurent first order in
  // 1/k_max.
  struct Case {
    int m, kind;
    const char* name;
  };
  for (const Case c : {Case{0, 1, "TM01"}, Case{1, 2, "HE11"}, Case{2, 2, "m=2"}}) {
    const double b = fiberBeta(c.m, c.kind, 4, 1, 2);
    ASSERT_TRUE(std::isfinite(b)) << c.name;
    const double li8 = betaError(c.m, 8, cyl::CylFactorization::LiR, b);
    const double li16 = betaError(c.m, 16, cyl::CylFactorization::LiR, b);
    const double la8 = betaError(c.m, 8, cyl::CylFactorization::Laurent, b);
    const double la16 = betaError(c.m, 16, cyl::CylFactorization::Laurent, b);
    std::printf("V21b %-5s beta=%.12f  Li: %.2e %.2e  Laurent: %.2e %.2e (k_max = 8, 16)\n", c.name,
                b, li8, li16, la8, la16);
    EXPECT_LT(li16, 2e-4) << c.name;
    EXPECT_GT(la8 / la16, 1.6) << c.name;  // first order
    EXPECT_LT(la8 / la16, 2.5) << c.name;
    EXPECT_LT(li16, 0.1 * la16) << c.name;
  }
}

TEST(V21b, TE01IndependentOfFactorization) {
  FG_SKIP_WITHOUT_STD_BESSEL();
  // TE01 has only E_phi (continuous at r = a): all factorizations agree (Plan, Befund)
  const double b = fiberBeta(0, 0, 4, 1, 2);
  const double li = betaError(0, 16, cyl::CylFactorization::LiR, b);
  const double la = betaError(0, 16, cyl::CylFactorization::Laurent, b);
  EXPECT_LT(li, 1e-4);  // truncation at k_max = 16 (4e-5)
  EXPECT_LT(std::abs(li - la), 1e-10);
}

// ---------------------------------------------------------------------------------------------------
// V21a: beam through a layer stack (Plan "Tests und Abnahme (M20)"); the stack is laterally
// invariant, so every spectral node is an independent thin-film problem: reference by the Cartesian
// solver (N = 0) at sin(phi) = k_j, weighted with the incident power per node and polarization
// (Plan eq:cyl_flux_closed).
namespace {

std::string beamConfig(const std::string& film2, const std::string& layers, const std::string& pol,
                       const std::string& fac = "li_r", int kappa = 2) {
  return R"({
    "cylindrical": {"spectrum": {"kappa_max": )" +
         std::to_string(kappa) + R"(}, "factorization": ")" + fac + R"("},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "film": {"n": 1.45}, "film2": )" +
         film2 + R"(},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.8,
                  "polarization": ")" +
         pol + R"("},
    "layers": [)" +
         layers + R"(]})";
}
const std::string kFilms = R"({"type": "homogeneous", "thickness_um": 0.3, "material": "film"},
                              {"type": "homogeneous", "thickness_um": 0.2, "material": "film2"})";

fg::Result solveConfig(const std::string& json, fg::SolverKind kind = fg::SolverKind::RCWA) {
  fg::io::RunConfig cfg = fg::io::parseConfig(json);
  cfg.solver.kind = kind;
  return fg::makeSolver(cfg.solver)->solve(cfg.stack, cfg.scan.front());
}

/// Thin-film power coefficients (R_s, R_p, T_s, T_p) at transverse wave number k (superstrate n =
/// 1).
std::array<double, 4> thinFilm(const std::string& film2, double k) {
  fg::io::RunConfig cfg = fg::io::parseConfig(R"({"dimension": "1D", "truncation": {"Nx": 0},
    "period_um": {"x": 1.0},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "film": {"n": 1.45}, "film2": )" +
                                              film2 + R"(},
    "superstrate": "air", "substrate": "glass", "incidence": {"lambda_um": 1.0},
    "layers": [)" + kFilms + "]}");
  std::array<double, 4> out{};
  for (double delta : {0.0, 90.0}) {
    const fg::Result r = fg::makeSolver(cfg.solver)
                             ->solve(cfg.stack, {1.0, std::asin(k), 0.0, delta * fg::PI / 180});
    for (const auto& o : r.reflected)
      if (o.n == 0) out[0] += o.eff_s, out[1] += o.eff_p;
    for (const auto& o : r.transmitted)
      if (o.n == 0) out[2] += o.eff_s, out[3] += o.eff_p;
  }
  return out;
}

/// Reference R, T of an on-axis beam: same spectral grid as the solver (no rings: R_s = 2 w0, k_max
/// = 2 * 1.5).
std::pair<double, double> beamReference(const std::string& film2, const std::string& pol) {
  const double k0 = 2 * fg::PI, w0 = 0.5 * k0;
  cyl::SpectrumOptions o;
  o.kmax = 2 * 1.5;
  o.Rs = 2 * 0.5 * k0;
  o.branchPoints = {1.5, 1.0};
  const auto grid = cyl::spectralGrid(o);
  double pTE = 0, pTM = 0, r = 0, t = 0;
  for (Eigen::Index j = 0; j < grid.k.size(); ++j) {
    const double k = grid.k(j).real(), w = grid.w(j).real();
    if (k >= 1.0) continue;
    const double beta = std::sqrt(1 - k * k), g = std::exp(-k * k * w0 * w0 / 4);
    // spectra: x -> |s| = |p| = g0 / 2 per harmonic; radial (TM) / azimuthal (TE): order-1 spectrum
    // k g
    double s2 = 0, p2 = 0;
    if (pol == "x") s2 = p2 = 0.25 * g * g;
    if (pol == "radial") p2 = k * k * g * g;
    if (pol == "azimuthal") s2 = k * k * g * g;
    const double te = w * k * beta * s2, tm = w * k * p2 / beta;
    const auto c = thinFilm(film2, k);
    pTE += te;
    pTM += tm;
    r += te * c[0] + tm * c[1];
    t += te * c[2] + tm * c[3];
  }
  return {r / (pTE + pTM), t / (pTE + pTM)};
}

}  // namespace

TEST(V21a, BeamThroughLayerStackEqualsThinFilmPerNode) {
  for (const std::string film2 :
       {std::string(R"({"n": 1.38})"), std::string(R"({"n": [1.38, 0.05]})")})
    for (const std::string pol : {"x", "radial", "azimuthal"}) {
      const fg::Result res = solveConfig(beamConfig(film2, kFilms, pol));
      ASSERT_TRUE(res.cylindrical.present);
      const auto [R, T] = beamReference(film2, pol);
      std::printf("V21a %-9s %s: R = %.15f (ref %.15f)  T = %.15f (ref %.15f)  N_k = %d\n",
                  pol.c_str(), film2.c_str(), res.cylindrical.reflected, R,
                  res.cylindrical.transmitted, T, res.cylindrical.nodes);
      EXPECT_NEAR(res.cylindrical.reflected, R, 1e-12) << pol << film2;
      EXPECT_NEAR(res.cylindrical.transmitted, T, 1e-12) << pol << film2;
      if (film2 == R"({"n": 1.38})")
        EXPECT_NEAR(res.cylindrical.reflected + res.cylindrical.transmitted, 1.0, 1e-12);
    }
}

TEST(V21a, RingsSolverIntegration) {
  // rings with one material = homogeneous layer (Li in r: frak A = 1 exactly)
  const std::string film2 = R"({"n": 1.38})";
  const std::string same = R"({"type": "homogeneous", "thickness_um": 0.3, "material": "film"},
    {"type": "rings", "thickness_um": 0.2, "radii_um": [0.3, 0.6], "materials": ["film2", "film2", "film2"]})";
  const fg::Result a = solveConfig(beamConfig(film2, kFilms, "x")),
                   b = solveConfig(beamConfig(film2, same, "x"));
  EXPECT_NEAR(a.cylindrical.reflected, b.cylindrical.reflected, 1e-12);
  EXPECT_NEAR(a.cylindrical.transmitted, b.cylindrical.transmitted, 1e-12);
  // a dielectric pillar: m = -1 and m = +1 carry the same power (mirror symmetry, Plan sec:cylexc);
  // Laurent conserves energy exactly (diag(wk) [[g]] symmetric), Li in r only to its truncation
  // error (7e-4 at kappa_max = 2, 1.7e-5 at 4, 1.5e-5 at 8); RK = RCWA (y-invariant layer)
  const std::string pillar =
      R"({"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["film", "air"]})";
  for (const std::string fac : {"laurent", "li_r"}) {
    const fg::Result p = solveConfig(beamConfig(film2, pillar, "x", fac, fac == "li_r" ? 4 : 2));
    const auto& c = p.cylindrical;
    ASSERT_EQ(c.harmonics.size(), 2u);
    std::printf("V21a pillar %-7s: R = %.12f T = %.12f R+T-1 = %.1e\n", fac.c_str(), c.reflected,
                c.transmitted, c.reflected + c.transmitted - 1);
    EXPECT_NEAR(c.reflectedPerHarmonic[0], c.reflectedPerHarmonic[1], 1e-12);
    EXPECT_NEAR(c.transmittedPerHarmonic[0], c.transmittedPerHarmonic[1], 1e-12);
    EXPECT_NEAR(c.reflected + c.transmitted, 1.0, fac == "laurent" ? 1e-12 : 5e-5);
    EXPECT_GT(std::abs(c.reflected - a.cylindrical.reflected), 1e-4);  // the pillar does something
    if (fac == "laurent") {
      const fg::Result q = solveConfig(beamConfig(film2, pillar, "x", fac), fg::SolverKind::RK);
      EXPECT_NEAR(q.cylindrical.reflected, c.reflected, 1e-8);
      EXPECT_NEAR(q.cylindrical.transmitted, c.transmitted, 1e-8);
    }
  }
}

// Plane waves (Plan sec:cylexc, "Ebene Welle (Streufeld)"): scattered field with the sources of Plan
// eq:cyl_scat; rotational symmetry of the cross sections, RK = RCWA, V21c (Mie sphere).
namespace {

std::string planeWaveConfig(const std::string& layers, double phi, double psi, double delta,
                            const std::string& fac = "li_r", double kappa = 2, const std::string& harmonics = "\"auto\"",
                            const std::string& sub = "glass") {
  char inc[160];
  std::snprintf(inc, sizeof inc, R"("phi_deg": %.17g, "psi_deg": %.17g, "delta_deg": %.17g)", phi, psi, delta);
  return R"({
    "cylindrical": {"spectrum": {"kappa_max": )" +
         std::to_string(kappa) + R"(}, "factorization": ")" + fac + R"(", "harmonics": )" + harmonics + R"(},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "film": {"n": 1.45}, "sph": {"n": 1.5}},
    "superstrate": "air", "substrate": ")" + sub + R"(",
    "incidence": {"type": "plane_wave", "lambda_um": 1.0, )" + inc + R"(},
    "layers": [)" + layers + R"(]})";
}
const std::string kPillar =
    R"({"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["film", "air"]})";

}  // namespace

TEST(V21PlaneWave, RotationalSymmetryAndRkEqualsRcwa) {
  // normal incidence: m = +-1, any polarization gives the same cross section
  const fg::Result a = solveConfig(planeWaveConfig(kPillar, 0, 0, 0)), b = solveConfig(planeWaveConfig(kPillar, 0, 0, 90));
  const auto& ca = a.cylindrical;
  ASSERT_TRUE(ca.present && ca.planeWave);
  ASSERT_EQ(ca.harmonics, (std::vector<int>{-1, 1}));
  std::printf("pillar, normal incidence: sigma_up = %.12f sigma_down = %.12f um^2\n", ca.scatteredUp_um2,
              ca.scatteredDown_um2);
  EXPECT_GT(ca.scattering_um2, 1e-3);
  EXPECT_NEAR(b.cylindrical.scattering_um2 / ca.scattering_um2, 1.0, 1e-12);
  EXPECT_NEAR(ca.reflectedPerHarmonic[0], ca.reflectedPerHarmonic[1], 1e-12 * ca.scattering_um2);
  EXPECT_NEAR(ca.transmittedPerHarmonic[0], ca.transmittedPerHarmonic[1], 1e-12 * ca.scattering_um2);
  // oblique: the azimuth of the plane of incidence only rotates the field (each harmonic by a phase);
  // "auto" would take |m| <= 8 here
  const std::string some = "[-2, -1, 0, 1, 2]";
  const fg::Result c = solveConfig(planeWaveConfig(kPillar, 30, 0, 20, "li_r", 2, some)),
                   d = solveConfig(planeWaveConfig(kPillar, 30, 40, 20, "li_r", 2, some));
  const auto& cc = c.cylindrical;
  const auto& cd = d.cylindrical;
  ASSERT_EQ(cc.harmonics, cd.harmonics);
  EXPECT_GT(cc.reflectedPerHarmonic[2], 1e-6 * cc.scattering_um2);  // m = 0 is excited
  for (std::size_t i = 0; i < cc.harmonics.size(); ++i) {
    EXPECT_NEAR(cc.reflectedPerHarmonic[i], cd.reflectedPerHarmonic[i], 1e-10 * cc.scattering_um2) << cc.harmonics[i];
    EXPECT_NEAR(cc.transmittedPerHarmonic[i], cd.transmittedPerHarmonic[i], 1e-10 * cc.scattering_um2) << cc.harmonics[i];
  }
  // rings of the background material scatter nothing
  const fg::Result z = solveConfig(planeWaveConfig(
      R"({"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["air", "air"]})", 30, 0, 0, "li_r", 2, some));
  EXPECT_LT(z.cylindrical.scattering_um2, 1e-28);
  // RK with the source term (augmented system) = particular RCWA solutions (Laurent, oblique)
  for (const std::string fac : {"laurent", "li_r"}) {
    const fg::Result p = solveConfig(planeWaveConfig(kPillar, 30, 0, 20, fac, 2, "[0, 1]")),
                     q = solveConfig(planeWaveConfig(kPillar, 30, 0, 20, fac, 2, "[0, 1]"), fg::SolverKind::RK);
    EXPECT_NEAR(q.cylindrical.scattering_um2 / p.cylindrical.scattering_um2, 1.0, 1e-8) << fac;
    EXPECT_NEAR(q.cylindrical.scatteredUp_um2 / p.cylindrical.scatteredUp_um2, 1.0, 1e-8) << fac;
  }
}

namespace {

/// Mie scattering efficiency Q_sca of a sphere (size parameter x, relative index m; Bohren-Huffman,
/// BHMIE: logarithmic derivative downward, Riccati-Bessel functions upward).
double mieQsca(double x, double m) {
  const int nstop = static_cast<int>(x + 4 * std::cbrt(x) + 2);
  const int nmx = static_cast<int>(std::max<double>(nstop, m * x)) + 16;
  std::vector<double> D(static_cast<std::size_t>(nmx) + 1, 0.0);
  for (int n = nmx; n > 0; --n) D[n - 1] = n / (m * x) - 1 / (D[n] + n / (m * x));
  double psi0 = std::cos(x), psi1 = std::sin(x), chi0 = -std::sin(x), chi1 = std::cos(x), q = 0;
  cplx xi1(psi1, -chi1);
  for (int n = 1; n <= nstop; ++n) {
    const double psi = (2 * n - 1) * psi1 / x - psi0, chi = (2 * n - 1) * chi1 / x - chi0;
    const cplx xi(psi, -chi);
    const double da = D[n] / m + n / x, db = m * D[n] + n / x;
    const cplx an = (da * psi - psi1) / (da * xi - xi1), bn = (db * psi - psi1) / (db * xi - xi1);
    q += (2 * n + 1) * (std::norm(an) + std::norm(bn));
    psi0 = psi1, psi1 = psi, chi0 = chi1, chi1 = chi, xi1 = xi;
  }
  return 2 * q / (x * x);
}

/// Sphere of radius a (k0 a = 3, lambda = 1 um) in air as N rings slices (radius at the slice centre).
double sphereSigma(int N, double kappa) {
  const double a = 3.0 / (2 * fg::PI), h = 2 * a / N;
  std::string layers;
  for (int i = 0; i < N; ++i) {
    const double yc = (i + 0.5) * h;
    char buf[200];
    std::snprintf(buf, sizeof buf, R"(%s{"type": "rings", "thickness_um": %.17g, "radii_um": [%.17g], "materials": ["sph", "air"]})",
                  i ? ", " : "", h, std::sqrt(a * a - (yc - a) * (yc - a)));
    layers += buf;
  }
  const fg::Result r = solveConfig(planeWaveConfig(layers, 0, 0, 0, "li_r", kappa, "\"auto\"", "air"));
  return r.cylindrical.scattering_um2;
}

}  // namespace

TEST(V21c, MieSphere) {
  // dielectric sphere n = 1.5, k0 a = 3, normal incidence, Li in r with the consistent source;
  // staircase of N slices: order p in N from kappa_max = 4 (N = 10, 20, 40), Richardson at
  // kappa_max = 8 (N = 20, 40), the rest in kappa_max as 1/kappa^2 (from N = 20)
  const double a = 3.0 / (2 * fg::PI), mie = mieQsca(3.0, 1.5) * fg::PI * a * a;
  const double s10 = sphereSigma(10, 4), s20 = sphereSigma(20, 4), s40 = sphereSigma(40, 4);
  const double t20 = sphereSigma(20, 8), t40 = sphereSigma(40, 8);
  const double p = std::log2((s10 - s20) / (s20 - s40));
  const double nInf = t40 - (t20 - t40) / (std::pow(2.0, p) - 1);
  const double est = nInf + (t20 - s20) / 3;  // kappa^-2: S(inf) - S(8) = (S(8) - S(4)) / 3
  std::printf("V21c Mie %.6f | kappa 4: %.6f %.6f %.6f (p = %.2f) | kappa 8: %.6f %.6f | extrapolated %.6f (%.1e)\n",
              mie, s10, s20, s40, p, t20, t40, est, est / mie - 1);
  EXPECT_NEAR(mie, 2.448002410791, 1e-9);  // reference implementation
  EXPECT_GT(p, 1.5);
  EXPECT_NEAR(est / mie, 1.0, 1e-3);
}

// Near field and volume integrals (Plan "M20: S-Matrix, Loeser und Nahfeld", "M20: Ausgaben"):
// focal plane of the beam = Gaussian, absorption inside the rings = flux balance, FFF reconstruction
// at a wall, optical theorem, gold sphere against Mie.
namespace {

fg::NearFieldResult nearFieldOf(const std::string& json, fg::SolverKind kind = fg::SolverKind::RCWA) {
  fg::io::RunConfig cfg = fg::io::parseConfig(json);
  cfg.solver.kind = kind;
  return fg::makeSolver(cfg.solver)->nearField(cfg.stack, cfg.scan.front(), cfg.field);
}

/// json with a "nearfield" block (inserted before the closing brace).
std::string withNearField(const std::string& json, const std::string& block) {
  const auto p = json.rfind('}');
  return json.substr(0, p) + R"(, "nearfield": {"enabled": true, )" + block + "}}";
}

std::string airBeam(const std::string& pol) {
  return R"({"cylindrical": {"spectrum": {"kappa_max": 2}},
    "materials": {"air": {"n": 1.0}}, "superstrate": "air", "substrate": "air",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 1.5, "focus_y_um": 0.5,
                  "polarization": ")" + pol + R"("},
    "layers": [{"type": "homogeneous", "thickness_um": 0.2, "material": "air"}]})";
}

/// Beam on a pillar of `pillar` (n or [n, k]) in air on glass.
std::string pillarBeam(const std::string& pillar, const std::string& fac, int kappa) {
  return R"({"cylindrical": {"spectrum": {"kappa_max": )" + std::to_string(kappa) + R"(}, "factorization": ")" + fac +
         R"("},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "p": {"n": )" + pillar + R"(}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.8, "polarization": "x"},
    "layers": [{"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["p", "air"]}]})";
}

/// Mie efficiencies (Q_ext, Q_sca) of a sphere with complex relative index m (Bohren-Huffman BHMIE).
std::pair<double, double> mieQ(double x, cplx m) {
  const int nstop = static_cast<int>(x + 4 * std::cbrt(x) + 2);
  const int nmx = static_cast<int>(std::max<double>(nstop, std::abs(m * x))) + 16;
  std::vector<cplx> D(static_cast<std::size_t>(nmx) + 1, 0.0);
  for (int n = nmx; n > 0; --n) D[n - 1] = double(n) / (m * x) - 1.0 / (D[n] + double(n) / (m * x));
  double psi0 = std::cos(x), psi1 = std::sin(x), chi0 = -std::sin(x), chi1 = std::cos(x), qe = 0, qs = 0;
  cplx xi1(psi1, -chi1);
  for (int n = 1; n <= nstop; ++n) {
    const double psi = (2 * n - 1) * psi1 / x - psi0, chi = (2 * n - 1) * chi1 / x - chi0;
    const cplx xi(psi, -chi);
    const cplx da = D[n] / m + double(n) / x, db = m * D[n] + double(n) / x;
    const cplx an = (da * psi - psi1) / (da * xi - xi1), bn = (db * psi - psi1) / (db * xi - xi1);
    qe += (2 * n + 1) * (an + bn).real();
    qs += (2 * n + 1) * (std::norm(an) + std::norm(bn));
    psi0 = psi1, psi1 = psi, chi0 = chi1, chi1 = chi, xi1 = xi;
  }
  return {2 * qe / (x * x), 2 * qs / (x * x)};
}

/// Sphere of radius a in air (N rings slices, radius at the slice centre), plane wave at normal incidence.
fg::CylindricalResult sphere(const std::string& n, double a, double lambda, int N, double kappa,
                             const std::string& fac = "li_r") {
  const double h = 2 * a / N;
  std::string layers;
  for (int i = 0; i < N; ++i) {
    const double yc = (i + 0.5) * h;
    char buf[200];
    std::snprintf(buf, sizeof buf, R"(%s{"type": "rings", "thickness_um": %.17g, "radii_um": [%.17g], "materials": ["sph", "air"]})",
                  i ? ", " : "", h, std::sqrt(a * a - (yc - a) * (yc - a)));
    layers += buf;
  }
  char head[400];
  std::snprintf(head, sizeof head, R"({"cylindrical": {"spectrum": {"kappa_max": %.17g}, "factorization": "%s", "optical_theorem": true},
    "materials": {"air": {"n": 1.0}, "sph": {"n": %s}}, "superstrate": "air", "substrate": "air",
    "incidence": {"type": "plane_wave", "lambda_um": %.17g, "phi_deg": 0, "psi_deg": 0, "delta_deg": 0}, "layers": [)",
                kappa, fac.c_str(), n.c_str(), lambda);
  return solveConfig(std::string(head) + layers + "]}").cylindrical;
}

}  // namespace

TEST(V21NearField, FocalPlaneOfTheBeamIsGaussian) {
  // all air, focus above the stack: the field in the focal plane is the prescribed beam, E_x = e^{-r^2/w0^2}
  // (x polarization), E_z = e^{-r^2/w0^2} (z), E_r = k0 r e^{-r^2/w0^2} (radial, the order-1 profile);
  // the aperture <= 1 cuts the spectrum at e^{-(k0 w0)^2/4} = 3e-10
  const double w0 = 1.5, k0 = 2 * fg::PI;
  for (const std::string pol : {"x", "z", "radial"}) {
    const fg::NearFieldResult nf = nearFieldOf(withNearField(
        airBeam(pol), R"("horizontal_y_um": 0.5, "nx": 13, "nz": 9, "ny": 3, "y_range_um": [0.5, 0.9], "extent_um": 2.5)"));
    const fg::FieldCut& c = nf.horizontal;
    ASSERT_EQ(c.E[0].size(), 13u * 9u);
    EXPECT_TRUE(nf.warnings.empty());
    double err = 0;
    for (std::size_t iz = 0; iz < 9; ++iz)
      for (std::size_t ix = 0; ix < 13; ++ix) {
        const double x = c.axis1_um[ix], z = c.axis2_um[iz], g = std::exp(-(x * x + z * z) / (w0 * w0));
        const cplx ex = pol == "z" ? 0.0 : pol == "x" ? g : k0 * x * g, ez = pol == "z" ? g : pol == "x" ? 0.0 : k0 * z * g;
        err = std::max({err, std::abs(c.E[0][iz * 13 + ix] - ex), std::abs(c.E[2][iz * 13 + ix] - ez)});
      }
    std::printf("V21 near field, focal plane %-6s: max |E - E_beam| = %.2e\n", pol.c_str(), err);
    EXPECT_LT(err, 1e-8) << pol;  // radial: amplitude up to 4
    // the vertical cut (z = 0) contains the focal plane as its first row
    for (std::size_t ix = 0; ix < 13; ++ix) EXPECT_LT(std::abs(nf.vertical.E[0][ix] - c.E[0][4 * 13 + ix]), 1e-12);
  }
}

TEST(V21NearField, AbsorptionInsideTheRingsEqualsTheFluxBalance) {
  // lossy pillar in air on glass: no guided modes, so A = 1 - R - T is the power absorbed in the rings,
  // here from the volume integral 1/2 int Im eps |E|^2 dV of the synthesized field. Laurent: the discrete
  // energy theorem holds, A_rings = A to quadrature accuracy (measured 4e-14 at kappa_max = 4). Li in r
  // (G_pm, Plan eq:cyl_G): A_rings - A is its discretization error and falls like kappa^-2 (2.3e-3,
  // 5.8e-4, 1.5e-4 at kappa_max = 2, 4, 8; the direct reconstruction only like 1/kappa: 7.5e-3, 3.5e-3, 1.9e-3)
  const fg::CylindricalResult l = solveConfig(pillarBeam("[1.5, 0.1]", "laurent", 4)).cylindrical;
  std::printf("V21 lossy pillar laurent: A = %.12f  A_rings = %.12f  (%.1e)\n", l.absorbed, l.absorbedRings,
              l.absorbedRings - l.absorbed);
  EXPECT_GT(l.absorbed, 0.1);
  EXPECT_NEAR(l.absorbedRings, l.absorbed, 1e-12);
  double d[2];
  for (int i = 0; i < 2; ++i) {
    const fg::CylindricalResult c = solveConfig(pillarBeam("[1.5, 0.1]", "li_r", 2 + 2 * i)).cylindrical;
    d[i] = c.absorbedRings - c.absorbed;
    std::printf("V21 lossy pillar li_r, kappa_max %d: A = %.12f  A_rings = %.12f  (%.1e)\n", 2 + 2 * i, c.absorbed,
                c.absorbedRings, d[i]);
  }
  EXPECT_LT(std::abs(d[1]), 4e-3 * l.absorbed);
  EXPECT_GT(d[0] / d[1], 3.0);  // second order
  // RK with the recorded heights = RCWA within the RK tolerance (kappa_max = 1: RK records are expensive)
  const fg::Result p = solveConfig(pillarBeam("[1.5, 0.1]", "laurent", 1)),
                   q = solveConfig(pillarBeam("[1.5, 0.1]", "laurent", 1), fg::SolverKind::RK);
  EXPECT_NEAR(q.cylindrical.absorbedRings, p.cylindrical.absorbedRings, 1e-8);
  // lossless rings: no volume integral
  EXPECT_EQ(solveConfig(pillarBeam("1.5", "li_r", 2)).cylindrical.absorbedRings, 0.0);
}

TEST(V21NearField, FFFReconstructionAtTheWall) {
  // at the wall r = R of the pillar D_r = eps E_r, E_phi and E_y are continuous (Plan sec:cylfact); the
  // FFF reconstruction (G_pm, Plan eq:cyl_G) keeps D_r continuous, the direct one has the Gibbs jump
  const double R = 0.4, d = 1e-6;
  const cplx epsIn = 1.5 * 1.5, epsOut = 1.0;
  double jump[2];
  int k = 0;
  for (const std::string rec : {"fff", "direct"}) {
    cplx D[2], Ez[2];
    for (int side = 0; side < 2; ++side) {
      // point (x, z = 0.05) at r = R -+ d: x = +extent of the vertical cut
      const double z = 0.05, r = R + (side ? d : -d), x = std::sqrt(r * r - z * z);
      char nfb[220];
      std::snprintf(nfb, sizeof nfb, R"("horizontal_y_um": 0.15, "nx": 2, "nz": 1, "ny": 1, "y_range_um": [0.15, 0.15], "extent_um": %.17g, "reconstruction": "%s", "vertical_z_um": %.17g)",
                    x, rec.c_str(), z);
      const fg::NearFieldResult nf = nearFieldOf(withNearField(pillarBeam("1.5", "li_r", 2), nfb));
      // E_r = (x E_x + z E_z) / r, E_phi = (x E_z - z E_x) / r
      const cplx Ex = nf.vertical.E[0][1], Ezz = nf.vertical.E[2][1];
      D[side] = (side ? epsOut : epsIn) * (x * Ex + z * Ezz) / r;
      Ez[side] = (x * Ezz - z * Ex) / r;  // E_phi
    }
    jump[k++] = std::abs(D[1] - D[0]) / std::abs(D[0]);
    std::printf("V21 wall (%s): |[D_r]| / |D_r| = %.2e, |[E_phi]| / |E_phi| = %.2e\n", rec.c_str(), jump[k - 1],
                std::abs(Ez[1] - Ez[0]) / std::abs(Ez[0]));
  }
  EXPECT_LT(jump[0], 1e-3);
  EXPECT_LT(10 * jump[0], jump[1]);
}

TEST(V21c, OpticalTheoremAndGoldSphere) {
  // dielectric sphere (V21c): sigma_abs = 0, the optical theorem (volume form, forward amplitude) equals
  // sigma_sca to the discretization error of the staircase
  const double a = 3.0 / (2 * fg::PI);
  const fg::CylindricalResult d = sphere("1.5", a, 1.0, 20, 4);
  std::printf("V21c dielectric: sigma_sca = %.6f  sigma_abs = %.1e  optical theorem %.6f (%.1e)\n", d.scattering_um2,
              d.absorption_um2, d.extinctionTheorem_um2, d.extinctionTheorem_um2 / d.scattering_um2 - 1);
  EXPECT_EQ(d.absorption_um2, 0.0);
  EXPECT_EQ(d.extinction_um2, d.scattering_um2);
  // without "optical_theorem" (lossless rings): no volume integrals, the theorem is NaN, sigma_sca unchanged
  const fg::CylindricalResult q = solveConfig(planeWaveConfig(kPillar, 0, 0, 0)).cylindrical;
  EXPECT_TRUE(std::isnan(q.extinctionTheorem_um2));
  EXPECT_EQ(q.extinction_um2, q.scattering_um2);
  EXPECT_NEAR(d.extinctionTheorem_um2 / d.scattering_um2, 1.0, 1e-2);
  // gold sphere (a = 50 nm, lambda = 550 nm, n = 0.43 + 2.45 i): k_max = kappa max Re n = kappa, so k_max a
  // needs kappa_max = 32 (k0 a = 0.57). Laurent conserves energy in the discrete sense: the optical theorem
  // equals sigma_sca + sigma_abs to rounding. Li in r, qualitative against Mie (Plan V21c): the staircase
  // of a metal sphere converges slowly in N and kappa_max (docs/validation.md)
  const cplx nAu(0.43, 2.45);
  const double aAu = 0.05, lam = 0.55, x = 2 * fg::PI * aAu / lam, geo = fg::PI * aAu * aAu;
  const auto [qe, qs] = mieQ(x, nAu);
  const fg::CylindricalResult l = sphere("[0.43, 2.45]", aAu, lam, 10, 16, "laurent");
  std::printf("V21c gold, Laurent: sigma_ext %.10e  theorem %.10e (%.1e)\n", l.extinction_um2, l.extinctionTheorem_um2,
              l.extinctionTheorem_um2 / l.extinction_um2 - 1);
  EXPECT_NEAR(l.extinctionTheorem_um2 / l.extinction_um2, 1.0, 1e-8);
  const fg::CylindricalResult g = sphere("[0.43, 2.45]", aAu, lam, 20, 32);
  std::printf("V21c gold: sigma_sca %.6e (Mie %.6e)  sigma_abs %.6e (Mie %.6e)  sigma_ext %.6e (Mie %.6e, theorem %.6e)\n",
              g.scattering_um2, qs * geo, g.absorption_um2, (qe - qs) * geo, g.extinction_um2, qe * geo,
              g.extinctionTheorem_um2);
  EXPECT_NEAR(qe * geo, 2.069503e-2, 1e-7);  // reference implementation (complex m)
  EXPECT_NEAR(g.scattering_um2 / (qs * geo), 1.0, 0.25);
  EXPECT_NEAR(g.absorption_um2 / ((qe - qs) * geo), 1.0, 0.25);
  EXPECT_NEAR(g.extinction_um2 / (qe * geo), 1.0, 0.05);
  EXPECT_NEAR(g.extinctionTheorem_um2 / g.extinction_um2, 1.0, 0.05);
}

// V21e (Plan "Tests und Abnahme (M20)"): a Gaussian beam in the scattered-field formulation. The background
// (every ring layer replaced by its outer material) lies on the Hankel grid, the source of a ring layer is
// sigma = (M - M_bg) F_bg (Plan sec:cylexc, "Konsistenz"); scattered field + background = total field: RCWA to
// 1e-10, RK to its tolerance, Laurent and Li in r, x (m = +-1) and radial (m = 0) polarization.
namespace {

/// Two ring layers (outer materials air and film) between films, lossy pillar.
std::string ringStack(const std::string& fac, const std::string& pol, int kappa, bool scattered,
                      const std::string& spectrum = "") {
  return R"({)" + std::string(scattered ? R"("scattered_field": true, )" : "") +
         R"("cylindrical": {"spectrum": {"kappa_max": )" + std::to_string(kappa) + spectrum + R"(}, "factorization": ")" + fac +
         R"("},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "film": {"n": 1.38}, "p": {"n": [1.7, 0.05]}, "q": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.9, "polarization": ")" +
         pol + R"("},
    "layers": [{"type": "homogeneous", "thickness_um": 0.15, "material": "film"},
               {"type": "rings", "thickness_um": 0.2, "radii_um": [0.2, 0.35], "materials": ["q", "p", "film"]},
               {"type": "rings", "thickness_um": 0.25, "radii_um": [0.3], "materials": ["p", "air"]},
               {"type": "homogeneous", "thickness_um": 0.1, "material": "air"}]})";
}

}  // namespace

TEST(V21e, BeamScatteredFieldEqualsTotalField) {
  for (const std::string fac : {"laurent", "li_r"})
    for (const std::string pol : {"x", "radial"}) {
      const fg::CylindricalResult t = solveConfig(ringStack(fac, pol, 2, false)).cylindrical,
                                  s = solveConfig(ringStack(fac, pol, 2, true)).cylindrical;
      std::printf("V21e %-7s %-6s: R = %.12f (%.1e)  T = %.12f (%.1e)  A_rings = %.12f (%.1e)\n", fac.c_str(),
                  pol.c_str(), t.reflected, s.reflected - t.reflected, t.transmitted, s.transmitted - t.transmitted,
                  t.absorbedRings, s.absorbedRings - t.absorbedRings);
      EXPECT_GT(t.absorbedRings, 1e-3);
      EXPECT_NEAR(s.reflected, t.reflected, 1e-10) << fac << " " << pol;
      EXPECT_NEAR(s.transmitted, t.transmitted, 1e-10) << fac << " " << pol;
      EXPECT_NEAR(s.absorbedRings, t.absorbedRings, 1e-10) << fac << " " << pol;
      ASSERT_EQ(s.reflectedPerHarmonic.size(), t.reflectedPerHarmonic.size());
      for (std::size_t i = 0; i < t.reflectedPerHarmonic.size(); ++i) {
        EXPECT_NEAR(s.reflectedPerHarmonic[i], t.reflectedPerHarmonic[i], 1e-10);
        EXPECT_NEAR(s.transmittedPerHarmonic[i], t.transmittedPerHarmonic[i], 1e-10);
      }
    }
  // near field (total field, both cuts through the rings and the half spaces)
  const std::string nfb =
      R"("horizontal_y_um": 0.3, "nx": 9, "nz": 5, "ny": 7, "y_range_um": [-0.1, 0.8], "extent_um": 0.6)";
  const fg::NearFieldResult a = nearFieldOf(withNearField(ringStack("li_r", "x", 2, false), nfb)),
                            b = nearFieldOf(withNearField(ringStack("li_r", "x", 2, true), nfb));
  double err = 0, ref = 0;
  for (auto cut : {&fg::NearFieldResult::horizontal, &fg::NearFieldResult::vertical})
    for (int c = 0; c < 3; ++c)
      for (std::size_t i = 0; i < (a.*cut).E[c].size(); ++i) {
        err = std::max({err, std::abs((a.*cut).E[c][i] - (b.*cut).E[c][i]), std::abs((a.*cut).H[c][i] - (b.*cut).H[c][i])});
        ref = std::max(ref, std::abs((a.*cut).E[c][i]));
      }
  std::printf("V21e near field: max |F_s + F_bg - F| = %.1e (max |E| = %.2f)\n", err, ref);
  EXPECT_LT(err, 1e-10 * ref);
  // RK with the source integrated along (kappa_max = 1 on the real path: RK records are expensive): within
  // the RK tolerance (measured 5e-9)
  const std::string real = R"(, "deform": false)";
  const fg::CylindricalResult p = solveConfig(ringStack("laurent", "x", 1, false, real)).cylindrical,
                              q = solveConfig(ringStack("laurent", "x", 1, true, real), fg::SolverKind::RK).cylindrical;
  EXPECT_NEAR(q.reflected, p.reflected, 1e-8);
  EXPECT_NEAR(q.transmitted, p.transmitted, 1e-8);
  EXPECT_NEAR(q.absorbedRings, p.absorbedRings, 1e-8);
}

// rings_profile (Plan "M20: Schichttypen"): radii as functions of y, solved as an RCWA staircase with the radii
// at the slice centres. The staircase equals the explicit stack of rings layers bit for bit; a constant profile
// is the rings layer; the staircase error of a cone falls like 1/K^2 (midpoint radii).
namespace {

std::string coneBeam(const std::string& layers) {
  return R"({"cylindrical": {"spectrum": {"kappa_max": 2}, "factorization": "li_r"},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "p": {"n": 1.6}, "q": {"n": 1.3}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.9, "polarization": "x"},
    "layers": [)" + layers + "]}";
}

std::string coneProfile(int K) {
  return R"({"type": "rings_profile", "thickness_um": 0.4, "slices": )" + std::to_string(K) +
         R"(, "radii": [{"table": {"y_um": [0, 0.4], "r_um": [0.1, 0.25]}}, {"cone": [0.45, 0.25]}],
            "materials": ["q", "p", "air"]})";
}

}  // namespace

TEST(V21RingsProfile, StaircaseOfRingsLayers) {
  // explicit staircase: radii at the slice centres (the inner ring, a table, meets the cone only at the top)
  const int K = 6;
  std::string layers;
  for (int k = 0; k < K; ++k) {
    const double y = (k + 0.5) * (0.4 / K), r1 = 0.1 + (0.25 - 0.1) * y / 0.4, r2 = 0.45 + (0.25 - 0.45) * y / 0.4;
    char buf[200];
    std::snprintf(buf, sizeof buf, R"(%s{"type": "rings", "thickness_um": %.17g, "radii_um": [%.17g, %.17g], "materials": ["q", "p", "air"]})",
                  k ? ", " : "", 0.4 / K, r1, r2);
    layers += buf;
  }
  const fg::CylindricalResult a = solveConfig(coneBeam(layers)).cylindrical, b = solveConfig(coneBeam(coneProfile(K))).cylindrical;
  EXPECT_EQ(a.reflected, b.reflected);
  EXPECT_EQ(a.transmitted, b.transmitted);
  // constant radii = rings layer (one slice: the same layer; three: three equal layers)
  const std::string rings = R"({"type": "rings", "thickness_um": 0.3, "radii_um": [0.2, 0.4], "materials": ["q", "p", "air"]})";
  const fg::CylindricalResult r = solveConfig(coneBeam(rings)).cylindrical;
  for (const int k : {1, 3}) {
    const fg::CylindricalResult c = solveConfig(coneBeam(
        R"({"type": "rings_profile", "thickness_um": 0.3, "slices": )" + std::to_string(k) +
        R"(, "radii": [0.2, 0.4], "materials": ["q", "p", "air"]})")).cylindrical;
    EXPECT_NEAR(c.reflected, r.reflected, k == 1 ? 0.0 : 1e-12);
    EXPECT_NEAR(c.transmitted, r.transmitted, k == 1 ? 0.0 : 1e-12);
  }
  // staircase error of the cone: second order in 1/K
  double R[3];
  for (int i = 0; i < 3; ++i) R[i] = solveConfig(coneBeam(coneProfile(3 << i))).cylindrical.reflected;
  std::printf("V21 rings_profile cone: R(3, 6, 12) = %.10f %.10f %.10f, ratio %.2f\n", R[0], R[1], R[2],
              (R[0] - R[1]) / (R[1] - R[2]));
  EXPECT_GT(std::abs(R[0] - R[1]), 1e-6);
  EXPECT_NEAR((R[0] - R[1]) / (R[1] - R[2]), 4.0, 1.0);
}

TEST(V21RingsProfile, SphereAndCapAndErrors) {
  // a sphere (centre y = a, radius a, thickness 2a) as rings_profile = the explicit staircase of V21c
  const double a = 3.0 / (2 * fg::PI);
  char sph[300];
  std::snprintf(sph, sizeof sph, R"({"type": "rings_profile", "thickness_um": %.17g, "slices": 10,
    "radii": [{"sphere": {"center_y_um": %.17g, "radius_um": %.17g}}], "materials": ["sph", "air"]})", 2 * a, a, a);
  const fg::Result p = solveConfig(planeWaveConfig(sph, 0, 0, 0, "li_r", 2, "\"auto\"", "air"));
  EXPECT_EQ(p.cylindrical.scattering_um2, sphereSigma(10, 2));
  // a cap: the upper slices have no ring left and become homogeneous layers of the outer material
  const std::string cap = R"({"type": "rings_profile", "thickness_um": 0.5, "slices": 10,
    "radii": [{"sphere": {"center_y_um": 0.0, "radius_um": 0.3}}], "materials": ["sph", "air"]})";
  const std::string capLow = R"({"type": "rings_profile", "thickness_um": 0.3, "slices": 6,
    "radii": [{"sphere": {"center_y_um": 0.0, "radius_um": 0.3}}], "materials": ["sph", "air"]},
    {"type": "homogeneous", "thickness_um": 0.2, "material": "air"})";
  const fg::Result c1 = solveConfig(planeWaveConfig(cap, 0, 0, 0, "li_r", 2)),
                   c2 = solveConfig(planeWaveConfig(capLow, 0, 0, 0, "li_r", 2));
  EXPECT_GT(c1.cylindrical.scattering_um2, 1e-3);
  EXPECT_NEAR(c1.cylindrical.scattering_um2, c2.cylindrical.scattering_um2, 1e-12);
  // crossing radii and wrong material counts are errors
  auto err = [](const std::string& layer) {
    try {
      solveConfig(planeWaveConfig(layer, 0, 0, 0));
    } catch (const fg::Error& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  EXPECT_NE(err(R"({"type": "rings_profile", "thickness_um": 0.3, "radii": [{"cone": [0.1, 0.4]}, 0.3],
                    "materials": ["sph", "film", "air"]})").find("must not cross"), std::string::npos);
  EXPECT_NE(err(R"({"type": "rings_profile", "thickness_um": 0.3, "radii": [0.3], "materials": ["sph"]})")
                .find("n + 1 materials"), std::string::npos);
  EXPECT_NE(err(R"({"type": "rings_profile", "thickness_um": 0.3, "radii": [{"cube": 1}], "materials": ["sph", "air"]})")
                .find("rings_profile"), std::string::npos);
}

// Angular spectrum and far field (Plan "M20: Ausgaben"): per harmonic and side dP/dtheta at the propagating
// nodes, split exactly into s ((e_+ + e_-)/2, TE) and p ((e_+ - e_-)/2, TM); with the quadrature weight in theta
// the sum is the power of the harmonic. Far field dP/dOmega on (theta nodes) x (phi): its phi average times
// 2 pi sin(theta) is the sum of the angular spectra; a small sphere gives the Mie differential cross sections.
namespace {

/// Mie amplitude functions S1, S2 at the scattering angle Theta (Bohren-Huffman BHMIE, real index m).
std::pair<cplx, cplx> mieS(double x, double m, double Theta) {
  const int nstop = static_cast<int>(x + 4 * std::cbrt(x) + 2);
  const int nmx = static_cast<int>(std::max<double>(nstop, m * x)) + 16;
  std::vector<double> D(static_cast<std::size_t>(nmx) + 1, 0.0);
  for (int n = nmx; n > 0; --n) D[n - 1] = n / (m * x) - 1 / (D[n] + n / (m * x));
  double psi0 = std::cos(x), psi1 = std::sin(x), chi0 = -std::sin(x), chi1 = std::cos(x);
  cplx xi1(psi1, -chi1), S1 = 0.0, S2 = 0.0;
  const double mu = std::cos(Theta);
  double pi0 = 0, pi1 = 1;
  for (int n = 1; n <= nstop; ++n) {
    const double psi = (2 * n - 1) * psi1 / x - psi0, chi = (2 * n - 1) * chi1 / x - chi0;
    const cplx xi(psi, -chi);
    const double da = D[n] / m + n / x, db = m * D[n] + n / x;
    const cplx an = (da * psi - psi1) / (da * xi - xi1), bn = (db * psi - psi1) / (db * xi - xi1);
    const double tau = n * mu * pi1 - (n + 1) * pi0, f = (2.0 * n + 1) / (n * (n + 1.0));
    S1 += f * (an * pi1 + bn * tau);
    S2 += f * (an * tau + bn * pi1);
    const double pi2 = ((2 * n + 1) * mu * pi1 - (n + 1) * pi0) / n;
    pi0 = pi1, pi1 = pi2;
    psi0 = psi1, psi1 = psi, chi0 = chi1, chi1 = chi, xi1 = xi;
  }
  return {S1, S2};
}

}  // namespace

TEST(V21FarField, AngularSpectrumAndFarFieldOfABeam) {
  // lossy pillar on glass, x polarization: per harmonic and side sum w (s + p) = C_Rm / C_Tm
  std::string json = pillarBeam("[1.5, 0.05]", "li_r", 2);
  json.insert(json.find("\"factorization\""), R"("far_field": {"n_phi": 24}, )");
  const fg::CylindricalResult c = solveConfig(json).cylindrical;
  ASSERT_EQ(c.angularSpectrum.size(), 4u);  // m = -1, 1; up, down
  double worst = 0;
  for (const auto& a : c.angularSpectrum) {
    ASSERT_FALSE(a.theta_deg.empty());
    EXPECT_LT(a.theta_deg.back(), 90.0);
    double P = 0;
    for (std::size_t j = 0; j < a.theta_deg.size(); ++j) P += a.weight[j] * (a.s[j] + a.p[j]);
    const auto i = static_cast<std::size_t>(std::find(c.harmonics.begin(), c.harmonics.end(), a.m) - c.harmonics.begin());
    const double ref = a.side > 0 ? c.reflectedPerHarmonic.at(i) : c.transmittedPerHarmonic.at(i);
    worst = std::max(worst, std::abs(P - ref));
  }
  std::printf("V21 angular spectrum: max |sum w (s + p) - C_m| = %.1e\n", worst);
  EXPECT_LT(worst, 1e-12);
  // far field: up and down; phi average x 2 pi sin(theta) = sum over m of the angular spectra
  ASSERT_EQ(c.farField.size(), 2u);
  for (const auto& f : c.farField) {
    const std::size_t nt = f.theta_deg.size(), np = f.phi_deg.size();
    ASSERT_EQ(np, 24u);
    ASSERT_EQ(f.s.size(), nt * np);
    double err = 0, mx = 0;
    for (std::size_t j = 0; j < nt; ++j) {
      double avg = 0, sum = 0;
      for (std::size_t l = 0; l < np; ++l) avg += (f.s[j * np + l] + f.p[j * np + l]) / np;
      for (const auto& a : c.angularSpectrum)
        if (a.side == f.side) sum += a.s[j] + a.p[j];
      const double th = f.theta_deg[j] * fg::PI / 180;
      err = std::max(err, std::abs(2 * fg::PI * std::sin(th) * avg - sum));
      mx = std::max(mx, sum);
    }
    EXPECT_LT(err, 1e-12 * std::max(1.0, mx)) << "side " << f.side;
    // x polarization: mirror symmetric about the x-y plane, dP(phi) = dP(-phi)
    for (std::size_t j = 0; j < nt; ++j)
      for (std::size_t l = 1; l < np; ++l)
        EXPECT_NEAR(f.s[j * np + l], f.s[j * np + (np - l)], 1e-12 * std::max(1.0, mx));
  }
}

TEST(V21FarField, SmallSphereGivesMieDifferentialCrossSections) {
  // dielectric sphere (n = 1.5, k0 a = 1) in air, plane wave along -y with E along z (delta = 0): in the plane
  // phi = 0 (x-y, perpendicular to E) dsigma/dOmega = |S1|^2 / k^2, in phi = 90 deg (y-z) |S2|^2 / k^2, at the
  // scattering angle pi - theta above and theta below. The deviation is the truncation in k_max (measured,
  // relative to the largest value: 2.2e-2 at kappa_max = 4, 4.7e-3 at 8; more slices make no difference)
  const double a = 1.0 / (2 * fg::PI), k = 2 * fg::PI;
  auto deviation = [&](double kappa, double& total, double& sigma) {
    char sph[300];
    std::snprintf(sph, sizeof sph, R"({"type": "rings_profile", "thickness_um": %.17g, "slices": 16,
      "radii": [{"sphere": {"center_y_um": %.17g, "radius_um": %.17g}}], "materials": ["sph", "air"]})", 2 * a, a, a);
    std::string json = planeWaveConfig(sph, 0, 0, 0, "li_r", kappa, "\"auto\"", "air");
    json.insert(json.find("\"factorization\""), R"("far_field": {"n_phi": 4}, )");
    const fg::CylindricalResult c = solveConfig(json).cylindrical;
    EXPECT_EQ(c.farField.size(), 2u);
    double worst = 0;
    for (const auto& f : c.farField) {
      double ref = 0, err = 0;
      for (std::size_t j = 0; j < f.theta_deg.size(); ++j) {
        const double th = f.theta_deg[j] * fg::PI / 180, Theta = f.side > 0 ? fg::PI - th : th;
        const auto [S1, S2] = mieS(k * a, 1.5, Theta);
        const double r0 = std::norm(S1) / (k * k), r90 = std::norm(S2) / (k * k);
        ref = std::max({ref, r0, r90});
        err = std::max({err, std::abs(f.s[j * 4] + f.p[j * 4] - r0), std::abs(f.s[j * 4 + 1] + f.p[j * 4 + 1] - r90)});
      }
      worst = std::max(worst, err / ref);
    }
    total = 0;
    for (const auto& s : c.angularSpectrum)
      for (std::size_t j = 0; j < s.theta_deg.size(); ++j) total += s.weight[j] * (s.s[j] + s.p[j]);
    sigma = c.scattering_um2;
    return worst;
  };
  double total = 0, sigma = 0;
  const double d4 = deviation(4, total, sigma), d8 = deviation(8, total, sigma);
  std::printf("V21 far field of a small sphere: deviation from Mie %.2e (kappa_max 4), %.2e (8); sigma %.6e = C_Ssca %.6e\n",
              d4, d8, total, sigma);
  EXPECT_NEAR(total, sigma, 1e-12 * sigma);
  EXPECT_LT(d8, 0.02);
  EXPECT_LT(d8, 0.7 * d4);
}

TEST(V21FarField, CsvFiles) {
  // fgsolve writes angular_spectrum.csv and far_field.csv with cylindrical.far_field (one row per node and
  // harmonic / per node and azimuth)
  std::string json = pillarBeam("1.5", "li_r", 2);
  json.insert(json.find("\"factorization\""), R"("far_field": {"n_phi": 6}, )");
  const fg::io::RunConfig cfg = fg::io::parseConfig(json);
  const std::vector<fg::Result> res = fg::solveScan(cfg.solver, cfg.stack, cfg.scan);
  const auto dir = std::filesystem::temp_directory_path() / "fg_v21_far_field";
  const auto files = fg::io::writeCylindricalSpectra(dir.string(), {cfg}, {res});
  ASSERT_EQ(files.size(), 2u);
  auto lines = [](const std::string& f) {
    std::ifstream in(f);
    std::string l;
    std::size_t n = 0;
    while (std::getline(in, l)) ++n;
    return n;
  };
  std::size_t na = 0, nf = 0;
  for (const auto& a : res[0].cylindrical.angularSpectrum) na += a.theta_deg.size();
  for (const auto& f : res[0].cylindrical.farField) nf += f.theta_deg.size() * f.phi_deg.size();
  EXPECT_EQ(lines(files[0]), na + 1);
  EXPECT_EQ(lines(files[1]), nf + 1);
  // without far_field nothing is written
  const fg::io::RunConfig plain = fg::io::parseConfig(pillarBeam("1.5", "li_r", 2));
  EXPECT_TRUE(fg::io::writeCylindricalSpectra(dir.string(), {plain}, {fg::solveScan(plain.solver, plain.stack, plain.scan)}).empty());
  std::filesystem::remove_all(dir);
}

// V21f (Plan "Tests und Abnahme (M20)"): symmetry and integration path, aperture (r = 0.25 um) in a metal film
// (0.1 um) on glass, beam at 633 nm. The deformed path only moves the nodes, so the results depend on its depth
// within the discretization error of the grid (Befund v2.31); the wrong side (delta < 0, past the plasmon pole)
// jumps.
namespace {

std::string aperture(const std::string& n, const std::string& deform, double dk = 1) {
  char head[200];
  std::snprintf(head, sizeof head, R"({"cylindrical": {"spectrum": {"kappa_max": 2, "dk_factor": %.17g, "deform": %s},)",
                dk, deform.c_str());
  return std::string(head) + R"( "factorization": "li_r"},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "metal": {"n": )" + n + R"(}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 0.633, "waist_um": 0.4, "focus_y_um": 0.1, "polarization": "x"},
    "layers": [{"type": "rings", "thickness_um": 0.1, "radii_um": [0.25], "materials": ["air", "metal"]}]})";
}

}  // namespace

TEST(V21f, SymmetryAndIntegrationPath) {
  const std::string ag = "[0.135, 3.99]";  // silver: eps = -15.9 + 1.08 i, plasmon pole close to the real axis
  const fg::CylindricalResult a = solveConfig(aperture(ag, R"({"delta": 0.05})")).cylindrical,
                              b = solveConfig(aperture(ag, R"({"delta": 0.1})")).cylindrical;
  // m <-> -m (x polarization: m = -1 and m = +1 carry the same power)
  EXPECT_NEAR(a.reflectedPerHarmonic[0], a.reflectedPerHarmonic[1], 1e-12);
  EXPECT_NEAR(a.transmittedPerHarmonic[0], a.transmittedPerHarmonic[1], 1e-12);
  // depth of the path: within the discretization error (dk_factor 1: 9.4e-5 / 4.6e-5, while R and T change by
  // 1e-4 / 2e-5 from dk_factor 1 to 2; V21h)
  std::printf("V21f silver: R = %.10f, T = %.10f; delta 0.05 -> 0.1: %.1e, %.1e\n", a.reflected, a.transmitted,
              b.reflected - a.reflected, b.transmitted - a.transmitted);
  EXPECT_LT(std::abs(b.reflected - a.reflected), 1e-4);
  EXPECT_LT(std::abs(b.transmitted - a.transmitted), 1e-3);
  // the wrong side (delta < 0, only through the C++ API) passes the plasmon pole: a jump far beyond that
  fg::io::RunConfig cfg = fg::io::parseConfig(aperture(ag, R"({"delta": 0.05})"));
  cfg.stack.cylindrical.deltaRel = -0.05;
  const fg::CylindricalResult w = fg::makeSolver(cfg.solver)->solve(cfg.stack, cfg.scan.front()).cylindrical;
  std::printf("V21f silver, wrong side: R = %.10f, T = %.10f\n", w.reflected, w.transmitted);
  EXPECT_GT(std::abs(w.transmitted - a.transmitted) + std::abs(w.reflected - a.reflected), 1e-2);
  // lossy metal (eps = -8.75 + 3 i, pole far from the axis): the real path is accurate and equals the deformed
  // one to quadrature accuracy (measured 2.3e-8 for |dR| + |dT|); with decreasing loss the pole approaches the
  // axis and the real path departs from the deformed one (2.7e-6, 1.2e-4; before the graded panels of v2.33
  // around the plasmon 2.3e-4, 4.6e-3)
  double d[3];
  int i = 0;
  for (const std::string n : {"[0.5, 3.0]", "[0.3, 3.5]", "[0.15, 3.9]"}) {
    const fg::CylindricalResult r = solveConfig(aperture(n, "false")).cylindrical,
                                p = solveConfig(aperture(n, R"({"delta": 0.05})")).cylindrical;
    d[i] = std::abs(r.reflected - p.reflected) + std::abs(r.transmitted - p.transmitted);
    std::printf("V21f metal n = %s: |real - deformed| = %.1e\n", n.c_str(), d[i]);
    ++i;
  }
  EXPECT_LT(d[0], 1e-7);
  EXPECT_LT(d[0], d[1]);
  EXPECT_LT(d[1], d[2]);
  // the configuration accepts only delta > 0
  EXPECT_THROW(fg::io::parseConfig(aperture(ag, R"({"delta": -0.05})")), fg::Error);
}

TEST(V21h, LeakyPlasmonBetweenTheBranchPoints) {
  // silver film between air and glass: the plasmon of the silver/air interface (k = 1.033 + 0.0023 i) lies
  // between the branch points 1 and 1.5 and leaks into the glass. Before v2.33, T oscillated in dk_factor
  // (0.416, 0.453, 0.537, 0.504, 0.477, 0.495 for 1 ... 6); graded panels around Re k_p (Plan sec:cylgrid)
  // make it converge: fgsolve gives T = 0.4940938 (dk_factor 1), 0.4940715 (2), 0.4940711 (3)
  const std::string ag = "[0.135, 3.99]";
  const fg::CylindricalResult h = solveConfig(aperture(ag, "\"auto\"", 0.5)).cylindrical,
                              f = solveConfig(aperture(ag, "\"auto\"", 1.0)).cylindrical;
  std::printf("V21h silver/air: dk 0.5: R = %.10f, T = %.10f; dk 1: R = %.10f, T = %.10f\n", h.reflected,
              h.transmitted, f.reflected, f.transmitted);
  EXPECT_NEAR(f.transmitted, 0.4940711, 1e-4);
  EXPECT_NEAR(f.reflected, 0.3504524, 2e-4);
  EXPECT_NEAR(h.transmitted, f.transmitted, 1e-3);
  EXPECT_NEAR(h.reflected, f.reflected, 2e-3);
}

