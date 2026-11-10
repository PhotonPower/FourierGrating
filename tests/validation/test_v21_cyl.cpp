// V21 (M20, Plan section "Tests und Abnahme (M20)"): rotationally symmetric structures.
// V21b: guided modes of a step-index fiber are eigenvalues +-i beta of Plan eq:cyl_M.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
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
    return std::cyl_bessel_j(double(std::abs(n)), x) * ((n < 0 && n % 2) ? -1 : 1);
  };
  auto K = [](int n, double x) { return std::cyl_bessel_k(double(std::abs(n)), x); };
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
