// V21 (M20, Plan section "Tests und Abnahme (M20)"): rotationally symmetric structures.
// V21b: guided modes of a step-index fiber are eigenvalues +-i beta of Plan eq:cyl_M.
#include <gtest/gtest.h>

#include "../support/std_bessel.hpp"

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
