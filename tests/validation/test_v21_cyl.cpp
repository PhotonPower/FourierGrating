// V21 (M20, Plan section "Tests und Abnahme (M20)"): rotationally symmetric structures.
// V21b: guided modes of a step-index fiber are eigenvalues +-i beta of Plan eq:cyl_M.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <functional>

#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"

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
