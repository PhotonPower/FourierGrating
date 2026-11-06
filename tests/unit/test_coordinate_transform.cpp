// Unit tests of the coordinate transformations (Plan section "M17", eq:ct_map, eq:ct_deriv, eq:asr,
// eq:ct_inc, eq:ct_rayleigh).
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "fg/fourier/contour.hpp"
#include "fg/geometry/coordinate_transform.hpp"
#include "fg/physics/modes.hpp"
#include "fg/physics/propagation_matrix.hpp"

using fg::cplx;
using fg::CoordinateMap;
using fg::MatC;
using fg::VecC;

TEST(CoordinateMap, AsrFixedPointsDerivativeInverse) {
  const CoordinateMap m(1.3, {0.95, 0.2, 1.3 + 0.6}, 0.8);  // reduced: 0.2, 0.6, 0.95
  ASSERT_EQ(m.breakpoints().size(), 3u);
  for (double b : m.breakpoints()) {
    EXPECT_NEAR(m.f(b), b, 1e-15);
    EXPECT_NEAR(m.df(b), 1.0 - 0.8, 1e-15);  // f' continuous, 1 - eta at every breakpoint
    EXPECT_NEAR(m.f(b + 2.6), b + 2.6, 1e-14);
  }
  for (double X = -2.0; X < 3.0; X += 0.0137) {
    EXPECT_NEAR(m.f(X + 1.3), m.f(X) + 1.3, 1e-14);  // Plan eq:ct_map
    EXPECT_NEAR(m.inverse(m.f(X)), X, 1e-14);
    const double h = 1e-6;
    EXPECT_NEAR((m.f(X + h) - m.f(X - h)) / (2 * h), m.df(X), 1e-8);
    EXPECT_GT(m.df(X), 0.0);
  }
}

// One breakpoint: 1/f' = 1/(1 - eta cos(K (X - x0))) has the coefficients r^|p| e^{-i p K x0} / sqrt(1 - eta^2),
// r = (1 - sqrt(1 - eta^2)) / eta.
TEST(CoordinateMap, InverseDerivativeCoefficientsAnalytic) {
  const double d = 0.9, x0 = 0.31, eta = 0.8, K = 2 * fg::PI / d;
  const CoordinateMap m(d, {x0}, eta);
  const auto c = m.inverseDerivativeCoefficients(12);
  const double s = std::sqrt(1 - eta * eta), r = (1 - s) / eta;
  for (int p = -12; p <= 12; ++p)
    EXPECT_NEAR(std::abs(c[static_cast<std::size_t>(p + 12)] - std::pow(r, std::abs(p)) / s * std::exp(-fg::I * (p * K * x0))),
                0.0, 1e-14)
        << p;
}

// Several breakpoints: against a fine DFT of 1/f'.
TEST(CoordinateMap, InverseDerivativeCoefficientsDFT) {
  const double d = 1.0;
  const CoordinateMap m(d, {0.0, 0.45, 0.7}, 0.6);
  const auto c = m.inverseDerivativeCoefficients(10);
  const int M = 1 << 16;
  for (int p = -10; p <= 10; ++p) {
    cplx ref = 0;
    for (int k = 0; k < M; ++k) ref += std::exp(-fg::I * (2 * fg::PI * p * k / M)) / m.df(d * k / M);
    ref /= static_cast<double>(M);
    EXPECT_NEAR(std::abs(c[static_cast<std::size_t>(p + 10)] - ref), 0.0, 1e-9) << p;
  }
}

// C D = 1 on the central orders (eq:ct_rayleigh is the inverse of eq:ct_inc for N -> infinity).
TEST(CoordinateMap, PhysicalConversionsAreInverse) {
  const double d = 0.8, alpha0 = 3.1;
  for (const auto& br : {std::vector<double>{0.25}, std::vector<double>{0.1, 0.5}}) {
    const CoordinateMap m(d, br, 0.5);
    const int N = 60, n0 = 4;
    const MatC C = m.toPhysical(alpha0, N), D = m.fromPhysical(alpha0, N);
    const MatC CD = C * D;
    double err = 0;
    for (int a = -n0; a <= n0; ++a)
      for (int b = -n0; b <= n0; ++b) err = std::max(err, std::abs(CD(a + N, b + N) - (a == b ? 1.0 : 0.0)));
    EXPECT_LT(err, br.size() == 1 ? 1e-13 : 1e-8);
  }
}

// The plane wave e^{i alpha_0 f(X)} expanded with eq:ct_inc: its X coefficients synthesize f.
TEST(CoordinateMap, IncidenceCoefficients) {
  const double d = 0.8, alpha0 = 3.1;
  const CoordinateMap m(d, {0.25}, 0.5);
  const int N = 40;
  const MatC D = m.fromPhysical(alpha0, N);
  for (double X : {0.0, 0.13, 0.4, 0.77}) {
    cplx u = 0;
    for (int k = -N; k <= N; ++k) u += D(k + N, N) * std::exp(fg::I * ((alpha0 + 2 * fg::PI * k / d) * X));
    EXPECT_NEAR(std::abs(u - std::exp(fg::I * (alpha0 * m.f(X)))), 0.0, 1e-12) << X;
  }
}

TEST(CoordinateMap, MapIntervals) {
  const CoordinateMap m(1.0, {0.3, 0.8}, 0.7);
  // breakpoints of the map: unchanged; values follow the physical intervals
  std::vector<double> br{0.0, 0.3, 0.8, 1.0};
  std::vector<int> v{1, 2, 3};
  fg::mapIntervals(&m, br, v);
  ASSERT_EQ(br.size(), v.size() + 1);
  EXPECT_EQ(br.front(), 0.0);
  EXPECT_EQ(br.back(), 1.0);
  // general breakpoints: the value at X is the physical value at f(X)
  std::vector<double> br2{0.0, 0.1, 0.55, 0.9, 1.0};
  std::vector<int> v2{1, 2, 3, 4};
  auto phys = [&](double x) {
    x -= std::floor(x);
    return x < 0.1 ? 1 : x < 0.55 ? 2 : x < 0.9 ? 3 : 4;
  };
  const auto b0 = br2;
  fg::mapIntervals(&m, br2, v2);
  for (std::size_t k = 0; k + 1 < br2.size(); ++k) {
    ASSERT_LT(br2[k], br2[k + 1]);
    const double X = 0.5 * (br2[k] + br2[k + 1]);
    EXPECT_EQ(v2[k], phys(m.f(X))) << X;
  }
  (void)b0;
}

// Plan "Ersetzungsregel": with A = diag(alpha), G = diag(gamma) the matrix path equals the diagonal
// path (1e-15).
TEST(PropagationMatrixTransformed, DiagonalAGEqualsDiagonalPath) {
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> u(-1, 1);
  const int n = 9;
  fg::QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Q(i, j) = MatC::Zero(n, n);
      for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) Q(i, j)(a, b) = cplx(u(rng), u(rng)) * (i == j ? 1.0 : 0.3);
      if (i == j) Q(i, j) += 3.0 * MatC::Identity(n, n);
    }
  VecC al(n), ga(n);
  for (int r = 0; r < n; ++r) {
    al(r) = 0.3 + 0.7 * (r - 4);
    ga(r) = 0.2 - 0.4 * (r % 3);
  }
  const MatC M0 = fg::buildPropagationMatrix(Q, al, ga);
  const MatC M1 = fg::buildPropagationMatrix(Q, MatC(al.asDiagonal()), MatC(ga.asDiagonal()));
  EXPECT_LT((M0 - M1).cwiseAbs().maxCoeff(), 1e-15 * M0.cwiseAbs().maxCoeff());
  VecC F(4 * n);
  for (int k = 0; k < 4 * n; ++k) F(k) = cplx(u(rng), u(rng));
  EXPECT_LT((fg::computeEy(Q, al, ga, F) - fg::computeEy(Q, MatC(al.asDiagonal()), MatC(ga.asDiagonal()), F)).norm(),
            1e-14);
  EXPECT_LT((fg::computeHy(al, ga, F) - fg::computeHy(MatC(al.asDiagonal()), MatC(ga.asDiagonal()), F)).norm(), 1e-15);
}

// Homogeneous medium in transformed coordinates: the eigenvalues of M approach +-i beta_n of the
// physical orders (central orders, 1D).
TEST(PropagationMatrixTransformed, HomogeneousEigenvalues) {
  const double d = 0.7, k0 = 2 * fg::PI / 0.6328, eta = 0.5;
  const CoordinateMap m(d, {0.2}, eta);
  const int N = 30, nf = 2 * N + 1;
  const double alpha0 = std::sin(20 * fg::PI / 180), K = 2 * fg::PI / (k0 * d);
  const auto c = m.inverseDerivativeCoefficients(2 * N);
  VecC al(nf), ga = VecC::Constant(nf, 0.1);
  for (int r = 0; r < nf; ++r) al(r) = alpha0 + (r - N) * K;
  MatC A(nf, nf), G = MatC(ga.asDiagonal());
  for (int r = 0; r < nf; ++r)
    for (int s = 0; s < nf; ++s) A(r, s) = c[static_cast<std::size_t>(r - s + 2 * N)] * al(s);
  const cplx eps = 2.25;
  const MatC M = fg::buildPropagationMatrix(fg::homogeneousIsotropicQ(eps, nf), A, G);
  Eigen::ComplexEigenSolver<MatC> es(M, false);
  for (int n = -2; n <= 2; ++n) {
    const cplx beta = fg::betaBranch(eps, al(n + N), ga(0));
    double best = 1e300;
    for (Eigen::Index k = 0; k < es.eigenvalues().size(); ++k)
      best = std::min(best, std::abs(es.eigenvalues()(k) - fg::I * beta));
    EXPECT_LT(best, 1e-8) << n;
  }
}

// A rectangle stays a rectangle under a separable map: mapped contour integral = exact rectangle
// coefficients of [f^{-1}(x1), f^{-1}(x2)] x [g^{-1}(z1), g^{-1}(z2)].
TEST(ContourTransformed, RectangleCoefficients) {
  const double dx = 1.0, dz = 0.8;
  const CoordinateMap fx(dx, {0.1, 0.45}, 0.7), fz(dz, {0.3}, 0.4);
  fg::geometry::Region r;
  r.outer = fg::geometry::Contour2D::polygon({{0.2, 0.1}, {0.6, 0.1}, {0.6, 0.5}, {0.2, 0.5}});
  r.normalize(1e-9);
  const auto g = fg::fourier::regionCoefficients(r, dx, dz, 6, 5, &fx, &fz);
  for (int p = -6; p <= 6; ++p)
    for (int q = -5; q <= 5; ++q) {
      const cplx ref = fg::fourier::rectangleCoefficient(fx.inverse(0.2), fx.inverse(0.6), fz.inverse(0.1),
                                                         fz.inverse(0.5), dx, dz, p, q);
      EXPECT_NEAR(std::abs(g(p, q) - ref), 0.0, 1e-13) << p << " " << q;
    }
}

