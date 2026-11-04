#include <gtest/gtest.h>

#include <cmath>

#include "fg/fourier/sampling.hpp"
#include "fg/fourier/toeplitz.hpp"
#include "fg/physics/qmatrix.hpp"

using fg::cplx;
using fg::MatC;

namespace {
fg::SliceCoefficients lamellar(const fg::IndexMap& map) {
  const int P = 2 * map.Nx();
  const auto ce = fg::fourier::piecewiseConstant1D({0.0, 0.4, 1.0}, {cplx(-20, 3), 1.0}, 1.0, P);
  const auto ci = fg::fourier::piecewiseConstant1D({0.0, 0.4, 1.0}, {1.0 / cplx(-20, 3), 1.0}, 1.0, P);
  fg::SliceCoefficients sc;
  sc.eps = fg::fourier::CoefficientGrid(P, 0);
  sc.invEps = fg::fourier::CoefficientGrid(P, 0);
  for (int p = -P; p <= P; ++p) {
    sc.eps(p, 0) = ce[static_cast<std::size_t>(p + P)];
    sc.invEps(p, 0) = ci[static_cast<std::size_t>(p + P)];
  }
  return sc;
}

fg::NormalFieldCoefficients constantNormal(int P, double nx, double ny, double nz) {
  fg::NormalFieldCoefficients nf;
  const double N[3] = {nx, ny, nz};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      auto& g = nf.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      g = fg::fourier::CoefficientGrid(P, 0);
      g(0, 0) = N[i] * N[j];
    }
  return nf;
}

double maxDiff(const fg::QBlocks& a, const fg::QBlocks& b) {
  double d = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) d = std::max(d, (a(i, j) - b(i, j)).cwiseAbs().maxCoeff());
  return d;
}
}  // namespace

TEST(QMatrix, FFFWithNormalAlongXEqualsLiLamellar) {
  const fg::IndexMap map(8, 0);
  auto sc = lamellar(map);
  const auto nf = constantNormal(16, 1, 0, 0);
  sc.normal = &nf;
  EXPECT_LT(maxDiff(fg::qFFFIsotropic(sc, map), fg::qLiLamellar(sc, map)), 1e-11);
}

TEST(QMatrix, FFFWithNormalAlongYIsClassicalLaurentExceptYY) {
  // N = y everywhere (book IV.2.3, a/d -> 0): Q_yy = [[1/eps]]^-1, all others Laurent
  const fg::IndexMap map(8, 0);
  auto sc = lamellar(map);
  const auto nf = constantNormal(16, 0, 1, 0);
  sc.normal = &nf;
  const auto Q = fg::qFFFIsotropic(sc, map);
  auto L = fg::qLaurent(sc, map);
  L(1, 1) = fg::fourier::toeplitz(map, sc.invEps).inverse();
  EXPECT_LT(maxDiff(Q, L), 1e-11);
}

TEST(Sampling, SmoothFunctionCoefficients) {
  // exp(cos x) has coefficients I_p(1) (modified Bessel), check I_0(1), I_1(1)
  const auto c = fg::fourier::sampledCoefficients1D(
      [](double x) { return cplx(std::exp(std::cos(x)), 0); }, 2 * fg::PI, 3, 64);
  EXPECT_NEAR(c[3].real(), 1.2660658777520082, 1e-14);
  EXPECT_NEAR(c[4].real(), 0.5651591039924851, 1e-14);
  EXPECT_NEAR(c[2].real(), 0.5651591039924851, 1e-14);
}

TEST(Sampling, TwoDimensionalSeparable) {
  const auto g = fg::fourier::sampledCoefficients2D(
      [](double x, double z) { return cplx(std::cos(2 * fg::PI * x) * std::sin(2 * fg::PI * z), 0); },
      1.0, 1.0, 2, 2, 16, 16);
  EXPECT_NEAR(std::abs(g(1, 1) - cplx(0, -0.25)), 0.0, 1e-14);
  EXPECT_NEAR(std::abs(g(0, 0)), 0.0, 1e-14);
}
