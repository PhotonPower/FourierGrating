#include <gtest/gtest.h>

#include <cmath>

#include "fg/fourier/coefficients.hpp"

using fg::cplx;

namespace {
// Midpoint-rule reference for the coefficient of the indicator of [x1, x2].
cplx numericInterval(double x1, double x2, double d, int p) {
  const int M = 200000;
  cplx s = 0;
  const double K = 2 * fg::PI / d;
  for (int k = 0; k < M; ++k) {
    const double x = (k + 0.5) * d / M;
    if (x >= x1 && x <= x2) s += std::exp(-fg::I * (p * K * x));
  }
  return s / double(M);
}
}  // namespace

TEST(Coefficients, IntervalMatchesQuadrature) {
  for (int p : {-7, -1, 0, 2, 5})
    EXPECT_NEAR(std::abs(fg::fourier::intervalCoefficient(0.13, 0.58, 1.0, p) -
                         numericInterval(0.13, 0.58, 1.0, p)),
                0.0, 2e-5);
}

TEST(Coefficients, ComplementaryRegionsSumToDelta) {
  // eps = 1 on [0, 0.3], 1 on [0.3, 1] -> coefficients delta_p0
  const auto c = fg::fourier::piecewiseConstant1D({0.0, 0.3, 1.0}, {1.0, 1.0}, 1.0, 6);
  for (int p = -6; p <= 6; ++p)
    EXPECT_NEAR(std::abs(c[static_cast<std::size_t>(p + 6)] - (p == 0 ? 1.0 : 0.0)), 0.0, 1e-14);
}

TEST(Coefficients, PolygonEqualsRectangle) {
  const double dx = 0.9, dz = 0.7;
  const std::vector<std::array<double, 2>> rect{{0.1, 0.2}, {0.5, 0.2}, {0.5, 0.45}, {0.1, 0.45}};
  for (int p = -4; p <= 4; ++p)
    for (int q = -3; q <= 3; ++q) {
      const cplx a = fg::fourier::polygonCoefficient(rect, dx, dz, p, q);
      const cplx b = fg::fourier::rectangleCoefficient(0.1, 0.5, 0.2, 0.45, dx, dz, p, q);
      EXPECT_NEAR(std::abs(a - b), 0.0, 1e-13) << "p=" << p << " q=" << q;
    }
}

TEST(Coefficients, PolygonTriangleArea) {
  const std::vector<std::array<double, 2>> tri{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  EXPECT_NEAR(std::real(fg::fourier::polygonCoefficient(tri, 2.0, 2.0, 0, 0)), 0.5 / 4.0, 1e-15);
}
