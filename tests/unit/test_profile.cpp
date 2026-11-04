#include <gtest/gtest.h>

#include <cmath>

#include "fg/fourier/sampling.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/geometry/profile.hpp"

using fg::Profile1DLayer;

TEST(Profile, SinusCrossingsAndIntervals) {
  Profile1DLayer p;
  p.depth_um = 0.2;
  p.x0_um = 0.1;
  const double d = 0.5;
  for (double y : {0.01, 0.1, 0.19}) {
    const auto xs = fg::profileCrossings(p, y, d);
    ASSERT_EQ(xs.size(), 2u);
    for (double x : xs) EXPECT_NEAR(fg::profileValue(p, x, d), y, 1e-13);
    std::vector<double> br;
    std::vector<bool> below;
    fg::profileIntervals(p, y, d, br, below);
    double fill = 0;
    for (std::size_t k = 0; k < below.size(); ++k)
      if (below[k]) fill += br[k + 1] - br[k];
    // fraction below the profile: g > y  <=>  |K(x - x0)| < acos(2y/a - 1)
    EXPECT_NEAR(fill / d, std::acos(2 * y / p.depth_um - 1) / fg::PI, 1e-13);
  }
}

TEST(Profile, TriangleCrossings) {
  Profile1DLayer p;
  p.shape = Profile1DLayer::Shape::PiecewiseLinear;
  p.depth_um = 0.3;
  p.xs_um = {0.0, 0.7};
  p.ys_um = {0.0, 0.3};
  const auto xs = fg::profileCrossings(p, 0.15, 1.0);
  ASSERT_EQ(xs.size(), 2u);
  EXPECT_NEAR(xs[0], 0.35, 1e-14);
  EXPECT_NEAR(xs[1], 0.85, 1e-14);
  EXPECT_NEAR(fg::profileSlope(p, 0.2, 1.0), 0.3 / 0.7, 1e-14);
  EXPECT_NEAR(fg::profileSlope(p, 0.9, 1.0), -0.3 / 0.3, 1e-14);
}

TEST(Profile, FunctionShapeMatchesSinus) {
  Profile1DLayer s, f;
  s.depth_um = f.depth_um = 0.2;
  f.shape = Profile1DLayer::Shape::Function;
  f.g = [](double x) { return 0.1 * (1 + std::cos(2 * fg::PI * x / 0.5)); };
  const auto a = fg::profileCrossings(s, 0.07, 0.5), b = fg::profileCrossings(f, 0.07, 0.5);
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t k = 0; k < a.size(); ++k) EXPECT_NEAR(a[k], b[k], 1e-12);
}

TEST(NormalField, PiecewiseLinearExactVsSampled) {
  // exact interval coefficients of N_i N_j vs. dense sampling (slowly converging, O(1/M))
  Profile1DLayer p;
  p.shape = Profile1DLayer::Shape::PiecewiseLinear;
  p.depth_um = 0.3;
  p.xs_um = {0.0, 0.7};
  p.ys_um = {0.0, 0.3};
  const fg::IndexMap map(5, 0);
  const auto nf = fg::normalField(p, 1.0, map);
  const auto ref = fg::fourier::sampledCoefficients1D(
      [&](double x) {
        const double s = fg::profileSlope(p, x, 1.0);
        return fg::cplx(s * s / (1 + s * s), 0);
      },
      1.0, 10, 200000);
  for (int q = -10; q <= 10; ++q) EXPECT_NEAR(std::abs(nf.NN[0][0](q, 0) - ref[static_cast<std::size_t>(q + 10)]), 0.0, 1e-4);
  // N_y^2 + N_x^2 = 1 everywhere -> coefficients delta_q0
  for (int q = -10; q <= 10; ++q)
    EXPECT_NEAR(std::abs(nf.NN[0][0](q, 0) + nf.NN[1][1](q, 0) - (q == 0 ? 1.0 : 0.0)), 0.0, 1e-14);
}
