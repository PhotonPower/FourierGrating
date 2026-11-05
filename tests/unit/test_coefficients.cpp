#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <random>

#include "fg/fourier/coefficients.hpp"
#include "fg/fourier/sampling.hpp"

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

// FFT of any length (radix 2 / Bluestein, Plan "Allgemeine Funktionen") against the direct DFT sum.
TEST(Sampling, FFTEqualsDirectDFT) {
  std::mt19937 rng(7);
  std::normal_distribution<double> nd;
  for (int n : {1, 2, 3, 7, 16, 97, 100, 336, 1024, 1025}) {
    std::vector<cplx> a(static_cast<std::size_t>(n));
    for (auto& v : a) v = cplx(nd(rng), nd(rng));
    const auto X = fg::fourier::fft(a);
    double err = 0, ref = 0;
    for (int p = 0; p < n; ++p) {
      cplx s = 0;
      for (int k = 0; k < n; ++k)
        s += a[static_cast<std::size_t>(k)] *
             std::exp(-fg::I *
                      (2 * fg::PI * static_cast<double>((static_cast<long>(p) * k) % n) / n));
      err = std::max(err, std::abs(X[static_cast<std::size_t>(p)] - s));
      ref = std::max(ref, std::abs(s));
    }
    EXPECT_LT(err, 1e-13 * ref * std::log2(n + 1.0)) << "n = " << n;
  }
  // large prime length (Bluestein, k^2 beyond 32 bit): selected outputs against the direct sum
  const int n = 100003;
  std::vector<cplx> a(static_cast<std::size_t>(n));
  for (auto& v : a) v = cplx(nd(rng), nd(rng));
  const auto X = fg::fourier::fft(a);
  for (int p : {0, 1, 777, n - 1}) {
    cplx s = 0;
    for (int k = 0; k < n; ++k)
      s += a[static_cast<std::size_t>(k)] *
           std::exp(-fg::I *
                    (2 * fg::PI * static_cast<double>((static_cast<std::int64_t>(p) * k) % n) / n));
    EXPECT_LT(std::abs(X[static_cast<std::size_t>(p)] - s),
              1e-9 * std::sqrt(static_cast<double>(n)))
        << "p = " << p;
  }
}

// Sampled coefficients of trigonometric polynomials are exact (no aliasing for M > 2 pmax),
// also for sample counts that are not powers of two; both paths (fft_sampling off / on).
TEST(Sampling, SampledCoefficientsExactForTrigPolynomials) {
  for (const bool on : {false, true}) {
    fg::fourier::setFFTSampling(on);
    const double d = 0.7, K = 2 * fg::PI / d;
    auto f1 = [&](double x) {
      return cplx(1.5, 0) + cplx(0.3, -0.2) * std::exp(fg::I * (3 * K * x)) +
             0.4 * std::cos(5 * K * x);
    };
    for (int M : {64, 97, 336}) {
      const auto c = fg::fourier::sampledCoefficients1D(f1, d, 6, M);
      EXPECT_NEAR(std::abs(c[6] - 1.5), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(c[6 + 3] - cplx(0.3, -0.2)), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(c[6 + 5] - 0.2), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(c[6 - 5] - 0.2), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(c[6 - 3]), 0.0, 1e-14);
    }
    const double dz = 0.4, Kz = 2 * fg::PI / dz;
    auto f2 = [&](double x, double z) {
      return 2.0 + std::exp(fg::I * (2 * K * x - 3 * Kz * z)) +
             cplx(0, 0.5) * std::exp(fg::I * (-K * x + Kz * z));
    };
    for (const auto& [Mx, Mz] : {std::pair{32, 64}, std::pair{45, 34}}) {
      const auto g = fg::fourier::sampledCoefficients2D(f2, d, dz, 4, 5, Mx, Mz);
      EXPECT_NEAR(std::abs(g(0, 0) - 2.0), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(g(2, -3) - 1.0), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(g(-1, 1) - cplx(0, 0.5)), 0.0, 1e-14);
      EXPECT_NEAR(std::abs(g(1, -1)), 0.0, 1e-14);
    }
  }
  fg::fourier::setFFTSampling(false);
}

// fft_sampling on = off (v0.10 DFT) up to rounding, for powers of two and other grid sizes.
TEST(Sampling, FFTSwitchMatchesDFT) {
  std::mt19937 rng(9);
  std::normal_distribution<double> nd;
  for (const auto& [Mx, Mz] : {std::pair{64, 32}, std::pair{84, 51}}) {
    std::vector<cplx> v(static_cast<std::size_t>(Mx * Mz));
    for (auto& x : v) x = cplx(nd(rng), nd(rng));
    fg::fourier::setFFTSampling(false);
    const auto a = fg::fourier::gridCoefficients2D(v, Mx, Mz, 10, 7);
    fg::fourier::setFFTSampling(true);
    const auto b = fg::fourier::gridCoefficients2D(v, Mx, Mz, 10, 7);
    for (int p = -10; p <= 10; ++p)
      for (int q = -7; q <= 7; ++q)
        EXPECT_NEAR(std::abs(a(p, q) - b(p, q)), 0.0, 1e-14) << Mx << " " << Mz;
  }
  auto f = [](double x) {
    return cplx(std::exp(std::cos(2 * fg::PI * x)), std::sin(4 * fg::PI * x));
  };
  for (int M : {1024, 1000}) {
    fg::fourier::setFFTSampling(false);
    const auto a = fg::fourier::sampledCoefficients1D(f, 1.0, 20, M);
    fg::fourier::setFFTSampling(true);
    const auto b = fg::fourier::sampledCoefficients1D(f, 1.0, 20, M);
    // v0.10 DFT uses the twiddle recurrence w^k (error ~ M eps)
    for (std::size_t k = 0; k < a.size(); ++k) EXPECT_NEAR(std::abs(a[k] - b[k]), 0.0, 1e-12) << M;
  }
  fg::fourier::setFFTSampling(false);
}
