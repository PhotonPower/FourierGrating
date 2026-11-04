#include <gtest/gtest.h>

#include <cmath>

#include "fg/fourier/coefficients.hpp"
#include "fg/fourier/toeplitz.hpp"

using fg::cplx;
using fg::MatC;

TEST(Toeplitz, BookExampleIV1) {
  // Book IV.1: f = 2 on [0, d/2), 1/2 on [d/2, d), g = 1/f, h = f g = 1.
  // N = 0: Laurent gives f0 g0 = 1.5625, inverse rule gives exactly 1.
  const auto f = fg::fourier::piecewiseConstant1D({0.0, 0.5, 1.0}, {2.0, 0.5}, 1.0, 0);
  const auto g = fg::fourier::piecewiseConstant1D({0.0, 0.5, 1.0}, {0.5, 2.0}, 1.0, 0);
  const MatC F = fg::fourier::toeplitz1D(f, 0), G = fg::fourier::toeplitz1D(g, 0);
  EXPECT_NEAR(std::real((F * G)(0, 0)), 1.5625, 1e-14);
  // inverse rule: [[1/f]]^-1 [[g]] with 1/f = g
  EXPECT_NEAR(std::real((G.inverse() * G)(0, 0)), 1.0, 1e-14);
}

TEST(Toeplitz, LaurentInteriorForSmoothFunctions) {
  // f = exp(cos x), g = 1 + 0.5 sin x: interior of [[f]][[g]] equals [[fg]] (Plan, Toeplitz remark)
  const int N = 20, P = 4 * N, Ms = 4096;
  auto coeffs = [&](auto fun, int pmax) {
    std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1), 0.0);
    for (int p = -pmax; p <= pmax; ++p) {
      cplx s = 0;
      for (int k = 0; k < Ms; ++k) {
        const double x = 2 * fg::PI * k / Ms;
        s += fun(x) * std::exp(-fg::I * (p * x));
      }
      c[static_cast<std::size_t>(p + pmax)] = s / double(Ms);
    }
    return c;
  };
  auto f = [](double x) { return std::exp(std::cos(x)); };
  auto g = [](double x) { return 1.0 + 0.5 * std::sin(x); };
  auto fg_ = [&](double x) { return f(x) * g(x); };
  auto cut = [&](const std::vector<cplx>& c) {  // keep p = -2N..2N from p = -P..P
    return std::vector<cplx>(c.begin() + (P - 2 * N), c.begin() + (P + 2 * N) + 1);
  };
  const MatC F = fg::fourier::toeplitz1D(cut(coeffs(f, P)), N);
  const MatC G = fg::fourier::toeplitz1D(cut(coeffs(g, P)), N);
  const MatC H = fg::fourier::toeplitz1D(cut(coeffs(fg_, P)), N);
  const MatC D = F * G - H;
  EXPECT_LT(D.block(N / 2, N / 2, N + 1, N + 1).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(Toeplitz, TwoDimensionalBlockStructure) {
  const fg::IndexMap map(2, 1);
  const MatC T = fg::fourier::toeplitz(map, [](int p, int q) { return cplx(p, q); });
  for (int r = 0; r < map.size(); ++r)
    for (int s = 0; s < map.size(); ++s) {
      const auto [n1, m1] = map.orders(r);
      const auto [n2, m2] = map.orders(s);
      EXPECT_EQ(T(r, s), cplx(n1 - n2, m1 - m2));
    }
}
