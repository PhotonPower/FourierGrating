#include <gtest/gtest.h>

#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>

#include "fg/solvers/butcher.hpp"

using fg::MatC;

class ButcherTest : public ::testing::TestWithParam<std::string> {};

TEST_P(ButcherTest, OrderConditionsQuadrature) {
  const auto t = fg::ButcherTableau::byName(GetParam());
  EXPECT_LT((t.A.rowwise().sum() - t.c).cwiseAbs().maxCoeff(), 1e-14);
  for (int q = 1; q <= t.order; ++q)
    EXPECT_NEAR(t.b.dot(t.c.array().pow(q - 1).matrix()), 1.0 / q, 1e-14) << "q=" << q;
}

TEST_P(ButcherTest, LocalErrorOrderForConstantM) {
  const auto t = fg::ButcherTableau::byName(GetParam());
  MatC M = MatC::Random(6, 6);
  double err[2];
  const double hs[2] = {0.1, 0.05};
  for (int k = 0; k < 2; ++k) {
    const std::vector<MatC> nodes(static_cast<std::size_t>(t.stages()), M);
    const MatC P = fg::rkStepPropagator(nodes, hs[k], t);
    const MatC E = (hs[k] * M).exp();
    err[k] = (P - E).norm();
  }
  const double observed = std::log2(err[0] / err[1]);
  EXPECT_NEAR(observed, t.order + 1, 0.3);
}

INSTANTIATE_TEST_SUITE_P(All, ButcherTest, ::testing::Values("gauss2", "gauss3", "radau3"));
