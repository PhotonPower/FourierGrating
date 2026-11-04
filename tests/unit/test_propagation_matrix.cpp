#include <gtest/gtest.h>

#include <random>

#include "fg/physics/modes.hpp"
#include "fg/physics/propagation_matrix.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;

TEST(PropagationMatrix, IsotropicModesAreEigenvectors) {
  const int n = 5;
  const cplx eps(2.3, 0.4);
  VecC alpha(n), gamma(n);
  for (int r = 0; r < n; ++r) {
    alpha(r) = 0.31 + 0.9 * (r - 2);
    gamma(r) = -0.52 + 0.2 * r;
  }
  const MatC M = fg::buildPropagationMatrix(fg::homogeneousIsotropicQ(eps, n), alpha, gamma);
  const auto mb = fg::isotropicModes(eps, alpha, gamma);
  EXPECT_LT((M * mb.Psi - mb.Psi * mb.lambda.asDiagonal()).cwiseAbs().maxCoeff(), 1e-12);
  EXPECT_LT((mb.Psi * mb.PsiInv - MatC::Identity(4 * n, 4 * n)).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(PropagationMatrix, MatchesDirectElimination) {
  // Random Q blocks: M F must equal the derivatives obtained directly from Maxwell (Plan 3.3/3.4).
  std::mt19937 rng(11);
  std::normal_distribution<double> nd;
  const int n = 4;
  auto rnd = [&](double diag) {
    MatC A(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) A(i, j) = cplx(nd(rng), nd(rng)) + (i == j ? diag : 0.0);
    return A;
  };
  fg::QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Q(i, j) = rnd(i == j ? 4.0 : 0.0);
  VecC alpha(n), gamma(n), F(4 * n);
  for (int r = 0; r < n; ++r) {
    alpha(r) = nd(rng);
    gamma(r) = nd(rng);
  }
  for (int r = 0; r < 4 * n; ++r) F(r) = cplx(nd(rng), nd(rng));
  const MatC M = fg::buildPropagationMatrix(Q, alpha, gamma);
  const VecC Ex = F.segment(0, n), Ez = F.segment(n, n), Hx = F.segment(2 * n, n),
             Hz = F.segment(3 * n, n);
  const VecC Ey = fg::computeEy(Q, alpha, gamma, F);
  const VecC Hy = fg::computeHy(alpha, gamma, F);
  const VecC Dx = Q(0, 0) * Ex + Q(0, 1) * Ey + Q(0, 2) * Ez;
  const VecC Dz = Q(2, 0) * Ex + Q(2, 1) * Ey + Q(2, 2) * Ez;
  const cplx i(0, 1);
  VecC dF(4 * n);
  dF << i * alpha.cwiseProduct(Ey) - i * Hz,  // dEx/dy = i alpha Ey - i Bz
      i * gamma.cwiseProduct(Ey) + i * Hx,     // dEz/dy = i gamma Ey + i Bx
      i * alpha.cwiseProduct(Hy) + i * Dz,     // dHx/dy = i alpha Hy + i Dz
      i * gamma.cwiseProduct(Hy) - i * Dx;     // dHz/dy = i gamma Hy - i Dx
  EXPECT_LT((M * F - dF).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(Modes, BetaBranch) {
  EXPECT_GE(fg::betaBranch(1.0, 1.5, 0.0).imag(), 0.0);           // evanescent: Im > 0
  EXPECT_GT(fg::betaBranch(2.25, 0.3, 0.2).real(), 0.0);          // propagating: Re > 0
  EXPECT_GT(std::abs(fg::betaBranch(1.0, 1.0, 0.0)), 0.0);        // Rayleigh anomaly regularized
}
