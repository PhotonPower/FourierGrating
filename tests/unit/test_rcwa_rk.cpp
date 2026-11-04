#include <gtest/gtest.h>

#include <random>
#include <unsupported/Eigen/MatrixFunctions>

#include "fg/fourier/coefficients.hpp"
#include "fg/physics/qmatrix.hpp"
#include "fg/solvers/rcwa.hpp"
#include "fg/solvers/rk.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;

namespace {
// Isotropic lamellar layer, conical incidence: M is block anti-diagonal.
struct Fixture {
  fg::IndexMap map{6, 0};
  VecC alpha, gamma;
  MatC M;
  fg::ModeBasis ref;
  Fixture() {
    const int nf = map.size(), P = 2 * map.Nx();
    alpha.resize(nf);
    gamma.resize(nf);
    for (int r = 0; r < nf; ++r) {
      alpha(r) = 0.3 + map.n(r) * 0.8;
      gamma(r) = 0.25;
    }
    const auto ce = fg::fourier::piecewiseConstant1D({0.0, 0.5, 1.0}, {cplx(4.0, 0.1), 1.0}, 1.0, P);
    const auto ci = fg::fourier::piecewiseConstant1D({0.0, 0.5, 1.0}, {1.0 / cplx(4.0, 0.1), 1.0}, 1.0, P);
    fg::SliceCoefficients sc;
    sc.eps = fg::fourier::CoefficientGrid(P, 0);
    sc.invEps = fg::fourier::CoefficientGrid(P, 0);
    for (int p = -P; p <= P; ++p) {
      sc.eps(p, 0) = ce[static_cast<std::size_t>(p + P)];
      sc.invEps(p, 0) = ci[static_cast<std::size_t>(p + P)];
    }
    M = fg::buildPropagationMatrix(fg::qLiLamellar(sc, map), alpha, gamma);
    ref = fg::isotropicModes(1.0, alpha, gamma);
  }
};
}  // namespace

TEST(RCWA, ReducedModesAreEigenvectors) {
  Fixture f;
  ASSERT_TRUE(fg::isBlockAntiDiagonal(f.M));
  const auto m = fg::layerModesReduced(f.M);
  EXPECT_LT((f.M * m.Wu - m.Wu * m.lu.asDiagonal()).cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT((f.M * m.Wd - m.Wd * m.ld.asDiagonal()).cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LE(m.lu.real().maxCoeff(), 1e-12);
  EXPECT_GE(m.ld.real().minCoeff(), -1e-12);
}

TEST(RCWA, StableSMatrixEqualsTransferSMatrixForThinLayer) {
  Fixture f;
  const double h = 0.2;
  const auto S1 = fg::stableLayerSMatrix(fg::layerModesReduced(f.M), h, f.ref);
  const auto S2 = fg::fromTransfer(f.ref.PsiInv * (h * f.M).exp() * f.ref.Psi);
  EXPECT_LT((S1.full() - S2.full()).cwiseAbs().maxCoeff(), 1e-9);
  const auto S3 = fg::stableLayerSMatrix(fg::layerModesGeneral(f.M), h, f.ref);
  EXPECT_LT((S1.full() - S3.full()).cwiseAbs().maxCoeff(), 1e-9);
}

TEST(RCWA, ThickLayerStaysFinite) {
  Fixture f;
  const auto S = fg::stableLayerSMatrix(fg::layerModesReduced(f.M), 200.0, f.ref);
  EXPECT_TRUE(S.full().allFinite());
}

TEST(RK, ConstantTransferMatchesExponential) {
  Fixture f;
  fg::RKStats st;
  const double h = 0.3;
  const MatC T = fg::rkTransferConstant(f.M, h, fg::ButcherTableau::gauss2(), 1e-12, 5.0, st);
  const MatC E = (h * f.M).exp();
  EXPECT_LT((T - E).norm() / E.norm(), 1e-10);
  fg::RKStats st2;
  const MatC T2 = fg::rkTransfer([&](double) { return f.M; }, 0.0, h, fg::ButcherTableau::gauss3(),
                                 1e-13, 5.0, st2);
  EXPECT_LT((T2 - E).norm() / E.norm(), 1e-10);
}

TEST(RK, VariableCoefficientsAgainstReference) {
  // dF/dy = (A + y B) F: compare adaptive Gauss-Legendre with a very fine reference
  std::mt19937 rng(5);
  std::normal_distribution<double> nd;
  MatC A(4, 4), B(4, 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      A(i, j) = cplx(nd(rng), nd(rng));
      B(i, j) = cplx(nd(rng), nd(rng));
    }
  auto M = [&](double y) { return MatC(A + y * B); };
  fg::RKStats st;
  const MatC T = fg::rkTransfer(M, 0.0, 1.0, fg::ButcherTableau::gauss2(), 1e-12, 5.0, st);
  MatC R = MatC::Identity(4, 4);
  const int n = 20000;
  for (int k = 0; k < n; ++k) R = (MatC(M((k + 0.5) / n)) / double(n)).exp() * R;  // midpoint Magnus
  EXPECT_LT((T - R).norm() / R.norm(), 1e-6);
}
