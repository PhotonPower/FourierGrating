#include <gtest/gtest.h>

#include <random>

#include "fg/physics/anisotropic.hpp"
#include "fg/physics/modes.hpp"
#include "fg/physics/propagation_matrix.hpp"

using fg::cplx;
using fg::Mat3C;
using fg::MatC;
using fg::Vec3D;
using fg::VecC;

namespace {
Mat3C randomTensor(std::mt19937& rng) {
  std::normal_distribution<double> nd;
  Mat3C e;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) e(i, j) = cplx(nd(rng), nd(rng)) + (i == j ? 3.0 : 0.0);
  return e;
}
}  // namespace

TEST(Anisotropic, CMatrixIsInverseOfAMatrix) {
  std::mt19937 rng(1);
  for (int k = 0; k < 20; ++k) {
    const Mat3C e = randomTensor(rng);
    const double t = 0.3 * k;
    // profile normal (N_z = 0) and y-invariant normal (N_y = 0, T1 = y, Plan eq:T_yinv)
    for (int variant = 0; variant < 2; ++variant) {
      Vec3D N, T1, T2;
      if (variant == 0) {
        N = Vec3D(-std::sin(t), std::cos(t), 0);
        fg::profileTangents(N, T1, T2);
      } else {
        N = Vec3D(std::cos(t), 0, std::sin(t));
        T1 = Vec3D(0, 1, 0);
        T2 = T1.cross(N);
      }
      ASSERT_LT((T2.cross(T1) - N).norm(), 1e-14);
      Mat3C A;
      A.row(0) = T1.cast<cplx>().transpose();
      A.row(1) = N.cast<cplx>().transpose() * e;
      A.row(2) = T2.cast<cplx>().transpose();
      EXPECT_LT((A * fg::cMatrix(e, N, T1, T2) - Mat3C::Identity()).cwiseAbs().maxCoeff(), 1e-13);
    }
  }
}

TEST(Anisotropic, CMatrixIsotropicLimit) {
  const cplx e(-11.5, 1.4);
  const Vec3D N(-0.6, 0.8, 0);
  Vec3D T1, T2;
  fg::profileTangents(N, T1, T2);
  Mat3C ref;
  ref.col(0) = T1.cast<cplx>();
  ref.col(1) = N.cast<cplx>() / e;
  ref.col(2) = T2.cast<cplx>();
  EXPECT_LT((fg::cMatrix(e * Mat3C::Identity(), N, T1, T2) - ref).cwiseAbs().maxCoeff(), 1e-15);
}

TEST(Anisotropic, ProjectorReconstruction) {
  // E = H1 e_t + d_n / xi, D = G1 e_t + G2 d_n  with e_t = P_T E, d_n = P_N D (Plan eq:projector)
  std::mt19937 rng(2);
  std::normal_distribution<double> nd;
  for (int k = 0; k < 10; ++k) {
    const Mat3C e = randomTensor(rng);
    Vec3D N(nd(rng), nd(rng), nd(rng));
    N.normalize();
    const fg::Vec3C E(cplx(nd(rng), nd(rng)), cplx(nd(rng), nd(rng)), cplx(nd(rng), nd(rng)));
    const fg::Vec3C D = e * E;
    const auto p = fg::projectorMatrices(e, N);
    const Mat3C PN = (N * N.transpose()).cast<cplx>();
    const fg::Vec3C et = p.PT * E, dn = PN * D;
    EXPECT_LT((p.H1 * et + p.invXi * dn - E).norm(), 1e-12);
    EXPECT_LT((p.G1 * et + p.G2 * dn - D).norm(), 1e-12);
  }
}

TEST(Anisotropic, RotatedTensorKeepsEigenvalues) {
  const Mat3C e = fg::rotatedTensor(2.0, 3.0, 5.0, 0.3, 1.1, -0.7);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(e.real());
  EXPECT_NEAR(es.eigenvalues()(0), 2.0, 1e-13);
  EXPECT_NEAR(es.eigenvalues()(2), 5.0, 1e-13);
  EXPECT_LT(e.imag().cwiseAbs().maxCoeff(), 1e-15);
}

TEST(AnisotropicModes, IsotropicTensorReproducesIsotropicModes) {
  const int n = 4;
  VecC al(n), ga(n);
  for (int r = 0; r < n; ++r) {
    al(r) = 0.2 + 0.9 * (r - 1);
    ga(r) = 0.3;
  }
  const cplx eps(2.25, 0.1);
  const auto a = fg::anisotropicModes(eps * Mat3C::Identity(), al, ga);
  const auto b = fg::isotropicModes(eps, al, ga);
  const MatC M = fg::buildPropagationMatrix(fg::homogeneousIsotropicQ(eps, n), al, ga);
  EXPECT_LT((M * a.Psi - a.Psi * a.lambda.asDiagonal()).cwiseAbs().maxCoeff(), 1e-12);
  EXPECT_LT((a.Psi * a.PsiInv - MatC::Identity(4 * n, 4 * n)).cwiseAbs().maxCoeff(), 1e-12);
  // same eigenvalues in the same (down | up) halves
  for (int r = 0; r < n; ++r)
    for (int c : {0, 1}) {
      EXPECT_NEAR(std::abs(a.lambda(c * n + r) - b.lambda(c * n + r)), 0.0, 1e-12);
      EXPECT_NEAR(std::abs(a.lambda((2 + c) * n + r) - b.lambda((2 + c) * n + r)), 0.0, 1e-12);
    }
}

TEST(AnisotropicModes, UniaxialAlongYAnalytic) {
  // eps = diag(eo, ee, eo): TE-like k_y^2 = eo - kappa^2, TM-like k_y^2 = eo (1 - alpha^2 / ee)
  const double eo = 2.4, ee = 3.1;
  Mat3C e = Mat3C::Zero();
  e(0, 0) = eo;
  e(1, 1) = ee;
  e(2, 2) = eo;
  VecC al(1), ga(1);
  al(0) = 0.7;
  ga(0) = 0.0;
  const auto m = fg::anisotropicModes(e, al, ga);
  std::vector<double> ky;
  for (int c = 2; c < 4; ++c) ky.push_back(std::abs(m.lambda(c)));  // up modes: lambda = i k_y
  std::sort(ky.begin(), ky.end());
  std::vector<double> ref{std::sqrt(eo - 0.49), std::sqrt(eo * (1 - 0.49 / ee))};
  std::sort(ref.begin(), ref.end());
  EXPECT_NEAR(ky[0], ref[0], 1e-12);
  EXPECT_NEAR(ky[1], ref[1], 1e-12);
  for (int c = 2; c < 4; ++c) {  // up-going modes carry flux in +y
    VecC F = m.Psi.col(c);
    EXPECT_GT(fg::poyntingY(F), 0.0);
  }
}
