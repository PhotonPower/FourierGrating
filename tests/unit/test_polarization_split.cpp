// Decoupled 2 n_f systems of the 1D in-plane problem (Plan eq:TE, eq:TM, "Spezialfaelle"): building blocks
// against the blocks of the full 4 n_f system, solver results split = full.
#include <gtest/gtest.h>

#include <random>
#include <string>

#include "fg/io/config.hpp"
#include "fg/physics/modes.hpp"
#include "fg/physics/propagation_matrix.hpp"
#include "fg/solvers/rcwa.hpp"

using fg::cplx;
using fg::MatC;
using fg::Polarization;
using fg::VecC;

namespace {
/// Blocks of a 4 n_f matrix on the field rows/columns of one polarization.
MatC fieldBlocks(const MatC& M, Polarization p) {
  const Eigen::Index n = M.rows() / 4;
  const auto b = fg::polarizationBlocks(p);
  MatC R(2 * n, 2 * n);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      R.block(i * n, j * n, n, n) = M.block(b.field[i] * n, b.field[j] * n, n, n);
  return R;
}

/// Random Q blocks with the in-plane isotropic structure (Q_xz = Q_zx = Q_yz = Q_zy = 0) and
/// Q_xy = Q_yx != 0 as for FFF profiles.
fg::QBlocks inPlaneQ(int n, unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> nd;
  auto rnd = [&](double diag) {
    MatC A(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) A(i, j) = 0.3 * cplx(nd(rng), nd(rng)) + (i == j ? diag : 0.0);
    return A;
  };
  fg::QBlocks Q = fg::homogeneousIsotropicQ(0.0, n);
  Q(0, 0) = rnd(3.0);
  Q(1, 1) = rnd(3.0);
  Q(2, 2) = rnd(3.0);
  Q(0, 1) = rnd(0.0);
  Q(1, 0) = Q(0, 1);
  return Q;
}

VecC inPlaneAlpha(int n) {
  VecC alpha(n);
  for (int r = 0; r < n; ++r) alpha(r) = 0.27 + 1.1 * (r - n / 2);
  return alpha;
}
}  // namespace

TEST(PolarizationSplit, ReducedMatricesAreBlocksOfFullM) {
  const int n = 7;
  const fg::QBlocks Q = inPlaneQ(n, 3);
  const VecC alpha = inPlaneAlpha(n), gamma = VecC::Zero(n);
  const MatC M = fg::buildPropagationMatrix(Q, alpha, gamma);
  for (const Polarization p : {Polarization::TE, Polarization::TM}) {
    const MatC Mp = fg::buildPropagationMatrix(Q, alpha, p);
    EXPECT_LT((Mp - fieldBlocks(M, p)).cwiseAbs().maxCoeff(), 1e-12 * M.cwiseAbs().maxCoeff());
  }
  // the other blocks couple nothing (exact decoupling)
  const auto te = fg::polarizationBlocks(Polarization::TE),
             tm = fg::polarizationBlocks(Polarization::TM);
  for (int i : te.field)
    for (int j : tm.field) {
      EXPECT_EQ(M.block(i * n, j * n, n, n).cwiseAbs().maxCoeff(), 0.0);
      EXPECT_EQ(M.block(j * n, i * n, n, n).cwiseAbs().maxCoeff(), 0.0);
    }
}

TEST(PolarizationSplit, ReducedBasisIsEigenbasis) {
  const int n = 6;
  const cplx eps(2.1, 0.3);
  const VecC alpha = inPlaneAlpha(n), gamma = VecC::Zero(n);
  const fg::ModeBasis full = fg::isotropicModes(eps, alpha, gamma);
  for (const Polarization p : {Polarization::TE, Polarization::TM}) {
    const fg::ModeBasis mb = fg::polarizationBasis(full, p);
    const MatC Mp = fg::buildPropagationMatrix(fg::homogeneousIsotropicQ(eps, n), alpha, p);
    EXPECT_LT((Mp * mb.Psi - mb.Psi * mb.lambda.asDiagonal()).cwiseAbs().maxCoeff(), 1e-12);
    EXPECT_LT((mb.PsiInv * mb.Psi - MatC::Identity(2 * n, 2 * n)).cwiseAbs().maxCoeff(), 1e-12);
  }
}

TEST(PolarizationSplit, GeneralTMBlockModes) {
  // Q_xy != 0: the TM matrix is not block anti-diagonal -> classified eigenproblem
  const int n = 6;
  const fg::QBlocks Q = inPlaneQ(n, 5);
  const MatC M = fg::buildPropagationMatrix(Q, inPlaneAlpha(n), Polarization::TM);
  ASSERT_FALSE(fg::isBlockAntiDiagonal(M));
  const fg::LayerModes m = fg::layerModesPolarization(M, Polarization::TM);
  ASSERT_EQ(m.Wu.rows(), 2 * n);
  ASSERT_EQ(m.lu.size(), n);
  EXPECT_LT((M * m.Wu - m.Wu * m.lu.asDiagonal()).cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT((M * m.Wd - m.Wd * m.ld.asDiagonal()).cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LE(m.lu.real().maxCoeff(), m.ld.real().minCoeff());
}

namespace {
const char* kBase = R"({
  "dimension": "1D", "solver": "SOLVER", "truncation": {"Nx": 8},
  "period_um": {"x": 0.55},
  "incidence": {"lambda_um": 0.63, "phi_deg": 17, "psi_deg": PSI, "delta_deg": DELTA},
  "rk": {"tol": 1e-9}, "rcwa": {"slices": 12},
  "materials": {"air": {"n": 1.0}, "glass": {"n": [1.5, 0.01]}, "al": {"n": [1.2, 7.0]}},
  "superstrate": "air", "substrate": "glass",
  "layers": [
    {"type": "homogeneous", "thickness_um": 0.07, "material": "glass"},
    {"type": "profile1d", "shape": "sinus", "depth_um": 0.12, "below": "al", "above": "air"},
    {"type": "lamellar1d", "thickness_um": 0.15, "breakpoints_um": [0.0, 0.2, 0.55], "materials": ["glass", "air"]}
  ]
})";

fg::Result run(const std::string& solver, const std::string& psi, bool split,
               const std::string& delta = "35") {
  std::string s = kBase;
  s.replace(s.find("SOLVER"), 6, solver);
  s.replace(s.find("PSI"), 3, psi);
  s.replace(s.find("DELTA"), 5, delta);
  auto c = fg::io::parseConfig(s);
  c.solver.splitPolarizations = split;
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
}

void expectSame(const fg::Result& a, const fg::Result& b, double tol) {
  ASSERT_EQ(a.reflected.size(), b.reflected.size());
  ASSERT_EQ(a.transmitted.size(), b.transmitted.size());
  for (std::size_t k = 0; k < a.reflected.size(); ++k) {
    EXPECT_EQ(a.reflected[k].n, b.reflected[k].n);
    EXPECT_NEAR(a.reflected[k].efficiency, b.reflected[k].efficiency, tol);
    EXPECT_NEAR(a.reflected[k].eff_s, b.reflected[k].eff_s, tol);
    EXPECT_NEAR(a.reflected[k].eff_p, b.reflected[k].eff_p, tol);
    EXPECT_NEAR(std::abs(a.reflected[k].amp_p - b.reflected[k].amp_p), 0.0, tol);
  }
  for (std::size_t k = 0; k < a.transmitted.size(); ++k) {
    EXPECT_NEAR(a.transmitted[k].efficiency, b.transmitted[k].efficiency, tol);
    EXPECT_NEAR(std::abs(a.transmitted[k].amp_s - b.transmitted[k].amp_s), 0.0, tol);
  }
  EXPECT_NEAR(a.absorption, b.absorption, tol);
}
}  // namespace

TEST(PolarizationSplit, RCWASplitEqualsFull) {
  const auto split = run("rcwa", "0", true), full = run("rcwa", "0", false);
  expectSame(split, full, 1e-10);
  EXPECT_GT(split.reflected.size(), 1u);
}

TEST(PolarizationSplit, PureTEAndTM) {
  // delta = 0 / 90 deg: only one polarization is excited, the other one is skipped
  for (const char* delta : {"0", "90"}) {
    const auto split = run("rcwa", "0", true, delta), full = run("rcwa", "0", false, delta);
    expectSame(split, full, 1e-10);
    double es = 0, ep = 0;
    for (const auto& o : split.reflected) {
      es += o.eff_s;
      ep += o.eff_p;
    }
    EXPECT_GT(es + ep, 0.01);
    EXPECT_LT(std::string(delta) == "0" ? ep : es, 1e-12) << "delta = " << delta;
  }
}

TEST(PolarizationSplit, RKSplitEqualsFull) {
  const auto split = run("rk", "0", true), full = run("rk", "0", false);
  expectSame(split, full, 1e-8);  // adaptive steps may differ between the systems
  EXPECT_LT(std::abs(split.diag.energyBalance + split.absorption - 1.0), 1e-12);
}

TEST(PolarizationSplit, ConicalUsesFullSystem) {
  // gamma != 0: no decoupling, the option must not change anything
  expectSame(run("rcwa", "20", true), run("rcwa", "20", false), 0.0);
}
