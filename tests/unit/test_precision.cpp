// Working precision at the backend (Plan section "M16: Waehlbare Genauigkeit"): backend operations
// in fp32 against fp64 (well-conditioned random matrices, n = 50..400): relative error <= 1e-5;
// iterative refinement (mixed) reaches FP64 accuracy (1e-12) for cond <= 1e6. CPU always, CUDA if
// available. Solver: settings, diagnostics, results close to fp64, FG_PRECISION.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <memory>
#include <random>

#include "fg/core/error.hpp"
#include "fg/io/config.hpp"
#include "fg/linalg/backend.hpp"

using fg::cplx;
using fg::MatC;
using fg::linalg::Precision;

namespace {

void setEnv(const char* key, const char* value) {
#ifdef _WIN32
  _putenv_s(key, value ? value : "");
#else
  if (value)
    setenv(key, value, 1);
  else
    unsetenv(key);
#endif
}

MatC random(Eigen::Index n, Eigen::Index m, unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> nd;
  MatC A(n, m);
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j < m; ++j) A(i, j) = cplx(nd(rng), nd(rng));
  return A;
}

/// Matrix with condition number cond (singular values logarithmically spaced).
MatC conditioned(Eigen::Index n, double cond, unsigned seed) {
  const MatC U = random(n, n, seed).householderQr().householderQ();
  const MatC V = random(n, n, seed + 1).householderQr().householderQ();
  fg::VecC s(n);
  for (Eigen::Index k = 0; k < n; ++k)
    s(k) = std::pow(cond, -static_cast<double>(k) / static_cast<double>(n - 1));
  return U * s.asDiagonal() * V.adjoint();
}

double rel(const MatC& a, const MatC& b) {
  return (a - b).norm() / b.norm();
}

std::vector<std::unique_ptr<fg::linalg::Backend>> backends() {
  std::vector<std::unique_ptr<fg::linalg::Backend>> v;
  v.push_back(fg::linalg::makeCpuBackend());
  if (fg::linalg::cudaAvailable()) v.push_back(fg::linalg::makeCudaBackend(0));
  return v;
}

}  // namespace

TEST(Precision, Names) {
  for (Precision p : {Precision::FP64, Precision::Mixed, Precision::FP32})
    EXPECT_EQ(fg::linalg::parsePrecision(fg::linalg::precisionName(p)), p);
  EXPECT_THROW(fg::linalg::parsePrecision("tf32"), fg::Error);
  EXPECT_EQ(fg::linalg::unitRoundoff(Precision::FP32), std::ldexp(1.0, -24));
}

TEST(Precision, Fp32OperationsAgainstFp64) {
  for (auto& b : backends())
    for (Eigen::Index n : {50, 150, 400}) {
      const MatC A = conditioned(n, 2.0, 1), B = random(n, 7, 2);
      b->setPrecision(Precision::FP64);
      const MatC G64 = b->gemm(A, B), X64 = b->solve(A, B);
      b->setPrecision(Precision::FP32);
      EXPECT_LT(rel(b->gemm(A, B), G64), 1e-5) << b->name() << " n = " << n;
      const MatC X32 = b->solve(A, B);  // backward error (forward error ~ n eps32 cond)
      EXPECT_LT((A * X32 - B).norm() / (A.norm() * X32.norm()), 1e-5) << b->name() << " n = " << n;
      EXPECT_LT(rel(X32, X64), 1e-4) << b->name() << " n = " << n;
      // eigenpairs: residual relative to ||A||, eigenvalues against FP64
      const MatC E = random(n, n, 3) / std::sqrt(static_cast<double>(n));
      const auto e32 = b->eig(E);
      b->setPrecision(Precision::FP64);
      const auto e64 = b->eig(E);
      double res = 0, dist = 0;
      for (Eigen::Index k = 0; k < n; ++k) {
        res = std::max(res, (E * e32.vectors.col(k) - e32.values(k) * e32.vectors.col(k)).norm());
        double best = 1e300;  // nearest FP64 eigenvalue
        for (Eigen::Index j = 0; j < n; ++j)
          best = std::min(best, std::abs(e32.values(k) - e64.values(j)));
        dist = std::max(dist, best);
        EXPECT_NEAR(e32.vectors.col(k).norm(), 1.0, 1e-6);
      }
      const double nE = E.norm();
      EXPECT_LT(res / nE, 1e-5) << b->name() << " n = " << n;
      EXPECT_LT(dist / nE, 1e-5) << b->name() << " n = " << n;
    }
}

TEST(Precision, MixedRefinementReachesFp64) {
  for (auto& b : backends())
    for (double cond : {1e2, 1e4, 1e6}) {
      const Eigen::Index n = 200;
      const MatC A = conditioned(n, cond, 5), Xtrue = random(n, 3, 6), B = A * Xtrue;
      b->setPrecision(Precision::FP64);
      const MatC X64 = b->solve(A, B);
      fg::linalg::RefinementCounter c;
      fg::linalg::setRefinementCounter(&c);
      b->setPrecision(Precision::Mixed);
      const MatC Xm = b->solve(A, B);
      fg::linalg::setRefinementCounter(nullptr);
      // as accurate as the FP64 solve (forward error ~ cond eps64)
      EXPECT_LT(rel(Xm, Xtrue), std::max(1e-12, 4 * rel(X64, Xtrue)))
          << b->name() << " cond = " << cond;
      EXPECT_GT(c.steps.load() + c.fallbacks.load(), 0) << b->name();
      b->setPrecision(Precision::FP32);
      EXPECT_GT(rel(b->solve(A, B), Xtrue), 1e-10) << "fp32 alone is not refined";
    }
}

TEST(Precision, Fp32ResidentChainOnGpu) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  auto b = fg::linalg::makeCudaBackend(0);
  b->setPrecision(Precision::FP32);
  EXPECT_TRUE(b->resident());
  const MatC A = conditioned(120, 10.0, 7), B = random(120, 120, 8);
  const auto dA = b->upload(A), dB = b->upload(B);
  EXPECT_EQ(dA.elementSize(), sizeof(std::complex<float>));
  EXPECT_LT(rel(b->download(b->multiply(dA, dB)), A * B), 1e-5);
  EXPECT_LT(rel(b->download(b->luSolve(dA, dB)), A.partialPivLu().solve(B)), 1e-5);
  EXPECT_LT(rel(b->download(b->add(1.0, dA, true, 2.0, dB, false)), MatC(A.transpose()) + 2.0 * B),
            1e-6);
  b->setPrecision(Precision::Mixed);
  EXPECT_FALSE(b->resident());  // mixed uses the host API
}

namespace {
fg::io::RunConfig cfg(const std::string& extra) {
  return fg::io::parseConfig(
      R"({"dimension": "2D", "solver": "rcwa", "truncation": {"Nx": 3, "Nz": 3},
    "period_um": {"x": 0.6, "z": 0.6}, "incidence": {"lambda_um": 0.633, "phi_deg": 10, "psi_deg": 20, "delta_deg": 45},
    "rcwa": {"slices": 6}, "backend": "cpu",
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
    "superstrate": "air", "substrate": "glass",
    "layers": [{"type": "profile2d", "shape": "egg_crate", "depth_um": 0.15, "below": "glass", "above": "air"}])" +
      extra + "}");
}
double maxDiff(const fg::Result& a, const fg::Result& b) {
  double d = 0;
  for (std::size_t k = 0; k < a.reflected.size(); ++k)
    d = std::max(d, std::abs(a.reflected[k].efficiency - b.reflected[k].efficiency));
  for (std::size_t k = 0; k < a.transmitted.size(); ++k)
    d = std::max(d, std::abs(a.transmitted[k].efficiency - b.transmitted[k].efficiency));
  return d;
}
fg::Result run(const fg::io::RunConfig& c) {
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
}
}  // namespace

TEST(Precision, SolverModesOnTheCpu) {
  const auto ref = run(cfg(""));
  EXPECT_EQ(ref.diag.precision, "fp64");  // auto on the CPU
  const auto mixed = run(cfg(R"(, "precision": "mixed")"));
  EXPECT_EQ(mixed.diag.precision, "mixed");
  EXPECT_LT(maxDiff(mixed, ref), 1e-5);  // eigenproblems in C dominate (Plan "M16", study)
  EXPECT_GT(mixed.diag.refinementSteps + mixed.diag.refinementFallbacks, 0);
  const auto fp32 = run(cfg(R"(, "precision": "fp32")"));
  EXPECT_EQ(fp32.diag.precision, "fp32");
  EXPECT_LT(maxDiff(fp32, ref), 1e-4);
  EXPECT_FALSE(fp32.diag.precisionFallback);  // lossless, energy conserved to ~1e-6
  EXPECT_NEAR(fp32.diag.energyBalance, 1.0, 1e-4);
  // FG_PRECISION overrides the configuration
  setEnv("FG_PRECISION", "fp32");
  const auto env = run(cfg(R"(, "precision": "fp64")"));
  setEnv("FG_PRECISION", nullptr);
  EXPECT_EQ(env.diag.precision, "fp32");
  EXPECT_THROW(cfg(R"(, "precision": "half")"), fg::Error);
}

TEST(Precision, AutoUsesTheDevicePolicyOnTheGpu) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  auto c = cfg(R"(, "backend": "cuda")");
  const auto r = run(c);
  EXPECT_EQ(r.diag.precision, fg::linalg::precisionName(fg::linalg::gpuAutoPrecision()));
  const auto ref = run(cfg(""));
  EXPECT_LT(maxDiff(r, ref), 1e-4);
}
