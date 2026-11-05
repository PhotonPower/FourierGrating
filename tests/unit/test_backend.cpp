// Linear algebra backends (Plan section "CUDA-Backend", M10): dispatch, GPU = CPU on 1e-10
// (skipped without a CUDA device), parallel wavelength scan.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::MatC;

namespace {
/// Portable setenv/unsetenv (MinGW/UCRT has no setenv); value nullptr removes the variable.
void setEnv(const char* key, const char* value) {
#ifdef _WIN32
  _putenv_s(key, value ? value : "");
#else
  if (value) setenv(key, value, 1);
  else unsetenv(key);
#endif
}

MatC randomMatrix(Eigen::Index n, Eigen::Index m, unsigned seed) {
  std::srand(seed);
  return MatC::Random(n, m) + 2.0 * MatC::Identity(n, m);
}

fg::Stack holes() {
  fg::Stack st;
  st.dx_um = st.dz_um = 0.9;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["ag"] = fg::Material::fromIndex("ag", cplx(0.1, 8.94));
  st.superstrate = "air";
  st.substrate = "glass";
  fg::Polygon2DLayer g;
  g.thickness_um = 0.2;
  g.background = "ag";
  g.inclusions.push_back({{{0.325, 0.325}, {0.575, 0.325}, {0.575, 0.575}, {0.325, 0.575}}, "air"});
  fg::LayerEntry le;
  le.spec = g;
  st.layers = {le};
  return st;
}
}  // namespace

TEST(Backend, DispatchSettings) {
  EXPECT_THROW(fg::linalg::setPreference("gpu"), fg::Error);
  fg::linalg::setPreference("cpu");
  EXPECT_FALSE(fg::linalg::offload(100000));
  fg::linalg::setPreference("auto", 1000);
  EXPECT_EQ(fg::linalg::offload(5000), fg::linalg::cudaAvailable());
  EXPECT_FALSE(fg::linalg::offload(10));
  if (!fg::linalg::cudaAvailable()) EXPECT_THROW(fg::linalg::setPreference("cuda"), fg::Error);
  fg::linalg::setPreference("auto");
}

TEST(Backend, EnvironmentOverride) {
  setEnv("FG_BACKEND", "cpu");
  fg::linalg::setPreference("auto", 0);
  EXPECT_FALSE(fg::linalg::offload(100000));
  setEnv("FG_BACKEND", "gpu");
  EXPECT_THROW(fg::linalg::setPreference("auto"), fg::Error);
  setEnv("FG_BACKEND", nullptr);
  setEnv("FG_BACKEND_THRESHOLD", "5");
  fg::linalg::setPreference("auto", 1000);
  EXPECT_EQ(fg::linalg::offload(10), fg::linalg::cudaAvailable());
  setEnv("FG_BACKEND_THRESHOLD", nullptr);
  fg::linalg::setPreference("auto");
}

TEST(Backend, CpuBackendOperations) {
  auto cpu = fg::linalg::makeCpuBackend();
  const MatC A = randomMatrix(40, 40, 1), B = randomMatrix(40, 7, 2);
  EXPECT_LT((A * cpu->solve(A, B) - B).norm(), 1e-12 * B.norm());
  EXPECT_LT((cpu->gemm(A, B) - A * B).norm(), 1e-13 * (A * B).norm());
  const auto e = cpu->eig(A);
  for (Eigen::Index k = 0; k < 40; ++k)
    EXPECT_LT((A * e.vectors.col(k) - e.values(k) * e.vectors.col(k)).norm(), 1e-11 * A.norm());
}

TEST(Backend, CudaEqualsCpu) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  auto cpu = fg::linalg::makeCpuBackend();
  auto gpu = fg::linalg::makeCudaBackend(0);
  const Eigen::Index n = 300;
  const MatC A = randomMatrix(n, n, 3), B = randomMatrix(n, 50, 4);
  EXPECT_LT((gpu->gemm(A, B) - cpu->gemm(A, B)).norm(), 1e-12 * (A * B).norm());
  EXPECT_LT((gpu->solve(A, B) - cpu->solve(A, B)).norm(), 1e-10 * cpu->solve(A, B).norm());
  const auto eg = gpu->eig(A), ec = cpu->eig(A);
  std::vector<cplx> vg(eg.values.data(), eg.values.data() + n), vc(ec.values.data(), ec.values.data() + n);
  auto lex = [](cplx a, cplx b) { return a.real() != b.real() ? a.real() < b.real() : a.imag() < b.imag(); };
  std::sort(vg.begin(), vg.end(), lex);
  std::sort(vc.begin(), vc.end(), lex);
  for (Eigen::Index k = 0; k < n; ++k) {
    EXPECT_LT(std::abs(vg[static_cast<std::size_t>(k)] - vc[static_cast<std::size_t>(k)]), 1e-10 * A.norm());
    EXPECT_LT((A * eg.vectors.col(k) - eg.values(k) * eg.vectors.col(k)).norm(), 1e-10 * A.norm());
  }
}

// Full 2D RCWA (Ag hole array, N = 4: matrices up to 4 n_f = 324) on the device = CPU on 1e-10.
TEST(Backend, CudaSolverEqualsCpu) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  fg::SolverOptions o;
  o.Nx = o.Nz = 4;
  o.backend = "cpu";
  const fg::Wavelength wl{1.388, 0.2, 0.3, 0.4};
  const auto a = fg::makeSolver(o)->solve(holes(), wl);
  o.backend = "cuda";  // every offloadable operation on the device
  const auto b = fg::makeSolver(o)->solve(holes(), wl);
  ASSERT_EQ(a.transmitted.size(), b.transmitted.size());
  for (std::size_t k = 0; k < a.transmitted.size(); ++k)
    EXPECT_NEAR(a.transmitted[k].efficiency, b.transmitted[k].efficiency,
                1e-10 * std::max(1e-3, a.transmitted[k].efficiency));
  EXPECT_NEAR(a.absorption, b.absorption, 1e-10);
}

// Parallel scan (batching over wavelengths) = sequential solves.
TEST(Backend, ParallelScanEqualsSequential) {
  fg::SolverOptions o;
  o.Nx = o.Nz = 2;
  std::vector<fg::Wavelength> scan;
  for (int k = 0; k < 5; ++k) scan.push_back({1.35 + 0.02 * k, 0.1, 0.0, 0.0});
  const auto par = fg::solveScan(o, holes(), scan);
  auto s = fg::makeSolver(o);
  for (std::size_t k = 0; k < scan.size(); ++k) {
    const auto r = s->solve(holes(), scan[k]);
    ASSERT_EQ(r.transmitted.size(), par[k].transmitted.size());
    for (std::size_t j = 0; j < r.transmitted.size(); ++j)
      EXPECT_EQ(r.transmitted[j].efficiency, par[k].transmitted[j].efficiency);
  }
}
