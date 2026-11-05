// Resident RCWA chain (Plan section "M15", step 5): slice modes, slice S-matrices and star products
// on DeviceMatrix handles. Host reference of the resident operations = host path to rounding (runs
// everywhere); CUDA backend = host path on 1e-10 (GPU only, skipped otherwise).
#include <gtest/gtest.h>

#include <cstdlib>
#include <optional>

#include "fg/linalg/backend.hpp"
#include "fg/physics/modes.hpp"
#include "fg/solvers/resident.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;

namespace {

const int kNf = 5;

MatC random(Eigen::Index r, Eigen::Index c, unsigned seed) {
  std::srand(seed);
  return MatC::Random(r, c);
}

/// Propagation-matrix shapes of the three mode paths (4 n_f rows).
enum class Shape { Reduced, Decoupled, General };
MatC matrix(Shape s, unsigned seed) {
  const int n = kNf, N = 4 * kNf;
  MatC M = MatC::Zero(N, N);
  if (s == Shape::Reduced) {  // [[0, P], [R, 0]] (isotropic layers)
    M.topRightCorner(2 * n, 2 * n) =
        random(2 * n, 2 * n, seed) + 2.0 * MatC::Identity(2 * n, 2 * n);
    M.bottomLeftCorner(2 * n, 2 * n) =
        random(2 * n, 2 * n, seed + 1) - 3.0 * MatC::Identity(2 * n, 2 * n);
  } else if (s == Shape::Decoupled) {  // (Ex, Hz) and (Ez, Hx) do not couple
    const int sets[2][2] = {{0, 3}, {1, 2}};
    unsigned k = seed;
    for (const auto& set : sets)
      for (int a : set)
        for (int b : set)
          M.block(a * n, b * n, n, n) =
              random(n, n, k++) + (a == b ? 0.5 : 0.0) * MatC::Identity(n, n);
  } else {
    M = random(N, N, seed) + cplx(0.3, 0.1) * MatC::Identity(N, N);
  }
  return M;
}

fg::ModeBasis reference(bool inPlane = false) {
  VecC alpha(kNf), gamma(kNf);
  for (int r = 0; r < kNf; ++r) {
    alpha(r) = 0.2 + 0.7 * (r - kNf / 2);
    gamma(r) = inPlane ? 0.0 : 0.1 - 0.3 * (r - kNf / 2);
  }
  return fg::isotropicModes(1.0, alpha, gamma);
}

double rel(const fg::SMatrix& a, const fg::SMatrix& b) {
  const MatC fa = a.full(), fb = b.full();
  return (fa - fb).cwiseAbs().maxCoeff() / fa.cwiseAbs().maxCoeff();
}

/// Slice S-matrices of the three shapes and their star products, resident (backend d) vs host.
void compare(fg::linalg::Backend& d, double tol) {
  const fg::ModeBasis ref = reference();
  const fg::linalg::DeviceMatrix psiInv = d.upload(ref.PsiInv);
  std::vector<fg::SMatrix> host;
  std::vector<fg::resident::SMatrixD> dev;
  unsigned seed = 10;
  for (Shape s : {Shape::Reduced, Shape::Decoupled, Shape::General}) {
    const MatC M = matrix(s, seed += 10);
    const double h = 0.37;
    host.push_back(fg::stableLayerSMatrix(fg::layerModesAuto(M), h, ref));
    dev.push_back(
        fg::resident::layerSMatrix(d, fg::resident::layerModes(d, M, std::nullopt), h, psiInv));
    EXPECT_LT(rel(host.back(), fg::resident::download(d, dev.back())), tol)
        << "shape " << static_cast<int>(s);
  }
  // one TE/TM block (2 n_f rows, split polarizations)
  {
    const MatC Mb = matrix(Shape::General, 77).topLeftCorner(2 * kNf, 2 * kNf);
    const fg::ModeBasis rb = fg::polarizationBasis(reference(true), fg::Polarization::TM);
    const fg::SMatrix a =
        fg::stableLayerSMatrix(fg::layerModesPolarization(Mb, fg::Polarization::TM), 0.2, rb);
    const fg::SMatrix b = fg::resident::download(
        d, fg::resident::layerSMatrix(d, fg::resident::layerModes(d, Mb, fg::Polarization::TM), 0.2,
                                      d.upload(rb.PsiInv)));
    EXPECT_LT(rel(a, b), tol) << "polarization block";
  }
  fg::SMatrix h = host[0];
  fg::resident::SMatrixD g = dev[0];
  for (std::size_t k = 1; k < host.size(); ++k) {
    h = fg::star(h, host[k]);
    g = fg::resident::star(d, g, dev[k]);
  }
  EXPECT_LT(rel(h, fg::resident::download(d, g)), tol) << "star chain";
}

}  // namespace

TEST(Resident, HostReferenceEqualsHostPath) {
  const auto cpu = fg::linalg::makeCpuBackend();
  EXPECT_FALSE(cpu->resident());  // the solver never takes the resident path on the CPU
  compare(*cpu, 1e-12);
}

TEST(Resident, CudaEqualsHostPath) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  const auto gpu = fg::linalg::makeCudaBackend(0);
  ASSERT_TRUE(gpu->resident());
  compare(*gpu, 1e-10);
}

// Whole solver: resident chain on the GPU (several streams, tree) = CPU host chain on 1e-10.
TEST(Resident, CudaSolverEqualsCpu) {
  if (!fg::linalg::cudaAvailable()) GTEST_SKIP() << "no CUDA device";
  fg::Stack st;
  st.dx_um = st.dz_um = 0.6;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["si"] = fg::Material::fromIndex("si", cplx(3.88, 0.02));
  st.superstrate = "air";
  st.substrate = "si";
  fg::Profile2DLayer p;
  p.shape = fg::Profile2DLayer::Shape::EggCrate;
  p.depth_um = 0.15;
  p.below = "si";
  p.above = "air";
  p.slices = 6;
  fg::LayerEntry le;
  le.spec = p;
  st.layers = {le};
  fg::SolverOptions o;
  o.Nx = o.Nz = 3;
  o.backend = "cpu";
  const fg::Wavelength wl{0.6328, 0.17, 0.35, 0.79};
  const auto a = fg::makeSolver(o)->solve(st, wl);
  o.backend = "cuda";
  o.backendThreshold = 0;
  for (bool tree : {false, true}) {
    o.treeReduction = tree;
    o.threadsSlices = 3;
    const auto b = fg::makeSolver(o)->solve(st, wl);
    ASSERT_EQ(a.reflected.size(), b.reflected.size());
    for (std::size_t k = 0; k < a.reflected.size(); ++k)
      EXPECT_NEAR(a.reflected[k].efficiency, b.reflected[k].efficiency, 1e-10) << tree;
    EXPECT_NEAR(a.absorption, b.absorption, 1e-10) << tree;
  }
}
