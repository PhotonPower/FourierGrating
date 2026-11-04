#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"

namespace fg::linalg {

namespace {
class CpuBackend final : public Backend {
 public:
  std::string name() const override { return "cpu-eigen"; }
  MatC gemm(const MatC& A, const MatC& B) override { return A * B; }
  MatC solve(const MatC& A, const MatC& B) override { return A.partialPivLu().solve(B); }
  EigResult eig(const MatC& A) override {
    Eigen::ComplexEigenSolver<MatC> es(A, true);
    FG_REQUIRE(es.info() == Eigen::Success, "eigen decomposition failed");
    return {es.eigenvalues(), es.eigenvectors()};
  }
};
}  // namespace

std::unique_ptr<Backend> makeCpuBackend() { return std::make_unique<CpuBackend>(); }

#ifndef FG_HAVE_CUDA
std::unique_ptr<Backend> makeCudaBackend(int) {
  throw Error("FourierGrating was built without CUDA (configure with -DFG_ENABLE_CUDA=ON)");
}
#endif

std::unique_ptr<Backend> makeBackend(const std::string& preference, Eigen::Index problemSize) {
#ifdef FG_HAVE_CUDA
  if (preference == "cuda" || (preference == "auto" && problemSize >= 1000))
    return makeCudaBackend(0);
#else
  (void)problemSize;
  FG_REQUIRE(preference != "cuda", "CUDA backend requested but not built");
#endif
  return makeCpuBackend();
}

}  // namespace fg::linalg
