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
bool cudaAvailable() { return false; }
#endif

std::unique_ptr<Backend> makeBackend(const std::string& preference, Eigen::Index problemSize) {
  FG_REQUIRE(preference == "auto" || preference == "cpu" || preference == "cuda",
             "backend must be \"auto\", \"cpu\" or \"cuda\"");
  if (preference == "cuda") return makeCudaBackend(0);
  if (preference == "auto" && problemSize >= 1000 && cudaAvailable()) return makeCudaBackend(0);
  return makeCpuBackend();
}

namespace {
struct DispatchState {
  std::string preference = "auto";
  Eigen::Index threshold = 1000;
  std::unique_ptr<Backend> dev;
  int available = -1;  // lazily probed
};
DispatchState& state() {
  thread_local DispatchState st;
  return st;
}
}  // namespace

void setPreference(const std::string& preference, Eigen::Index threshold) {
  FG_REQUIRE(preference == "auto" || preference == "cpu" || preference == "cuda",
             "backend must be \"auto\", \"cpu\" or \"cuda\"");
  if (preference == "cuda") FG_REQUIRE(cudaAvailable(), "CUDA backend requested but not available");
  state().preference = preference;
  state().threshold = threshold;
}

bool offload(Eigen::Index n) {
  DispatchState& st = state();
  if (st.preference == "cpu") return false;
  if (st.available < 0) st.available = cudaAvailable() ? 1 : 0;
  if (!st.available) return false;
  return st.preference == "cuda" || n >= st.threshold;
}

Backend& device() {
  DispatchState& st = state();
  if (!st.dev) st.dev = makeCudaBackend(0);
  return *st.dev;
}

}  // namespace fg::linalg
