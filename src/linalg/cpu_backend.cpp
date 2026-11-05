#include <algorithm>
#include <cstdlib>

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

// ---- host reference of the resident operations (M15, step 5) ----
namespace {
Eigen::Map<MatC> map(const DeviceMatrix& A) { return Eigen::Map<MatC>(A.data(), A.rows(), A.cols()); }
}  // namespace

DeviceMatrix Backend::hostMatrix(Eigen::Index rows, Eigen::Index cols) {
  std::shared_ptr<void> mem(new cplx[static_cast<std::size_t>(std::max<Eigen::Index>(rows * cols, 1))],
                            [](void* p) { delete[] static_cast<cplx*>(p); });
  return DeviceMatrix(std::move(mem), rows, cols);
}
DeviceMatrix Backend::upload(const MatC& A) {
  DeviceMatrix D = hostMatrix(A.rows(), A.cols());
  map(D) = A;
  return D;
}
MatC Backend::download(const DeviceMatrix& A) { return map(A); }
DeviceMatrix Backend::zeros(Eigen::Index rows, Eigen::Index cols) {
  DeviceMatrix D = hostMatrix(rows, cols);
  map(D).setZero();
  return D;
}
DeviceMatrix Backend::identity(Eigen::Index n) {
  DeviceMatrix D = hostMatrix(n, n);
  map(D).setIdentity();
  return D;
}
DeviceMatrix Backend::multiply(const DeviceMatrix& A, const DeviceMatrix& B, cplx alpha, const DeviceMatrix* C,
                               cplx beta) {
  FG_REQUIRE(A.cols() == B.rows() && (!C || (C->rows() == A.rows() && C->cols() == B.cols())), "gemm: size mismatch");
  DeviceMatrix D = hostMatrix(A.rows(), B.cols());
  map(D).noalias() = alpha * (map(A) * map(B));
  if (C) map(D) += beta * map(*C);
  return D;
}
DeviceMatrix Backend::luSolve(const DeviceMatrix& A, const DeviceMatrix& B) {
  FG_REQUIRE(A.rows() == A.cols() && A.rows() == B.rows(), "solve: size mismatch");
  DeviceMatrix D = hostMatrix(B.rows(), B.cols());
  map(D) = map(A).partialPivLu().solve(map(B));
  return D;
}
DeviceMatrix Backend::add(cplx alpha, const DeviceMatrix& A, bool transA, cplx beta, const DeviceMatrix& B, bool transB) {
  const MatC a = transA ? MatC(map(A).transpose()) : MatC(map(A));
  DeviceMatrix D = hostMatrix(a.rows(), a.cols());
  map(D) = alpha * a;
  if (!B.empty()) {
    const MatC b = transB ? MatC(map(B).transpose()) : MatC(map(B));
    FG_REQUIRE(b.rows() == a.rows() && b.cols() == a.cols(), "add: size mismatch");
    map(D) += beta * b;
  }
  return D;
}
DeviceMatrix Backend::scaleRows(const VecC& d, const DeviceMatrix& A) {
  FG_REQUIRE(d.size() == A.rows(), "scaleRows: size mismatch");
  DeviceMatrix D = hostMatrix(A.rows(), A.cols());
  map(D) = d.asDiagonal() * map(A);
  return D;
}
DeviceMatrix Backend::scaleCols(const DeviceMatrix& A, const VecC& d) {
  FG_REQUIRE(d.size() == A.cols(), "scaleCols: size mismatch");
  DeviceMatrix D = hostMatrix(A.rows(), A.cols());
  map(D) = map(A) * d.asDiagonal();
  return D;
}
DeviceMatrix Backend::block(const DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, Eigen::Index nr, Eigen::Index nc) {
  FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + nr <= A.rows() && c0 + nc <= A.cols(), "block: out of range");
  DeviceMatrix D = hostMatrix(nr, nc);
  map(D) = map(A).block(r0, c0, nr, nc);
  return D;
}
void Backend::setBlock(DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, const DeviceMatrix& B) {
  FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + B.rows() <= A.rows() && c0 + B.cols() <= A.cols(), "setBlock: out of range");
  map(A).block(r0, c0, B.rows(), B.cols()) = map(B);
}
DeviceMatrix Backend::columns(const DeviceMatrix& A, const std::vector<Eigen::Index>& idx) {
  DeviceMatrix D = hostMatrix(A.rows(), static_cast<Eigen::Index>(idx.size()));
  for (std::size_t k = 0; k < idx.size(); ++k) map(D).col(static_cast<Eigen::Index>(k)) = map(A).col(idx[k]);
  return D;
}
DeviceEig Backend::eigResident(const DeviceMatrix& A) {
  EigResult e = eig(download(A));
  return {std::move(e.values), upload(e.vectors)};
}

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

void setPreference(const std::string& requested, Eigen::Index threshold) {
  // environment override (whole test suite on the GPU: FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest)
  std::string preference = requested;
  if (const char* env = std::getenv("FG_BACKEND"); env && *env) preference = env;
  if (const char* env = std::getenv("FG_BACKEND_THRESHOLD"); env && *env) threshold = std::atol(env);
  FG_REQUIRE(preference == "auto" || preference == "cpu" || preference == "cuda",
             "backend must be \"auto\", \"cpu\" or \"cuda\" (also FG_BACKEND)");
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

int deviceStreams() {
  if (const char* env = std::getenv("FG_GPU_STREAMS"); env && *env && std::atoi(env) > 0) return std::atoi(env);
  return 4;
}

Backend* detail::currentDevice() { return state().dev.get(); }

}  // namespace fg::linalg
