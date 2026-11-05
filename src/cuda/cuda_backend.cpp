// FourierGrating - CUDA backend (milestone M10, Plan section "CUDA-Backend").
// Host API of cuBLAS (Zgemm) and cuSOLVER (64-bit Xgetrf/Xgetrs, Xgeev; CUDA >= 12.6 Update 2);
// no device kernels, so a plain C++ compiler suffices. Column-major data (Eigen default).
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <vector>

#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"

namespace fg::linalg {

namespace {

#define FG_CUDA(call)                                                                                    \
  do {                                                                                                   \
    const cudaError_t e_ = (call);                                                                       \
    if (e_ != cudaSuccess) throw Error(std::string(#call) + ": " + cudaGetErrorString(e_));              \
  } while (0)
#define FG_CUBLAS(call)                                                                                  \
  do {                                                                                                   \
    const cublasStatus_t e_ = (call);                                                                    \
    if (e_ != CUBLAS_STATUS_SUCCESS)                                                                     \
      throw Error(std::string(#call) + ": cuBLAS status " + std::to_string(static_cast<int>(e_)));       \
  } while (0)
#define FG_CUSOLVER(call)                                                                                \
  do {                                                                                                   \
    const cusolverStatus_t e_ = (call);                                                                  \
    if (e_ != CUSOLVER_STATUS_SUCCESS)                                                                   \
      throw Error(std::string(#call) + ": cuSOLVER status " + std::to_string(static_cast<int>(e_)));     \
  } while (0)

/// RAII device memory (grows, never shrinks).
template <class T>
class DeviceBuffer {
 public:
  DeviceBuffer() = default;
  ~DeviceBuffer() { release(); }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  void resize(std::size_t n) {
    if (n <= n_) return;
    release();
    FG_CUDA(cudaMalloc(reinterpret_cast<void**>(&p_), n * sizeof(T)));
    n_ = n;
  }
  T* get() { return p_; }

 private:
  void release() {
    if (p_) cudaFree(p_);
    p_ = nullptr;
    n_ = 0;
  }
  T* p_ = nullptr;
  std::size_t n_ = 0;
};

/// RAII pinned host memory (Xgeev host workspace).
class PinnedBuffer {
 public:
  PinnedBuffer() = default;
  ~PinnedBuffer() {
    if (p_) cudaFreeHost(p_);
  }
  PinnedBuffer(const PinnedBuffer&) = delete;
  PinnedBuffer& operator=(const PinnedBuffer&) = delete;
  void* get(std::size_t bytes) {
    if (bytes > n_) {
      if (p_) cudaFreeHost(p_);
      p_ = nullptr;
      FG_CUDA(cudaMallocHost(&p_, bytes));
      n_ = bytes;
    }
    return p_;
  }

 private:
  void* p_ = nullptr;
  std::size_t n_ = 0;
};

static_assert(sizeof(cplx) == sizeof(cuDoubleComplex), "complex layout");

/// Device memory of the DeviceMatrix handles of one backend (M15, step 5): released buffers are kept
/// for reuse (no cudaFree, which would synchronize the device). A buffer released by another thread
/// is reused only after that thread's queued work has finished (see release()).
class MemoryPool {
 public:
  explicit MemoryPool(const Backend* owner) : owner_(owner) {}
  ~MemoryPool() {
    for (auto& [bytes, list] : free_)
      for (void* p : list) cudaFree(p);
  }
  MemoryPool(const MemoryPool&) = delete;
  MemoryPool& operator=(const MemoryPool&) = delete;
  void* acquire(std::size_t bytes) {
    {
      const std::lock_guard<std::mutex> lock(mu_);
      auto it = free_.find(bytes);
      if (it != free_.end() && !it->second.empty()) {
        void* p = it->second.back();
        it->second.pop_back();
        return p;
      }
    }
    void* p = nullptr;
    if (cudaMalloc(&p, bytes) != cudaSuccess) {
      (void)cudaGetLastError();
      trim();  // cached buffers of other sizes back to the device, then once more
      FG_CUDA(cudaMalloc(&p, bytes));
    }
    return p;
  }
  /// Work queued by the releasing thread on its own stream may still read the buffer: wait for it
  /// unless the buffer goes back to the same backend (its stream orders the reuse).
  void release(void* p, std::size_t bytes) {
    Backend* cur = detail::currentDevice();
    if (cur && cur != owner_) cur->synchronize();
    const std::lock_guard<std::mutex> lock(mu_);
    free_[bytes].push_back(p);
  }
  void trim() {
    FG_CUDA(cudaDeviceSynchronize());
    const std::lock_guard<std::mutex> lock(mu_);
    for (auto& [bytes, list] : free_)
      for (void* p : list) cudaFree(p);
    free_.clear();
  }

 private:
  const Backend* owner_;
  std::mutex mu_;
  std::map<std::size_t, std::vector<void*>> free_;
};

class CudaBackend final : public Backend {
 public:
  explicit CudaBackend(int device) {
    int count = 0;
    FG_CUDA(cudaGetDeviceCount(&count));
    FG_REQUIRE(device < count, "CUDA device index out of range");
    FG_CUDA(cudaSetDevice(device));
    cudaDeviceProp prop{};
    FG_CUDA(cudaGetDeviceProperties(&prop, device));
    std::ostringstream os;
    os << "cuda:" << prop.name << " (sm_" << prop.major << prop.minor << ", " << (prop.totalGlobalMem >> 20)
       << " MiB)";
    name_ = os.str();
    FG_CUDA(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
    FG_CUBLAS(cublasCreate(&blas_));
    FG_CUBLAS(cublasSetStream(blas_, stream_));
    FG_CUSOLVER(cusolverDnCreate(&solver_));
    FG_CUSOLVER(cusolverDnSetStream(solver_, stream_));
    FG_CUSOLVER(cusolverDnCreateParams(&params_));
  }
  ~CudaBackend() override {
    cusolverDnDestroyParams(params_);
    cusolverDnDestroy(solver_);
    cublasDestroy(blas_);
    cudaStreamDestroy(stream_);
  }
  CudaBackend(const CudaBackend&) = delete;
  CudaBackend& operator=(const CudaBackend&) = delete;
  std::string name() const override { return name_; }

  MatC gemm(const MatC& A, const MatC& B) override {
    FG_REQUIRE(A.cols() == B.rows(), "gemm: size mismatch");
    const std::int64_t m = A.rows(), k = A.cols(), n = B.cols();
    if (m == 0 || n == 0 || k == 0) return MatC::Zero(m, n);
    a_.resize(static_cast<std::size_t>(m * k));
    b_.resize(static_cast<std::size_t>(k * n));
    c_.resize(static_cast<std::size_t>(m * n));
    upload(a_.get(), A);
    upload(b_.get(), B);
    const cuDoubleComplex one{1.0, 0.0}, zero{0.0, 0.0};
    FG_CUBLAS(cublasZgemm(blas_, CUBLAS_OP_N, CUBLAS_OP_N, static_cast<int>(m), static_cast<int>(n),
                          static_cast<int>(k), &one, a_.get(), static_cast<int>(m), b_.get(), static_cast<int>(k),
                          &zero, c_.get(), static_cast<int>(m)));
    MatC C(m, n);
    download(C, c_.get());
    return C;
  }

  MatC solve(const MatC& A, const MatC& B) override {
    FG_REQUIRE(A.rows() == A.cols() && A.rows() == B.rows(), "solve: size mismatch");
    const std::int64_t n = A.rows(), nrhs = B.cols();
    if (n == 0 || nrhs == 0) return MatC::Zero(n, nrhs);
    a_.resize(static_cast<std::size_t>(n * n));
    b_.resize(static_cast<std::size_t>(n * nrhs));
    ipiv_.resize(static_cast<std::size_t>(n));
    info_.resize(1);
    upload(a_.get(), A);
    upload(b_.get(), B);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgetrf_bufferSize(solver_, params_, n, n, CUDA_C_64F, a_.get(), n, CUDA_C_64F, &dBytes,
                                            &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    std::vector<char> hwork(hBytes > 0 ? hBytes : 1);
    FG_CUSOLVER(cusolverDnXgetrf(solver_, params_, n, n, CUDA_C_64F, a_.get(), n, ipiv_.get(), CUDA_C_64F,
                                 work_.get(), dBytes, hwork.data(), hBytes, info_.get()));
    checkInfo("cusolverDnXgetrf");
    FG_CUSOLVER(cusolverDnXgetrs(solver_, params_, CUBLAS_OP_N, n, nrhs, CUDA_C_64F, a_.get(), n, ipiv_.get(),
                                 CUDA_C_64F, b_.get(), n, info_.get()));
    checkInfo("cusolverDnXgetrs");
    MatC X(n, nrhs);
    download(X, b_.get());
    return X;
  }

  EigResult eig(const MatC& A) override {
    FG_REQUIRE(A.rows() == A.cols(), "eig: matrix must be square");
    const std::int64_t n = A.rows();
    EigResult r;
    if (n == 0) return r;
    a_.resize(static_cast<std::size_t>(n * n));
    c_.resize(static_cast<std::size_t>(n * n));  // right eigenvectors
    w_.resize(static_cast<std::size_t>(n));
    info_.resize(1);
    upload(a_.get(), A);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgeev_bufferSize(solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR, n,
                                           CUDA_C_64F, a_.get(), n, CUDA_C_64F, w_.get(), CUDA_C_64F, nullptr, n,
                                           CUDA_C_64F, c_.get(), n, CUDA_C_64F, &dBytes, &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    FG_CUSOLVER(cusolverDnXgeev(solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR, n, CUDA_C_64F,
                                a_.get(), n, CUDA_C_64F, w_.get(), CUDA_C_64F, nullptr, n, CUDA_C_64F, c_.get(), n,
                                CUDA_C_64F, work_.get(), dBytes, pinned_.get(hBytes > 0 ? hBytes : 1), hBytes,
                                info_.get()));
    checkInfo("cusolverDnXgeev");
    r.values.resize(n);
    r.vectors.resize(n, n);
    FG_CUDA(cudaMemcpyAsync(r.values.data(), w_.get(), static_cast<std::size_t>(n) * sizeof(cplx),
                            cudaMemcpyDeviceToHost, stream_));
    download(r.vectors, c_.get());
    for (Eigen::Index k = 0; k < n; ++k) r.vectors.col(k).normalize();  // unit norm as Eigen
    return r;
  }

  // ---- resident operations (M15, step 5): everything on this backend's stream ----
  bool resident() const override { return true; }
  void synchronize() override { FG_CUDA(cudaStreamSynchronize(stream_)); }

  DeviceMatrix upload(const MatC& A) override {
    DeviceMatrix D = alloc(A.rows(), A.cols());
    if (A.size())
      FG_CUDA(cudaMemcpyAsync(D.data(), A.data(), bytes(D), cudaMemcpyHostToDevice, stream_));
    synchronize();  // A may be a temporary of the caller
    return D;
  }
  MatC download(const DeviceMatrix& A) override {
    MatC M(A.rows(), A.cols());
    if (A.size()) FG_CUDA(cudaMemcpyAsync(M.data(), A.data(), bytes(A), cudaMemcpyDeviceToHost, stream_));
    synchronize();
    return M;
  }
  DeviceMatrix zeros(Eigen::Index rows, Eigen::Index cols) override {
    DeviceMatrix D = alloc(rows, cols);
    if (D.size()) FG_CUDA(cudaMemsetAsync(D.data(), 0, bytes(D), stream_));
    return D;
  }
  DeviceMatrix identity(Eigen::Index n) override {
    DeviceMatrix D = zeros(n, n);
    if (n == 0) return D;
    const std::vector<cplx> ones(static_cast<std::size_t>(n), cplx(1.0, 0.0));
    // the diagonal: n elements with pitch (n + 1) elements
    FG_CUDA(cudaMemcpy2DAsync(D.data(), static_cast<std::size_t>(n + 1) * sizeof(cplx), ones.data(), sizeof(cplx),
                              sizeof(cplx), static_cast<std::size_t>(n), cudaMemcpyHostToDevice, stream_));
    synchronize();  // ones is a local host buffer
    return D;
  }
  DeviceMatrix multiply(const DeviceMatrix& A, const DeviceMatrix& B, cplx alpha, const DeviceMatrix* C,
                        cplx beta) override {
    FG_REQUIRE(A.cols() == B.rows() && (!C || (C->rows() == A.rows() && C->cols() == B.cols())),
               "multiply: size mismatch");
    DeviceMatrix D = alloc(A.rows(), B.cols());
    if (D.empty()) return D;
    cuDoubleComplex b{0.0, 0.0};
    if (C) {
      FG_CUDA(cudaMemcpyAsync(D.data(), C->data(), bytes(D), cudaMemcpyDeviceToDevice, stream_));
      b = cz(beta);
    }
    if (A.cols() == 0) {
      if (!C) FG_CUDA(cudaMemsetAsync(D.data(), 0, bytes(D), stream_));
      return D;
    }
    const cuDoubleComplex a = cz(alpha);
    FG_CUBLAS(cublasZgemm(blas_, CUBLAS_OP_N, CUBLAS_OP_N, i32(A.rows()), i32(B.cols()), i32(A.cols()), &a, z(A),
                          i32(A.rows()), z(B), i32(B.rows()), &b, z(D), i32(D.rows())));
    return D;
  }
  DeviceMatrix luSolve(const DeviceMatrix& A, const DeviceMatrix& B) override {
    FG_REQUIRE(A.rows() == A.cols() && A.rows() == B.rows(), "luSolve: size mismatch");
    const std::int64_t n = A.rows(), nrhs = B.cols();
    DeviceMatrix X = alloc(n, nrhs);
    if (n == 0 || nrhs == 0) return X;
    DeviceMatrix LU = alloc(n, n);
    FG_CUDA(cudaMemcpyAsync(LU.data(), A.data(), bytes(A), cudaMemcpyDeviceToDevice, stream_));
    FG_CUDA(cudaMemcpyAsync(X.data(), B.data(), bytes(B), cudaMemcpyDeviceToDevice, stream_));
    ipiv_.resize(static_cast<std::size_t>(n));
    info_.resize(1);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgetrf_bufferSize(solver_, params_, n, n, CUDA_C_64F, LU.data(), n, CUDA_C_64F, &dBytes,
                                            &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    void* hwork = pinned_.get(hBytes > 0 ? hBytes : 1);
    FG_CUSOLVER(cusolverDnXgetrf(solver_, params_, n, n, CUDA_C_64F, LU.data(), n, ipiv_.get(), CUDA_C_64F,
                                 work_.get(), dBytes, hwork, hBytes, info_.get()));
    checkInfo("cusolverDnXgetrf");
    FG_CUSOLVER(cusolverDnXgetrs(solver_, params_, CUBLAS_OP_N, n, nrhs, CUDA_C_64F, LU.data(), n, ipiv_.get(),
                                 CUDA_C_64F, X.data(), n, info_.get()));
    checkInfo("cusolverDnXgetrs");
    return X;
  }
  DeviceMatrix add(cplx alpha, const DeviceMatrix& A, bool transA, cplx beta, const DeviceMatrix& B,
                   bool transB) override {
    const Eigen::Index m = transA ? A.cols() : A.rows(), n = transA ? A.rows() : A.cols();
    DeviceMatrix D = alloc(m, n);
    if (D.empty()) return D;
    const bool useB = !B.empty();
    if (useB) FG_REQUIRE((transB ? B.cols() : B.rows()) == m && (transB ? B.rows() : B.cols()) == n, "add: size mismatch");
    const cuDoubleComplex a = cz(alpha), b = useB ? cz(beta) : cuDoubleComplex{0.0, 0.0};
    const DeviceMatrix& Bm = useB ? B : A;  // not read for beta = 0
    const bool tb = useB ? transB : transA;
    FG_CUBLAS(cublasZgeam(blas_, transA ? CUBLAS_OP_T : CUBLAS_OP_N, tb ? CUBLAS_OP_T : CUBLAS_OP_N, i32(m), i32(n),
                          &a, z(A), i32(A.rows()), &b, z(Bm), i32(Bm.rows()), z(D), i32(m)));
    return D;
  }
  DeviceMatrix scaleRows(const VecC& d, const DeviceMatrix& A) override { return dgmm(CUBLAS_SIDE_LEFT, A, d); }
  DeviceMatrix scaleCols(const DeviceMatrix& A, const VecC& d) override { return dgmm(CUBLAS_SIDE_RIGHT, A, d); }
  DeviceMatrix block(const DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, Eigen::Index nr, Eigen::Index nc) override {
    FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + nr <= A.rows() && c0 + nc <= A.cols(), "block: out of range");
    DeviceMatrix D = alloc(nr, nc);
    if (D.size())
      FG_CUDA(cudaMemcpy2DAsync(D.data(), static_cast<std::size_t>(nr) * sizeof(cplx), A.data() + c0 * A.rows() + r0,
                                static_cast<std::size_t>(A.rows()) * sizeof(cplx),
                                static_cast<std::size_t>(nr) * sizeof(cplx), static_cast<std::size_t>(nc),
                                cudaMemcpyDeviceToDevice, stream_));
    return D;
  }
  void setBlock(DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, const DeviceMatrix& B) override {
    FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + B.rows() <= A.rows() && c0 + B.cols() <= A.cols(), "setBlock: out of range");
    if (B.size())
      FG_CUDA(cudaMemcpy2DAsync(A.data() + c0 * A.rows() + r0, static_cast<std::size_t>(A.rows()) * sizeof(cplx),
                                B.data(), static_cast<std::size_t>(B.rows()) * sizeof(cplx),
                                static_cast<std::size_t>(B.rows()) * sizeof(cplx), static_cast<std::size_t>(B.cols()),
                                cudaMemcpyDeviceToDevice, stream_));
  }
  DeviceMatrix columns(const DeviceMatrix& A, const std::vector<Eigen::Index>& idx) override {
    DeviceMatrix D = alloc(A.rows(), static_cast<Eigen::Index>(idx.size()));
    const std::size_t col = static_cast<std::size_t>(A.rows()) * sizeof(cplx);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      FG_REQUIRE(idx[k] >= 0 && idx[k] < A.cols(), "columns: index out of range");
      FG_CUDA(cudaMemcpyAsync(D.data() + static_cast<Eigen::Index>(k) * A.rows(), A.data() + idx[k] * A.rows(), col,
                              cudaMemcpyDeviceToDevice, stream_));
    }
    return D;
  }
  DeviceEig eigResident(const DeviceMatrix& A) override {
    FG_REQUIRE(A.rows() == A.cols(), "eig: matrix must be square");
    const std::int64_t n = A.rows();
    DeviceEig r;
    if (n == 0) return r;
    DeviceMatrix W = alloc(n, n), VR = alloc(n, n), w = alloc(n, 1);
    FG_CUDA(cudaMemcpyAsync(W.data(), A.data(), bytes(A), cudaMemcpyDeviceToDevice, stream_));  // Xgeev overwrites A
    info_.resize(1);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgeev_bufferSize(solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR, n,
                                           CUDA_C_64F, W.data(), n, CUDA_C_64F, w.data(), CUDA_C_64F, nullptr, n,
                                           CUDA_C_64F, VR.data(), n, CUDA_C_64F, &dBytes, &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    FG_CUSOLVER(cusolverDnXgeev(solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR, n, CUDA_C_64F,
                                W.data(), n, CUDA_C_64F, w.data(), CUDA_C_64F, nullptr, n, CUDA_C_64F, VR.data(), n,
                                CUDA_C_64F, work_.get(), dBytes, pinned_.get(hBytes > 0 ? hBytes : 1), hBytes,
                                info_.get()));
    checkInfo("cusolverDnXgeev");
    r.values = download(w);
    // unit-norm columns as Eigen: norms on the device, scaling with dgmm
    norms_.resize(static_cast<std::size_t>(n));
    FG_CUBLAS(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_DEVICE));
    for (std::int64_t k = 0; k < n; ++k)
      FG_CUBLAS(cublasDznrm2(blas_, static_cast<int>(n), z(VR) + k * n, 1, norms_.get() + k));
    FG_CUBLAS(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_HOST));
    std::vector<double> nrm(static_cast<std::size_t>(n));
    FG_CUDA(cudaMemcpyAsync(nrm.data(), norms_.get(), static_cast<std::size_t>(n) * sizeof(double),
                            cudaMemcpyDeviceToHost, stream_));
    synchronize();
    VecC inv(n);
    for (std::int64_t k = 0; k < n; ++k) inv(k) = nrm[static_cast<std::size_t>(k)] > 0 ? 1.0 / nrm[static_cast<std::size_t>(k)] : 1.0;
    r.vectors = scaleCols(VR, inv);
    return r;
  }

 private:
  static std::size_t bytes(const DeviceMatrix& A) { return static_cast<std::size_t>(A.size()) * sizeof(cplx); }
  static int i32(Eigen::Index v) { return static_cast<int>(v); }
  static cuDoubleComplex cz(cplx v) { return {v.real(), v.imag()}; }
  static cuDoubleComplex* z(const DeviceMatrix& A) { return reinterpret_cast<cuDoubleComplex*>(A.data()); }

  DeviceMatrix alloc(Eigen::Index rows, Eigen::Index cols) {
    const std::size_t n = static_cast<std::size_t>(std::max<Eigen::Index>(rows * cols, 1)) * sizeof(cplx);
    std::shared_ptr<MemoryPool> pool = pool_;
    void* p = pool->acquire(n);
    return DeviceMatrix(std::shared_ptr<void>(p, [pool, n](void* q) { pool->release(q, n); }), rows, cols);
  }
  DeviceMatrix dgmm(cublasSideMode_t side, const DeviceMatrix& A, const VecC& d) {
    FG_REQUIRE(d.size() == (side == CUBLAS_SIDE_LEFT ? A.rows() : A.cols()), "scale: size mismatch");
    DeviceMatrix D = alloc(A.rows(), A.cols());
    if (D.empty()) return D;
    const DeviceMatrix x = upload(d);
    FG_CUBLAS(cublasZdgmm(blas_, side, i32(A.rows()), i32(A.cols()), z(A), i32(A.rows()), z(x), 1, z(D), i32(A.rows())));
    return D;
  }

  void upload(cuDoubleComplex* d, const MatC& A) {
    FG_CUDA(cudaMemcpyAsync(d, A.data(), static_cast<std::size_t>(A.size()) * sizeof(cplx), cudaMemcpyHostToDevice,
                            stream_));
  }
  void download(MatC& A, const cuDoubleComplex* d) {
    FG_CUDA(cudaMemcpyAsync(A.data(), d, static_cast<std::size_t>(A.size()) * sizeof(cplx), cudaMemcpyDeviceToHost,
                            stream_));
    FG_CUDA(cudaStreamSynchronize(stream_));
  }
  void checkInfo(const char* what) {
    int info = 0;
    FG_CUDA(cudaMemcpyAsync(&info, info_.get(), sizeof(int), cudaMemcpyDeviceToHost, stream_));
    FG_CUDA(cudaStreamSynchronize(stream_));
    if (info != 0) throw Error(std::string(what) + ": info = " + std::to_string(info));
  }

  std::string name_;
  cudaStream_t stream_ = nullptr;
  cublasHandle_t blas_ = nullptr;
  cusolverDnHandle_t solver_ = nullptr;
  cusolverDnParams_t params_ = nullptr;
  DeviceBuffer<cuDoubleComplex> a_, b_, c_, w_;
  DeviceBuffer<std::int64_t> ipiv_;
  DeviceBuffer<int> info_;
  DeviceBuffer<char> work_;
  DeviceBuffer<double> norms_;
  PinnedBuffer pinned_;
  std::shared_ptr<MemoryPool> pool_ = std::make_shared<MemoryPool>(this);
};

}  // namespace

std::unique_ptr<Backend> makeCudaBackend(int device) { return std::make_unique<CudaBackend>(device); }

bool cudaAvailable() {
  int count = 0;
  return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

}  // namespace fg::linalg
