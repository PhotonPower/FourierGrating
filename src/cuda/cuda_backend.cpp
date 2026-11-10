// FourierGrating - CUDA backend (milestone M10, Plan section "CUDA-Backend").
// Host API of cuBLAS and cuSOLVER (64-bit Xgetrf/Xgetrs, Xgeev; CUDA >= 12.6 Update 2); no device
// kernels, so a plain C++ compiler suffices. Column-major data (Eigen default).
// M15: resident operations on DeviceMatrix handles (memory pool, one stream per backend); the host
// API (gemm / solve / eig on MatC) is built on them.
// M16: working precision. Device matrices hold complex<double> (Z) or complex<float> (C); every
// operation runs in the type of its operands (Cgemm, Xgetrf / Xgeev with CUDA_C_32F, no TF32).
// FP32 converts at upload / download (host, O(n^2)); Mixed: GEMM and eigenproblems in C, solves
// with cusolverDnZCgesv (C factorization, FP64 iterative refinement).
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <algorithm>
#include <complex>
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

#define FG_CUDA(call)                                                                       \
  do {                                                                                      \
    const cudaError_t e_ = (call);                                                          \
    if (e_ != cudaSuccess) throw Error(std::string(#call) + ": " + cudaGetErrorString(e_)); \
  } while (0)
#define FG_CUBLAS(call)                                                                            \
  do {                                                                                             \
    const cublasStatus_t e_ = (call);                                                              \
    if (e_ != CUBLAS_STATUS_SUCCESS)                                                               \
      throw Error(std::string(#call) + ": cuBLAS status " + std::to_string(static_cast<int>(e_))); \
  } while (0)
#define FG_CUSOLVER(call)                                     \
  do {                                                        \
    const cusolverStatus_t e_ = (call);                       \
    if (e_ != CUSOLVER_STATUS_SUCCESS)                        \
      throw Error(std::string(#call) + ": cuSOLVER status " + \
                  std::to_string(static_cast<int>(e_)));      \
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

/// RAII pinned host memory (host workspaces of Xgetrf / Xgeev).
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
static_assert(sizeof(std::complex<float>) == sizeof(cuComplex), "complex layout");

/// Device memory of the DeviceMatrix handles of one backend (M15, step 5): released buffers are
/// kept for reuse (no cudaFree, which would synchronize the device). A buffer released by another
/// thread is reused only after that thread's queued work has finished (see release()).
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

/// cuBLAS / CUDA types of the two working precisions.
template <class T>
struct Cu;
template <>
struct Cu<cuDoubleComplex> {
  using Host = cplx;
  static constexpr cudaDataType type = CUDA_C_64F;
  static cuDoubleComplex value(cplx v) { return {v.real(), v.imag()}; }
  static cublasStatus_t gemm(cublasHandle_t h, int m, int n, int k, const cuDoubleComplex* a,
                             const cuDoubleComplex* A, const cuDoubleComplex* B,
                             const cuDoubleComplex* b, cuDoubleComplex* C) {
    return cublasZgemm(h, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, a, A, m, B, k, b, C, m);
  }
  static cublasStatus_t geam(cublasHandle_t h, cublasOperation_t ta, cublasOperation_t tb, int m,
                             int n, const cuDoubleComplex* a, const cuDoubleComplex* A, int lda,
                             const cuDoubleComplex* b, const cuDoubleComplex* B, int ldb,
                             cuDoubleComplex* C) {
    return cublasZgeam(h, ta, tb, m, n, a, A, lda, b, B, ldb, C, m);
  }
  static cublasStatus_t dgmm(cublasHandle_t h, cublasSideMode_t s, int m, int n,
                             const cuDoubleComplex* A, const cuDoubleComplex* x,
                             cuDoubleComplex* C) {
    return cublasZdgmm(h, s, m, n, A, m, x, 1, C, m);
  }
  static cublasStatus_t nrm2(cublasHandle_t h, int n, const cuDoubleComplex* x, double* r) {
    return cublasDznrm2(h, n, x, 1, r);
  }
  using Real = double;
};
template <>
struct Cu<cuComplex> {
  using Host = std::complex<float>;
  static constexpr cudaDataType type = CUDA_C_32F;
  static cuComplex value(cplx v) {
    return {static_cast<float>(v.real()), static_cast<float>(v.imag())};
  }
  static cublasStatus_t gemm(cublasHandle_t h, int m, int n, int k, const cuComplex* a,
                             const cuComplex* A, const cuComplex* B, const cuComplex* b,
                             cuComplex* C) {
    return cublasCgemm(h, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, a, A, m, B, k, b, C, m);
  }
  static cublasStatus_t geam(cublasHandle_t h, cublasOperation_t ta, cublasOperation_t tb, int m,
                             int n, const cuComplex* a, const cuComplex* A, int lda,
                             const cuComplex* b, const cuComplex* B, int ldb, cuComplex* C) {
    return cublasCgeam(h, ta, tb, m, n, a, A, lda, b, B, ldb, C, m);
  }
  static cublasStatus_t dgmm(cublasHandle_t h, cublasSideMode_t s, int m, int n, const cuComplex* A,
                             const cuComplex* x, cuComplex* C) {
    return cublasCdgmm(h, s, m, n, A, m, x, 1, C, m);
  }
  static cublasStatus_t nrm2(cublasHandle_t h, int n, const cuComplex* x, float* r) {
    return cublasScnrm2(h, n, x, 1, r);
  }
  using Real = float;
};

constexpr std::size_t kZ = sizeof(cuDoubleComplex), kC = sizeof(cuComplex);

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
    os << "cuda:" << prop.name << " (sm_" << prop.major << prop.minor << ", "
       << (prop.totalGlobalMem >> 20) << " MiB)";
    name_ = os.str();
    FG_CUDA(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
    FG_CUBLAS(cublasCreate(&blas_));
    FG_CUBLAS(cublasSetStream(blas_, stream_));
    FG_CUBLAS(cublasSetMathMode(blas_, CUBLAS_DEFAULT_MATH));  // no TF32 (Plan section "M16")
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

  // ---- host API: upload, resident operation, download (precision of the backend, M16) ----
  MatC gemm(const MatC& A, const MatC& B) override {
    FG_REQUIRE(A.cols() == B.rows(), "gemm: size mismatch");
    const std::size_t e = precision_ == Precision::FP64 ? kZ : kC;  // GEMM: C in Mixed and FP32
    return download(multiply(uploadAs(A, e), uploadAs(B, e), 1.0, nullptr, 0.0));
  }
  MatC solve(const MatC& A, const MatC& B) override {
    FG_REQUIRE(A.rows() == A.cols() && A.rows() == B.rows(), "solve: size mismatch");
    const std::size_t e =
        precision_ == Precision::FP32 ? kC : kZ;  // Mixed: Z data, C factorization
    return download(luSolve(uploadAs(A, e), uploadAs(B, e)));
  }
  EigResult eig(const MatC& A) override {
    FG_REQUIRE(A.rows() == A.cols(), "eig: matrix must be square");
    const std::size_t e = precision_ == Precision::FP64 ? kZ : kC;
    DeviceEig r = eigResident(uploadAs(A, e));
    return {std::move(r.values), download(r.vectors)};
  }

  // ---- resident operations (M15, step 5): everything on this backend's stream ----
  /// Mixed keeps Z storage for the refined solves; its resident chain would need device-side type
  /// conversion (kernels), so Mixed uses the host API.
  bool resident() const override { return precision_ != Precision::Mixed; }
  void synchronize() override { FG_CUDA(cudaStreamSynchronize(stream_)); }

  DeviceMatrix upload(const MatC& A) override {
    return uploadAs(A, precision_ == Precision::FP32 ? kC : kZ);
  }
  MatC download(const DeviceMatrix& A) override {
    MatC M(A.rows(), A.cols());
    if (A.empty()) return M;
    if (A.elementSize() == kZ) {
      FG_CUDA(cudaMemcpyAsync(M.data(), A.raw(), bytes(A), cudaMemcpyDeviceToHost, stream_));
      synchronize();
      return M;
    }
    Eigen::MatrixXcf F(A.rows(), A.cols());
    FG_CUDA(cudaMemcpyAsync(F.data(), A.raw(), bytes(A), cudaMemcpyDeviceToHost, stream_));
    synchronize();
    return F.cast<cplx>();
  }
  DeviceMatrix zeros(Eigen::Index rows, Eigen::Index cols) override {
    return zerosAs(rows, cols, elem());
  }
  DeviceMatrix identity(Eigen::Index n) override { return identityAs(n, elem()); }
  DeviceMatrix multiply(const DeviceMatrix& A, const DeviceMatrix& B, cplx alpha,
                        const DeviceMatrix* C, cplx beta) override {
    FG_REQUIRE(A.cols() == B.rows() && (!C || (C->rows() == A.rows() && C->cols() == B.cols())),
               "multiply: size mismatch");
    same(A, B);
    if (C) same(A, *C);
    return A.elementSize() == kZ ? multiplyT<cuDoubleComplex>(A, B, alpha, C, beta)
                                 : multiplyT<cuComplex>(A, B, alpha, C, beta);
  }
  DeviceMatrix luSolve(const DeviceMatrix& A, const DeviceMatrix& B) override {
    FG_REQUIRE(A.rows() == A.cols() && A.rows() == B.rows(), "luSolve: size mismatch");
    same(A, B);
    if (A.elementSize() == kC) return luSolveT<cuComplex>(A, B);
    if (precision_ == Precision::Mixed) return refinedSolve(A, B);
    return luSolveT<cuDoubleComplex>(A, B);
  }
  DeviceMatrix add(cplx alpha, const DeviceMatrix& A, bool transA, cplx beta, const DeviceMatrix& B,
                   bool transB) override {
    if (!B.empty()) same(A, B);
    return A.elementSize() == kZ ? addT<cuDoubleComplex>(alpha, A, transA, beta, B, transB)
                                 : addT<cuComplex>(alpha, A, transA, beta, B, transB);
  }
  DeviceMatrix scaleRows(const VecC& d, const DeviceMatrix& A) override {
    return dgmm(CUBLAS_SIDE_LEFT, A, d);
  }
  DeviceMatrix scaleCols(const DeviceMatrix& A, const VecC& d) override {
    return dgmm(CUBLAS_SIDE_RIGHT, A, d);
  }
  DeviceMatrix block(const DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, Eigen::Index nr,
                     Eigen::Index nc) override {
    FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + nr <= A.rows() && c0 + nc <= A.cols(),
               "block: out of range");
    const std::size_t e = A.elementSize();
    DeviceMatrix D = alloc(nr, nc, e);
    if (D.size())
      FG_CUDA(cudaMemcpy2DAsync(D.raw(), static_cast<std::size_t>(nr) * e, at(A, r0, c0),
                                static_cast<std::size_t>(A.rows()) * e,
                                static_cast<std::size_t>(nr) * e, static_cast<std::size_t>(nc),
                                cudaMemcpyDeviceToDevice, stream_));
    return D;
  }
  void setBlock(DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, const DeviceMatrix& B) override {
    FG_REQUIRE(r0 >= 0 && c0 >= 0 && r0 + B.rows() <= A.rows() && c0 + B.cols() <= A.cols(),
               "setBlock: out of range");
    same(A, B);
    const std::size_t e = A.elementSize();
    if (B.size())
      FG_CUDA(cudaMemcpy2DAsync(
          at(A, r0, c0), static_cast<std::size_t>(A.rows()) * e, B.raw(),
          static_cast<std::size_t>(B.rows()) * e, static_cast<std::size_t>(B.rows()) * e,
          static_cast<std::size_t>(B.cols()), cudaMemcpyDeviceToDevice, stream_));
  }
  DeviceMatrix columns(const DeviceMatrix& A, const std::vector<Eigen::Index>& idx) override {
    const std::size_t e = A.elementSize();
    DeviceMatrix D = alloc(A.rows(), static_cast<Eigen::Index>(idx.size()), e);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      FG_REQUIRE(idx[k] >= 0 && idx[k] < A.cols(), "columns: index out of range");
      FG_CUDA(cudaMemcpyAsync(at(D, 0, static_cast<Eigen::Index>(k)), at(A, 0, idx[k]),
                              static_cast<std::size_t>(A.rows()) * e, cudaMemcpyDeviceToDevice,
                              stream_));
    }
    return D;
  }
  DeviceEig eigResident(const DeviceMatrix& A) override {
    FG_REQUIRE(A.rows() == A.cols(), "eig: matrix must be square");
    return A.elementSize() == kZ ? eigT<cuDoubleComplex>(A) : eigT<cuComplex>(A);
  }

 private:
  std::size_t elem() const { return precision_ == Precision::FP32 ? kC : kZ; }
  static std::size_t bytes(const DeviceMatrix& A) {
    return static_cast<std::size_t>(A.size()) * A.elementSize();
  }
  static int i32(Eigen::Index v) { return static_cast<int>(v); }
  static void same(const DeviceMatrix& A, const DeviceMatrix& B) {
    FG_REQUIRE(A.elementSize() == B.elementSize(), "device matrices of different precision");
  }
  static void* at(const DeviceMatrix& A, Eigen::Index r, Eigen::Index c) {
    return static_cast<char*>(A.raw()) +
           static_cast<std::size_t>(c * A.rows() + r) * A.elementSize();
  }
  template <class T>
  static T* ptr(const DeviceMatrix& A) {
    return static_cast<T*>(A.raw());
  }

  DeviceMatrix alloc(Eigen::Index rows, Eigen::Index cols, std::size_t e) {
    const std::size_t n = static_cast<std::size_t>(std::max<Eigen::Index>(rows * cols, 1)) * e;
    std::shared_ptr<MemoryPool> pool = pool_;
    void* p = pool->acquire(n);
    return DeviceMatrix(std::shared_ptr<void>(p, [pool, n](void* q) { pool->release(q, n); }), rows,
                        cols, e);
  }
  DeviceMatrix uploadAs(const MatC& A, std::size_t e) {
    DeviceMatrix D = alloc(A.rows(), A.cols(), e);
    if (A.size() == 0) return D;
    if (e == kZ) {
      FG_CUDA(cudaMemcpyAsync(D.raw(), A.data(), bytes(D), cudaMemcpyHostToDevice, stream_));
    } else {
      const Eigen::MatrixXcf F =
          A.cast<std::complex<float>>();  // conversion at the boundary, O(n^2)
      FG_CUDA(cudaMemcpyAsync(D.raw(), F.data(), bytes(D), cudaMemcpyHostToDevice, stream_));
      synchronize();  // F is local
      return D;
    }
    synchronize();  // A may be a temporary of the caller
    return D;
  }
  DeviceMatrix uploadVector(const VecC& d, std::size_t e) {
    return uploadAs(MatC(Eigen::Map<const MatC>(d.data(), d.size(), 1)), e);
  }
  DeviceMatrix zerosAs(Eigen::Index rows, Eigen::Index cols, std::size_t e) {
    DeviceMatrix D = alloc(rows, cols, e);
    if (D.size()) FG_CUDA(cudaMemsetAsync(D.raw(), 0, bytes(D), stream_));
    return D;
  }
  DeviceMatrix identityAs(Eigen::Index n, std::size_t e) {
    DeviceMatrix D = zerosAs(n, n, e);
    if (n == 0) return D;
    std::vector<char> ones(static_cast<std::size_t>(n) * e);
    for (Eigen::Index k = 0; k < n; ++k) {
      char* p = ones.data() + static_cast<std::size_t>(k) * e;
      if (e == kZ)
        *reinterpret_cast<cplx*>(p) = 1.0;
      else
        *reinterpret_cast<std::complex<float>*>(p) = 1.0f;
    }
    // the diagonal: n elements with pitch (n + 1) elements
    FG_CUDA(cudaMemcpy2DAsync(D.raw(), static_cast<std::size_t>(n + 1) * e, ones.data(), e, e,
                              static_cast<std::size_t>(n), cudaMemcpyHostToDevice, stream_));
    synchronize();  // ones is a local host buffer
    return D;
  }

  template <class T>
  DeviceMatrix multiplyT(const DeviceMatrix& A, const DeviceMatrix& B, cplx alpha,
                         const DeviceMatrix* C, cplx beta) {
    DeviceMatrix D = alloc(A.rows(), B.cols(), sizeof(T));
    if (D.empty()) return D;
    T b = Cu<T>::value(0.0);
    if (C) {
      FG_CUDA(cudaMemcpyAsync(D.raw(), C->raw(), bytes(D), cudaMemcpyDeviceToDevice, stream_));
      b = Cu<T>::value(beta);
    }
    if (A.cols() == 0) {
      if (!C) FG_CUDA(cudaMemsetAsync(D.raw(), 0, bytes(D), stream_));
      return D;
    }
    const T a = Cu<T>::value(alpha);
    FG_CUBLAS(Cu<T>::gemm(blas_, i32(A.rows()), i32(B.cols()), i32(A.cols()), &a, ptr<T>(A),
                          ptr<T>(B), &b, ptr<T>(D)));
    return D;
  }

  template <class T>
  DeviceMatrix luSolveT(const DeviceMatrix& A, const DeviceMatrix& B) {
    const std::int64_t n = A.rows(), nrhs = B.cols();
    DeviceMatrix X = alloc(n, nrhs, sizeof(T));
    if (n == 0 || nrhs == 0) return X;
    DeviceMatrix LU = alloc(n, n, sizeof(T));
    FG_CUDA(cudaMemcpyAsync(LU.raw(), A.raw(), bytes(A), cudaMemcpyDeviceToDevice, stream_));
    FG_CUDA(cudaMemcpyAsync(X.raw(), B.raw(), bytes(B), cudaMemcpyDeviceToDevice, stream_));
    ipiv_.resize(static_cast<std::size_t>(n));
    info_.resize(1);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgetrf_bufferSize(solver_, params_, n, n, Cu<T>::type, LU.raw(), n,
                                            Cu<T>::type, &dBytes, &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    void* hwork = pinned_.get(hBytes > 0 ? hBytes : 1);
    FG_CUSOLVER(cusolverDnXgetrf(solver_, params_, n, n, Cu<T>::type, LU.raw(), n, ipiv_.get(),
                                 Cu<T>::type, work_.get(), dBytes, hwork, hBytes, info_.get()));
    checkInfo("cusolverDnXgetrf");
    FG_CUSOLVER(cusolverDnXgetrs(solver_, params_, CUBLAS_OP_N, n, nrhs, Cu<T>::type, LU.raw(), n,
                                 ipiv_.get(), Cu<T>::type, X.raw(), n, info_.get()));
    checkInfo("cusolverDnXgetrs");
    return X;
  }

  /// Mixed (Plan section "M16"): complex-single factorization with FP64 iterative refinement
  /// (cusolverDnZCgesv; it falls back to an FP64 factorization itself if the refinement fails).
  DeviceMatrix refinedSolve(const DeviceMatrix& A, const DeviceMatrix& B) {
    const int n = i32(A.rows()), nrhs = i32(B.cols());
    DeviceMatrix X = alloc(n, nrhs, kZ);
    if (n == 0 || nrhs == 0) return X;
    DeviceMatrix Ac = alloc(n, n, kZ), Bc = alloc(n, nrhs, kZ);
    FG_CUDA(cudaMemcpyAsync(Ac.raw(), A.raw(), bytes(A), cudaMemcpyDeviceToDevice, stream_));
    FG_CUDA(cudaMemcpyAsync(Bc.raw(), B.raw(), bytes(B), cudaMemcpyDeviceToDevice, stream_));
    ipiv32_.resize(static_cast<std::size_t>(n));
    info_.resize(1);
    std::size_t wBytes = 0;
    FG_CUSOLVER(cusolverDnZCgesv_bufferSize(solver_, n, nrhs, ptr<cuDoubleComplex>(Ac), n,
                                            ipiv32_.get(), ptr<cuDoubleComplex>(Bc), n,
                                            ptr<cuDoubleComplex>(X), n, nullptr, &wBytes));
    work_.resize(wBytes > 0 ? wBytes : 1);
    int iter = 0;
    FG_CUSOLVER(cusolverDnZCgesv(solver_, n, nrhs, ptr<cuDoubleComplex>(Ac), n, ipiv32_.get(),
                                 ptr<cuDoubleComplex>(Bc), n, ptr<cuDoubleComplex>(X), n,
                                 work_.get(), wBytes, &iter, info_.get()));
    checkInfo("cusolverDnZCgesv");
    if (RefinementCounter* c = refinementCounter()) {
      if (iter >= 0)
        c->steps += iter;
      else
        ++c->fallbacks;
    }
    return X;
  }

  template <class T>
  DeviceMatrix addT(cplx alpha, const DeviceMatrix& A, bool transA, cplx beta,
                    const DeviceMatrix& B, bool transB) {
    const Eigen::Index m = transA ? A.cols() : A.rows(), n = transA ? A.rows() : A.cols();
    DeviceMatrix D = alloc(m, n, sizeof(T));
    if (D.empty()) return D;
    const bool useB = !B.empty();
    if (useB)
      FG_REQUIRE((transB ? B.cols() : B.rows()) == m && (transB ? B.rows() : B.cols()) == n,
                 "add: size mismatch");
    const T a = Cu<T>::value(alpha), b = Cu<T>::value(useB ? beta : 0.0);
    const DeviceMatrix& Bm = useB ? B : A;  // not read for beta = 0
    const bool tb = useB ? transB : transA;
    FG_CUBLAS(Cu<T>::geam(blas_, transA ? CUBLAS_OP_T : CUBLAS_OP_N, tb ? CUBLAS_OP_T : CUBLAS_OP_N,
                          i32(m), i32(n), &a, ptr<T>(A), i32(A.rows()), &b, ptr<T>(Bm),
                          i32(Bm.rows()), ptr<T>(D)));
    return D;
  }

  DeviceMatrix dgmm(cublasSideMode_t side, const DeviceMatrix& A, const VecC& d) {
    FG_REQUIRE(d.size() == (side == CUBLAS_SIDE_LEFT ? A.rows() : A.cols()),
               "scale: size mismatch");
    const std::size_t e = A.elementSize();
    DeviceMatrix D = alloc(A.rows(), A.cols(), e);
    if (D.empty()) return D;
    const DeviceMatrix x = uploadVector(d, e);
    if (e == kZ)
      FG_CUBLAS(Cu<cuDoubleComplex>::dgmm(blas_, side, i32(A.rows()), i32(A.cols()),
                                          ptr<cuDoubleComplex>(A), ptr<cuDoubleComplex>(x),
                                          ptr<cuDoubleComplex>(D)));
    else
      FG_CUBLAS(Cu<cuComplex>::dgmm(blas_, side, i32(A.rows()), i32(A.cols()), ptr<cuComplex>(A),
                                    ptr<cuComplex>(x), ptr<cuComplex>(D)));
    return D;
  }

  template <class T>
  DeviceEig eigT(const DeviceMatrix& A) {
    using Real = typename Cu<T>::Real;
    const std::int64_t n = A.rows();
    DeviceEig r;
    if (n == 0) return r;
    DeviceMatrix W = alloc(n, n, sizeof(T)), VR = alloc(n, n, sizeof(T)),
                 w = alloc(n, 1, sizeof(T));
    FG_CUDA(cudaMemcpyAsync(W.raw(), A.raw(), bytes(A), cudaMemcpyDeviceToDevice,
                            stream_));  // Xgeev overwrites A
    info_.resize(1);
    std::size_t dBytes = 0, hBytes = 0;
    FG_CUSOLVER(cusolverDnXgeev_bufferSize(
        solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR, n, Cu<T>::type,
        W.raw(), n, Cu<T>::type, w.raw(), Cu<T>::type, nullptr, n, Cu<T>::type, VR.raw(), n,
        Cu<T>::type, &dBytes, &hBytes));
    work_.resize(dBytes > 0 ? dBytes : 1);
    {
      // cusolverDnXgeev fails with CUSOLVER_STATUS_INTERNAL_ERROR (7) when two threads call it at
      // once, although each has its own handle, stream and workspaces (CUDA 12.9, RTX 6000 Ada;
      // reproduced with FG_NUM_THREADS >= 2 and FG_GPU_STREAMS >= 2). The lock covers the call
      // and checkInfo (which waits for the stream); the rest of the pipeline stays concurrent.
      static std::mutex xgeevMu;
      const std::lock_guard<std::mutex> lock(xgeevMu);
      FG_CUSOLVER(cusolverDnXgeev(solver_, params_, CUSOLVER_EIG_MODE_NOVECTOR,
                                  CUSOLVER_EIG_MODE_VECTOR, n, Cu<T>::type, W.raw(), n,
                                  Cu<T>::type, w.raw(), Cu<T>::type, nullptr, n, Cu<T>::type,
                                  VR.raw(), n, Cu<T>::type, work_.get(), dBytes,
                                  pinned_.get(hBytes > 0 ? hBytes : 1), hBytes, info_.get()));
      checkInfo("cusolverDnXgeev");
    }
    const MatC vals = download(w);
    r.values = vals.col(0);
    // unit-norm columns as Eigen: norms on the device, scaling with dgmm
    norms_.resize(static_cast<std::size_t>(n) * sizeof(Real));
    Real* nd = reinterpret_cast<Real*>(norms_.get());
    FG_CUBLAS(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_DEVICE));
    for (std::int64_t k = 0; k < n; ++k)
      FG_CUBLAS(Cu<T>::nrm2(blas_, static_cast<int>(n), ptr<T>(VR) + k * n, nd + k));
    FG_CUBLAS(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_HOST));
    std::vector<Real> nrm(static_cast<std::size_t>(n));
    FG_CUDA(cudaMemcpyAsync(nrm.data(), nd, static_cast<std::size_t>(n) * sizeof(Real),
                            cudaMemcpyDeviceToHost, stream_));
    synchronize();
    VecC inv(n);
    for (std::int64_t k = 0; k < n; ++k) {
      const double v = static_cast<double>(nrm[static_cast<std::size_t>(k)]);
      inv(k) = v > 0 ? 1.0 / v : 1.0;
    }
    r.vectors = scaleCols(VR, inv);
    return r;
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
  DeviceBuffer<std::int64_t> ipiv_;
  DeviceBuffer<int> ipiv32_;
  DeviceBuffer<int> info_;
  DeviceBuffer<char> work_;
  DeviceBuffer<char> norms_;
  PinnedBuffer pinned_;
  std::shared_ptr<MemoryPool> pool_ = std::make_shared<MemoryPool>(this);
};

}  // namespace

std::unique_ptr<Backend> makeCudaBackend(int device) {
  return std::make_unique<CudaBackend>(device);
}

bool cudaAvailable() {
  int count = 0;
  return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

}  // namespace fg::linalg
