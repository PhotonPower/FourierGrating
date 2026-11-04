// FourierGrating - CUDA backend (milestone M10, Plan section "CUDA-Backend").
// Host API of cuBLAS (Zgemm) and cuSOLVER (64-bit Xgetrf/Xgetrs, Xgeev; CUDA >= 12.6 Update 2);
// no device kernels, so a plain C++ compiler suffices. Column-major data (Eigen default).
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <cstdint>
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

 private:
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
  PinnedBuffer pinned_;
};

}  // namespace

std::unique_ptr<Backend> makeCudaBackend(int device) { return std::make_unique<CudaBackend>(device); }

bool cudaAvailable() {
  int count = 0;
  return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

}  // namespace fg::linalg
