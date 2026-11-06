// FourierGrating - linear algebra backend abstraction (Plan section "CUDA-Backend").
// Physics code calls only this interface; CPU (Eigen/LAPACK) is the reference implementation.
#pragma once

#include <atomic>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fg/core/types.hpp"

namespace fg::linalg {

struct EigResult {
  VecC values;
  MatC vectors;
};

/// Working precision of the O(n^3) operations (Plan section "M16: Waehlbare Genauigkeit"): the
/// physics code always works with MatC (FP64); the backend converts at its boundary (O(n^2)).
///   FP64   Z everywhere (default, bit-identical)
///   Mixed  GEMM and eigenproblems in C, LU in C with FP64 iterative refinement
///   FP32   GEMM, LU and eigenproblems in C
enum class Precision { FP64, Mixed, FP32 };
std::string precisionName(Precision p);
/// "fp64" | "mixed" | "fp32"
Precision parsePrecision(const std::string& name);
/// Unit roundoff of the working precision (2^-53 or 2^-24).
double unitRoundoff(Precision p);

/// Complex column-major matrix in the memory of a backend (Plan section "M15", step 5: GPU-resident
/// chains). RAII handle: the memory goes back to the pool of the creating backend when the last
/// handle is dropped (after the dropping thread's queued work has finished).
class DeviceMatrix {
 public:
  DeviceMatrix() = default;
  /// elem: bytes per element (complex<double>; complex<float> on a GPU backend in FP32, M16)
  DeviceMatrix(std::shared_ptr<void> mem, Eigen::Index rows, Eigen::Index cols,
               std::size_t elem = sizeof(cplx))
      : mem_(std::move(mem)), rows_(rows), cols_(cols), elem_(elem) {}
  Eigen::Index rows() const { return rows_; }
  Eigen::Index cols() const { return cols_; }
  Eigen::Index size() const { return rows_ * cols_; }
  bool empty() const { return size() == 0; }
  std::size_t elementSize() const { return elem_; }
  void* raw() const { return mem_.get(); }
  /// complex<double> storage only (host reference)
  cplx* data() const { return static_cast<cplx*>(mem_.get()); }

 private:
  std::shared_ptr<void> mem_;
  Eigen::Index rows_ = 0, cols_ = 0;
  std::size_t elem_ = sizeof(cplx);
};

/// Eigenvalues (host) and unit-norm right eigenvectors (device) of a general complex matrix.
struct DeviceEig {
  VecC values;
  DeviceMatrix vectors;
};

class Backend {
 public:
  virtual ~Backend() = default;
  virtual std::string name() const = 0;
  virtual MatC gemm(const MatC& A, const MatC& B) = 0;
  virtual MatC solve(const MatC& A, const MatC& B) = 0;  ///< LU with partial pivoting
  virtual EigResult eig(const MatC& A) = 0;              ///< general complex eigenproblem
  /// Working precision of gemm / solve / eig and of the resident operations (M16).
  virtual void setPrecision(Precision p) { precision_ = p; }
  Precision precision() const { return precision_; }

  // ---- resident operations on DeviceMatrix handles (M15, step 5) ----
  // The default implementations keep the matrices in host memory (Eigen): reference for the tests
  // of the resident physics path. The CUDA backend keeps them on the GPU (resident() = true); only
  // then does the solver use the resident path.
  virtual bool resident() const { return false; }
  virtual DeviceMatrix upload(const MatC& A);
  virtual MatC download(const DeviceMatrix& A);
  virtual DeviceMatrix zeros(Eigen::Index rows, Eigen::Index cols);
  virtual DeviceMatrix identity(Eigen::Index n);
  /// alpha A B + beta C (C may be null: beta ignored)
  virtual DeviceMatrix multiply(const DeviceMatrix& A, const DeviceMatrix& B, cplx alpha = 1.0,
                                const DeviceMatrix* C = nullptr, cplx beta = 0.0);
  /// A^-1 B (LU with partial pivoting)
  virtual DeviceMatrix luSolve(const DeviceMatrix& A, const DeviceMatrix& B);
  /// alpha op(A) + beta op(B), op = transpose if the flag is set; B may be empty (alpha op(A))
  virtual DeviceMatrix add(cplx alpha, const DeviceMatrix& A, bool transA, cplx beta,
                           const DeviceMatrix& B, bool transB);
  virtual DeviceMatrix scaleRows(const VecC& d, const DeviceMatrix& A);  ///< diag(d) A
  virtual DeviceMatrix scaleCols(const DeviceMatrix& A, const VecC& d);  ///< A diag(d)
  virtual DeviceMatrix block(const DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0,
                             Eigen::Index nr, Eigen::Index nc);
  virtual void setBlock(DeviceMatrix& A, Eigen::Index r0, Eigen::Index c0, const DeviceMatrix& B);
  /// Columns idx of A (in this order).
  virtual DeviceMatrix columns(const DeviceMatrix& A, const std::vector<Eigen::Index>& idx);
  virtual DeviceEig eigResident(const DeviceMatrix& A);
  /// Waits for all queued work of this backend (results visible to other threads).
  virtual void synchronize() {}

 protected:
  /// Host storage for the default implementations.
  static DeviceMatrix hostMatrix(Eigen::Index rows, Eigen::Index cols);
  Precision precision_ = Precision::FP64;
};

/// CPU backend (Eigen). Always available.
std::unique_ptr<Backend> makeCpuBackend();
/// CUDA backend (cuBLAS/cuSOLVER host API: Zgemm, Xgetrf/Xgetrs, Xgeev); throws if built without
/// FG_ENABLE_CUDA or if no device is available.
std::unique_ptr<Backend> makeCudaBackend(int device = 0);
/// Picks CUDA above a size threshold if available, CPU otherwise.
std::unique_ptr<Backend> makeBackend(const std::string& preference, Eigen::Index problemSize);
/// True if this build contains the CUDA backend and a device is usable.
bool cudaAvailable();

// ---- dispatch for the physics code (Plan section "CUDA-Backend", M10) ----
// The O(n^3) operations (RCWA eigenproblems, stable layer S-matrix, star products, T -> S, RK step
// propagators) ask offload(n); if false they run their original Eigen code (CPU reference, bitwise
// unchanged), otherwise they use device(). Settings are thread-local (one solver per thread).
/// preference: "cpu" | "cuda" | "auto" (CUDA from matrix size `threshold` on, if available).
/// The environment variables FG_BACKEND and FG_BACKEND_THRESHOLD override both arguments
/// (e.g. the whole validation suite on the GPU: FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest).
void setPreference(const std::string& preference, Eigen::Index threshold = 1000);
/// Working precision of this thread (M16; FP64 by default). With FP32/Mixed every O(n^3) operation
/// goes through a backend (the GPU above the threshold, else the CPU backend in that precision).
void setPrecision(Precision p);
Precision precision();
/// True if an n x n operation runs on the GPU (preference, threshold, availability).
bool gpuSelected(Eigen::Index n);
/// True if the n x n operation goes through device() (GPU, or a CPU backend in FP32/Mixed); false:
/// the original Eigen code in FP64.
bool offload(Eigen::Index n);
/// The backend chosen by the last offload(n) of this thread.
Backend& device();
/// A B through device() if offload(A.rows()), else Eigen in FP64 (unchanged CPU path).
MatC product(const MatC& A, const MatC& B);
/// Thread-local GPU backend (created on first use; resident chain).
Backend& gpu();
/// Counts the iterative-refinement steps and FP64 fallbacks of mixed-precision solves (M16). Each
/// thread reports to the counter set with setRefinementCounter (none: not counted).
struct RefinementCounter {
  std::atomic<long> steps{0}, fallbacks{0};
};
void setRefinementCounter(RefinementCounter* c);
RefinementCounter* refinementCounter();
/// Precision of "precision": "auto" on this GPU (src/cuda/device_policy.cpp): "fp32" on a GeForce
/// RTX 3090, else "fp64"; "fp64" without CUDA.
Precision gpuAutoPrecision();
/// Threads (streams) for independent slices on the device (FG_GPU_STREAMS, default 4; M15).
int deviceStreams();

namespace detail {
/// The calling thread's device backend if it has one (else nullptr): a DeviceMatrix released by
/// another thread first waits for that thread's queued work (resident chain, M15).
Backend* currentDevice();
}  // namespace detail

}  // namespace fg::linalg
