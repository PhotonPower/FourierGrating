// FourierGrating - linear algebra backend abstraction (Plan section "CUDA-Backend").
// Physics code calls only this interface; CPU (Eigen/LAPACK) is the reference implementation.
#pragma once

#include <memory>
#include <string>

#include "fg/core/types.hpp"

namespace fg::linalg {

struct EigResult {
  VecC values;
  MatC vectors;
};

class Backend {
 public:
  virtual ~Backend() = default;
  virtual std::string name() const = 0;
  virtual MatC gemm(const MatC& A, const MatC& B) = 0;
  virtual MatC solve(const MatC& A, const MatC& B) = 0;   ///< LU with partial pivoting
  virtual EigResult eig(const MatC& A) = 0;               ///< general complex eigenproblem
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
bool offload(Eigen::Index n);
/// Thread-local device backend (created on first use).
Backend& device();

}  // namespace fg::linalg
