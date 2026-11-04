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
/// CUDA backend (cuBLAS/cuSOLVER); throws if built without FG_ENABLE_CUDA. TODO(M10).
std::unique_ptr<Backend> makeCudaBackend(int device = 0);
/// Picks CUDA above a size threshold if available, CPU otherwise.
std::unique_ptr<Backend> makeBackend(const std::string& preference, Eigen::Index problemSize);

}  // namespace fg::linalg
