// FourierGrating - CUDA backend (milestone M10).
// Planned: cuBLAS ZGEMM, cuSOLVER Xgetrf/Xgetrs (64-bit API) and Xgeev (CUDA >= 12.6 Update 2,
// hybrid CPU/GPU, pinned host memory), cuFFT for field synthesis. See Plan section "CUDA-Backend".
#include <cuda_runtime.h>

#include <sstream>

#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"

namespace fg::linalg {

namespace {
void checkCuda(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw Error(std::string(what) + ": " + cudaGetErrorString(e));
}

class CudaBackend final : public Backend {
 public:
  explicit CudaBackend(int device) : device_(device) {
    int count = 0;
    checkCuda(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
    if (device >= count) throw Error("CUDA device index out of range");
    checkCuda(cudaSetDevice(device), "cudaSetDevice");
    cudaDeviceProp prop{};
    checkCuda(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
    std::ostringstream os;
    os << "cuda:" << prop.name << " (sm_" << prop.major << prop.minor << ", "
       << (prop.totalGlobalMem >> 20) << " MiB)";
    name_ = os.str();
  }
  std::string name() const override { return name_; }
  MatC gemm(const MatC&, const MatC&) override { throw NotImplemented("CudaBackend::gemm (M10)"); }
  MatC solve(const MatC&, const MatC&) override { throw NotImplemented("CudaBackend::solve (M10)"); }
  EigResult eig(const MatC&) override { throw NotImplemented("CudaBackend::eig (M10)"); }

 private:
  int device_;
  std::string name_;
};
}  // namespace

std::unique_ptr<Backend> makeCudaBackend(int device) { return std::make_unique<CudaBackend>(device); }

}  // namespace fg::linalg
