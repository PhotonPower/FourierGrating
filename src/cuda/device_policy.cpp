// FourierGrating - default working precision per GPU (Plan section "M16: Standardwahl", decision of
// the project owner): "precision": "auto" means FP32 with FP64 fallback on a GeForce RTX 3090 (user
// installation, FP64 at 1/64 of the FP32 rate), FP64 elsewhere (server, RTX 6000 Ada). Add further
// cards to the table.
#include <cuda_runtime.h>

#include <cstring>

#include "fg/linalg/backend.hpp"

namespace fg::linalg {

namespace {

struct Rule {
  const char* name;  ///< substring of cudaDeviceProp::name
  int major, minor;  ///< compute capability
  Precision precision;
};

constexpr Rule kRules[] = {
    {"GeForce RTX 3090", 8, 6, Precision::FP32},
};

Precision lookup() {
  int dev = 0;
  cudaDeviceProp p{};
  if (cudaGetDevice(&dev) != cudaSuccess || cudaGetDeviceProperties(&p, dev) != cudaSuccess)
    return Precision::FP64;
  for (const Rule& r : kRules)
    if (std::strstr(p.name, r.name) && p.major == r.major && p.minor == r.minor) return r.precision;
  return Precision::FP64;
}

}  // namespace

Precision gpuAutoPrecision() {
  static const Precision p = lookup();
  return p;
}

}  // namespace fg::linalg
