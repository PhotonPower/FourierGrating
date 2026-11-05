#include "fg/core/parallel.hpp"

#include <cstdlib>

namespace fg {

int& detail::parallelDepth() {
  thread_local int depth = 0;
  return depth;
}

int availableThreads() {
  if (const char* env = std::getenv("FG_NUM_THREADS"); env && *env) {
    const int v = std::atoi(env);
    if (v > 0) return v;
  }
#ifdef FG_HAVE_OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
#endif
}

}  // namespace fg
