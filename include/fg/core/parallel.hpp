// FourierGrating - thread model (Plan section "M15: Profiling und Parallelisierung", step 2).
// parallelFor runs independent, index-addressed tasks; results are stored by index, so they do not
// depend on the thread count. OpenMP if available; otherwise std::thread (also used by the
// ThreadSanitizer build, which cannot see inside an uninstrumented OpenMP runtime).
#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <thread>
#include <vector>

#ifdef FG_HAVE_OPENMP
#include <omp.h>
#endif

namespace fg {

/// Threads available to FourierGrating: FG_NUM_THREADS, else the OpenMP maximum
/// (OMP_NUM_THREADS), else the hardware concurrency.
int availableThreads();

namespace detail {
/// Nesting depth of parallelFor on the calling thread.
int& parallelDepth();
struct DepthGuard {
  DepthGuard() { ++parallelDepth(); }
  ~DepthGuard() { --parallelDepth(); }
  DepthGuard(const DepthGuard&) = delete;
  DepthGuard& operator=(const DepthGuard&) = delete;
};
}  // namespace detail

/// Threads for inner loops (normal-field raster, contour pieces, section nodes): all available
/// threads if the caller is not itself a task of parallelFor, else 1 (no oversubscription).
inline int innerThreads() {
  return detail::parallelDepth() > 0 ? 1 : availableThreads();
}

/// Calls f(i) for i = 0..n-1 on up to `threads` threads (dynamic schedule). The first exception
/// (lowest index) is rethrown after all tasks have finished.
template <typename F>
void parallelFor(int n, int threads, F&& f) {
  const int t = std::min(threads, n);
  if (t <= 1) {
    for (int i = 0; i < n; ++i) f(i);
    return;
  }
  std::vector<std::exception_ptr> err(static_cast<std::size_t>(n));
#ifdef FG_HAVE_OPENMP
#pragma omp parallel for num_threads(t) schedule(dynamic, 1)
  for (int i = 0; i < n; ++i) {
    try {
      const detail::DepthGuard guard;
      f(i);
    } catch (...) {
      err[static_cast<std::size_t>(i)] = std::current_exception();
    }
  }
#else
  std::atomic<int> next{0};
  auto work = [&] {
    const detail::DepthGuard guard;
    for (int i = next++; i < n; i = next++) {
      try {
        f(i);
      } catch (...) {
        err[static_cast<std::size_t>(i)] = std::current_exception();
      }
    }
  };
  std::vector<std::thread> pool;
  for (int k = 1; k < t; ++k) pool.emplace_back(work);
  work();
  for (auto& th : pool) th.join();
#endif
  for (const auto& e : err)
    if (e) std::rethrow_exception(e);
}

}  // namespace fg
