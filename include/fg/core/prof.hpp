// FourierGrating - scope profiler (Plan section "M15: Profiling und Parallelisierung", step 1).
// FG_PROF_SCOPE("rcwa.eig") measures the enclosing scope when profiling is enabled; the times are
// aggregated per name over all threads (sum, count, maximum). Disabled: one atomic load per scope.
#pragma once

#include <atomic>
#include <chrono>
#include <string>
#include <vector>

namespace fg::prof {

struct Entry {
  std::string name;
  double total_s = 0;  ///< summed over all threads (busy time, not wall time)
  long count = 0;
  double max_s = 0;
};

void enable(bool on);
bool enabled();
/// Clears all entries.
void reset();
/// Entries sorted by total time, descending.
std::vector<Entry> snapshot();
/// {"name": {"total_s": .., "count": .., "max_s": ..}, ...}
std::string json();
/// Adds one measurement (thread-safe).
void record(const char* name, double seconds);

namespace detail {
extern std::atomic<bool> on;
}

class Scope {
 public:
  explicit Scope(const char* name)
      : name_(detail::on.load(std::memory_order_relaxed) ? name : nullptr) {
    if (name_) t0_ = std::chrono::steady_clock::now();
  }
  ~Scope() {
    if (name_)
      record(name_, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count());
  }
  Scope(const Scope&) = delete;
  Scope& operator=(const Scope&) = delete;

 private:
  const char* name_;
  std::chrono::steady_clock::time_point t0_{};
};

}  // namespace fg::prof

#define FG_PROF_CAT2(a, b) a##b
#define FG_PROF_CAT(a, b) FG_PROF_CAT2(a, b)
#define FG_PROF_SCOPE(name) const ::fg::prof::Scope FG_PROF_CAT(fg_prof_scope_, __LINE__)(name)
