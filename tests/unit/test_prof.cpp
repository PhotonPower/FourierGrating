// fg::prof (Plan section "M15: Profiling und Parallelisierung", step 1): scope timers, thread-safe
// aggregation, JSON output, negligible cost when disabled.
#include <gtest/gtest.h>

#include <chrono>
#include <thread>

#include "fg/core/prof.hpp"

namespace {

void busy(double seconds) {
  const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
  while (std::chrono::steady_clock::now() < end) {
  }
}

}  // namespace

TEST(Prof, DisabledRecordsNothing) {
  fg::prof::enable(false);
  fg::prof::reset();
  {
    FG_PROF_SCOPE("test.off");
    busy(1e-3);
  }
  EXPECT_TRUE(fg::prof::snapshot().empty());
}

TEST(Prof, AggregatesSumCountMaxOverThreads) {
  fg::prof::enable(true);
  fg::prof::reset();
  auto work = [] {
    for (int k = 0; k < 5; ++k) {
      FG_PROF_SCOPE("test.inner");
      busy(2e-3);
    }
  };
  std::thread a(work), b(work);
  a.join();
  b.join();
  {
    FG_PROF_SCOPE("test.outer");
    busy(1e-3);
  }
  const auto s = fg::prof::snapshot();
  fg::prof::enable(false);
  ASSERT_EQ(s.size(), 2u);
  // sorted by total time, descending
  EXPECT_EQ(s[0].name, "test.inner");
  EXPECT_EQ(s[0].count, 10);
  EXPECT_GE(s[0].total_s, 10 * 2e-3);
  EXPECT_GE(s[0].max_s, 2e-3);
  EXPECT_LE(s[0].max_s, s[0].total_s);
  EXPECT_EQ(s[1].name, "test.outer");
  EXPECT_EQ(s[1].count, 1);
  const std::string j = fg::prof::json();
  EXPECT_NE(j.find("\"test.inner\""), std::string::npos);
  EXPECT_NE(j.find("\"count\": 10"), std::string::npos);
  fg::prof::reset();
  EXPECT_TRUE(fg::prof::snapshot().empty());
}

TEST(Prof, DisabledCostIsNegligible) {
  // < 1 % of a cheap operation scale: a disabled scope is one relaxed atomic load
  fg::prof::enable(false);
  const int n = 10'000'000;
  volatile double acc = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < n; ++k) {
    FG_PROF_SCOPE("test.cost");
    acc = acc + 1.0;
  }
  const double per =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / n;
  EXPECT_LT(per, 20e-9);  // the instrumented scopes cost >= 10 us each (eig, star, coefficients)
}
