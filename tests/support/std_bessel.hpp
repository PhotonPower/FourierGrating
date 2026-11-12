// The C++17 special math functions std::cyl_bessel_j / std::cyl_bessel_k serve as independent
// references in some tests. libc++ (Apple clang) does not provide them; there those comparisons
// are skipped with FG_SKIP_WITHOUT_STD_BESSEL() and still run with libstdc++ and MSVC.
#pragma once

#include <cmath>
#include <limits>

#if defined(__cpp_lib_math_special_functions)
#define FG_HAS_STD_BESSEL 1
#else
#define FG_HAS_STD_BESSEL 0
#endif

namespace fgtest {

inline double cyl_bessel_j(double nu, double x) {
#if FG_HAS_STD_BESSEL
  return std::cyl_bessel_j(nu, x);
#else
  (void)nu;
  (void)x;
  return std::numeric_limits<double>::quiet_NaN();
#endif
}

inline double cyl_bessel_k(double nu, double x) {
#if FG_HAS_STD_BESSEL
  return std::cyl_bessel_k(nu, x);
#else
  (void)nu;
  (void)x;
  return std::numeric_limits<double>::quiet_NaN();
#endif
}

}  // namespace fgtest

#define FG_SKIP_WITHOUT_STD_BESSEL() \
  if (!FG_HAS_STD_BESSEL) GTEST_SKIP() << "std::cyl_bessel_j/k not available (libc++)"
