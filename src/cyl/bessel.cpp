#include "fg/cyl/bessel.hpp"

#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"

namespace fg::cyl {

std::vector<cplx> besselJAll(int nmax, cplx z) {
  FG_REQUIRE(nmax >= 0, "besselJ: nmax >= 0");
  FG_REQUIRE(std::abs(z.imag()) <= 700.0, "besselJ: |Im z| <= 700");
  std::vector<cplx> J(static_cast<std::size_t>(nmax) + 1, 0.0);
  const double az = std::abs(z);
  if (az < 1e-6) {  // two terms of the power series (the recurrence would overflow per step)
    const cplx h = 0.5 * z, h2 = h * h;
    cplx t = 1.0;
    for (int n = 0; n <= nmax; ++n) {
      J[static_cast<std::size_t>(n)] = t * (1.0 - h2 / double(n + 1));
      t *= h / double(n + 1);
    }
    return J;
  }
  // Start index beyond the Airy transition n ~ |z|: J_n ~ (2/x)^(1/3) Ai((n - x)(2/x)^(1/3))
  // is below 1e-17 of its peak for n - x >= 15.2 (x/2)^(1/3).
  int N = static_cast<int>(std::ceil(std::max<double>(nmax, az) + 16.0 * std::cbrt(az / 2) + 20.0));
  N += N % 2;
  // Normalization exp(-s i z) = J_0 + 2 sum_n (-s i)^n J_n with s = sign(Im z) (grows like J);
  // real z: 1 = J_0 + 2 sum_n J_2n, which keeps the result exactly real.
  const bool real = z.imag() == 0.0;
  const cplx unit = z.imag() >= 0 ? cplx(0, -1) : cplx(0, 1);
  cplx pw[4] = {1.0, unit, unit * unit, unit * unit * unit};
  if (real) pw[1] = pw[3] = 0.0, pw[2] = 1.0;
  cplx fp1 = 0.0, f = 1e-30, sum = 0.0;
  const cplx zinv = 1.0 / z;
  for (int n = N; n >= 0; --n) {
    if (n <= nmax) J[static_cast<std::size_t>(n)] = f;
    sum += (n == 0 ? 1.0 : 2.0) * pw[n % 4] * f;
    if (n == 0) break;
    const cplx fm1 = 2.0 * n * zinv * f - fp1;
    fp1 = f;
    f = fm1;
    if (std::abs(f) > 1e200) {  // rescale to avoid overflow
      const double s = 1e-200;
      f *= s;
      fp1 *= s;
      sum *= s;
      for (int k = n; k <= nmax; ++k) J[static_cast<std::size_t>(k)] *= s;
    }
  }
  const cplx norm = (real ? cplx(1.0) : std::exp(unit * z)) / sum;
  for (auto& v : J) v *= norm;
  return J;
}

cplx besselJ(int n, cplx z) {
  const int a = std::abs(n);
  const cplx v = besselJAll(a, z)[static_cast<std::size_t>(a)];
  return (n < 0 && a % 2 == 1) ? -v : v;
}

cplx besselJp(int n, cplx z) {
  return 0.5 * (besselJ(n - 1, z) - besselJ(n + 1, z));
}

}  // namespace fg::cyl
