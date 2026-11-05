#include "fg/core/quadrature.hpp"

#include <cmath>
#include <map>

#include "fg/core/error.hpp"
#include "fg/core/types.hpp"

namespace fg::quadrature {

const Rule& gaussLegendre(int n) {
  FG_REQUIRE(n >= 1, "Gauss-Legendre: n >= 1");
  thread_local std::map<int, Rule> cache;
  auto it = cache.find(n);
  if (it != cache.end()) return it->second;
  Rule r;
  r.x.resize(static_cast<std::size_t>(n));
  r.w.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < (n + 1) / 2; ++i) {
    double z = std::cos(PI * (i + 0.75) / (n + 0.5)), dp = 1;
    for (int iter = 0; iter < 100; ++iter) {
      double p0 = 1, p1 = z;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1);
      const double dz = p1 / dp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    const double w = 2.0 / ((1 - z * z) * dp * dp);
    // symmetric pairs (exact symmetry of the rule)
    r.x[static_cast<std::size_t>(i)] = -z;
    r.x[static_cast<std::size_t>(n - 1 - i)] = z;
    r.w[static_cast<std::size_t>(i)] = r.w[static_cast<std::size_t>(n - 1 - i)] = w;
  }
  if (n % 2 == 1) r.x[static_cast<std::size_t>(n / 2)] = 0.0;
  return cache.emplace(n, std::move(r)).first->second;
}

}  // namespace fg::quadrature
