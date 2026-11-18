#include "fg/cyl/field.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/cyl/bessel.hpp"

namespace fg::cyl {

MatC besselTable(const VecC& k, int nmax, double r) {
  FG_REQUIRE(nmax >= 0, "besselTable: nmax >= 0");
  MatC t(nmax + 1, k.size());
  for (Eigen::Index j = 0; j < k.size(); ++j) {
    const std::vector<cplx> v = besselJAll(nmax, k(j) * r);
    for (int n = 0; n <= nmax; ++n) t(n, j) = v[static_cast<std::size_t>(n)];
  }
  return t;
}

Synthesis synthesis(const SpectralGrid& grid, int m, const std::vector<double>& r) {
  const Eigen::Index nk = grid.k.size(), nr = static_cast<Eigen::Index>(r.size());
  const int nmax = std::abs(m) + 1;
  const VecC wk = grid.w.cwiseProduct(grid.k);
  Synthesis S{MatC(nr, nk), MatC(nr, nk), MatC(nr, nk)};
  for (Eigen::Index s = 0; s < nr; ++s) {
    const MatC t = besselTable(grid.k, nmax, r[static_cast<std::size_t>(s)]);
    for (Eigen::Index j = 0; j < nk; ++j) {
      S.Tm(s, j) = wk(j) * tableJ(t, m - 1, j);
      S.T0(s, j) = wk(j) * tableJ(t, m, j);
      S.Tp(s, j) = wk(j) * tableJ(t, m + 1, j);
    }
  }
  return S;
}

RadialRule radialRule(const std::vector<double>& radii, double kmax, int nodes) {
  FG_REQUIRE(!radii.empty() && kmax > 0, "radialRule: radii and k_max > 0");
  const auto& gl = quadrature::gaussLegendre(nodes);
  const double hmax = PI / (2 * kmax);
  RadialRule q;
  double r0 = 0;
  for (std::size_t i = 0; i < radii.size(); ++i) {
    const double r1 = radii[i];
    FG_REQUIRE(r1 > r0, "radialRule: radii increasing and > 0");
    const int panels = std::max(1, static_cast<int>(std::ceil((r1 - r0) / hmax - 1e-12)));
    const double h = (r1 - r0) / panels;
    for (int p = 0; p < panels; ++p)
      for (std::size_t s = 0; s < gl.x.size(); ++s) {
        const double r = r0 + h * (p + 0.5 * (gl.x[s] + 1));
        q.r.push_back(r);
        q.w.push_back(0.5 * h * gl.w[s] * r);
        q.ring.push_back(i);
      }
    r0 = r1;
  }
  return q;
}

}  // namespace fg::cyl
