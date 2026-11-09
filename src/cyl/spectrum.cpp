#include "fg/cyl/spectrum.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"

namespace fg::cyl {

SpectralGrid spectralGrid(const SpectrumOptions& opt) {
  FG_REQUIRE(opt.kmax > 0 && opt.Rs > 0 && opt.dkFactor > 0,
             "spectralGrid: k_max, R_s, dk_factor > 0");
  FG_REQUIRE(opt.nodesPerPanel >= 2, "spectralGrid: at least 2 nodes per panel");
  const bool deform = opt.delta > 0;
  if (deform)
    FG_REQUIRE(0 <= opt.deformA && opt.deformA < opt.deformB && opt.deformB <= opt.kmax,
               "spectralGrid: 0 <= k_a < k_b <= k_max");
  const double tol = 1e-12 * opt.kmax;
  std::vector<double> branch;
  for (double b : opt.branchPoints)
    if (b > tol && b < opt.kmax - tol) branch.push_back(b);
  auto isBranch = [&](double x) {
    return std::any_of(branch.begin(), branch.end(),
                       [&](double b) { return std::abs(b - x) <= tol; });
  };
  std::vector<double> edges{0.0, opt.kmax};
  edges.insert(edges.end(), branch.begin(), branch.end());
  if (deform) {
    edges.push_back(opt.deformA);
    edges.push_back(opt.deformB);
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(
      std::unique(edges.begin(), edges.end(), [&](double a, double b) { return b - a <= tol; }),
      edges.end());

  const double hmax = PI / (2 * opt.Rs * opt.dkFactor);
  const auto& gl = quadrature::gaussLegendre(opt.nodesPerPanel);
  std::vector<cplx> ks, ws;
  for (std::size_t e = 0; e + 1 < edges.size(); ++e) {
    const double a = edges[e], b = edges[e + 1];
    // a panel between two branch points gets two sub-panels (one-sided substitutions; a single
    // k = c + (h/2) sin t over [-pi/2, pi/2] is 1e3 times less accurate with 8 nodes)
    const int nsub = std::max(isBranch(a) && isBranch(b) ? 2 : 1,
                              static_cast<int>(std::ceil((b - a) / hmax - 1e-12)));
    const double h = (b - a) / nsub;
    for (int s = 0; s < nsub; ++s) {
      const double ua = a + s * h, ub = ua + h;
      // square-root substitution on sub-panels that touch a branch point: du/dt ~ sqrt(|u - k_h|)
      const bool lo = s == 0 && isBranch(a), hi = s == nsub - 1 && isBranch(b);
      const double t0 = (lo || hi) ? 0.0 : -1.0, t1 = (lo || hi) ? PI / 2 : 1.0;
      for (std::size_t q = 0; q < gl.x.size(); ++q) {
        const double t = 0.5 * (t1 - t0) * gl.x[q] + 0.5 * (t1 + t0),
                     wt = 0.5 * (t1 - t0) * gl.w[q];
        double u, du;
        if (lo) {
          u = ua + h * (1 - std::cos(t));
          du = h * std::sin(t);
        } else if (hi) {
          u = ua + h * std::sin(t);
          du = h * std::cos(t);
        } else {
          u = 0.5 * (ua + ub) + 0.5 * h * t;
          du = 0.5 * h;
        }
        cplx k = u, dk = du;
        if (deform && u > opt.deformA && u < opt.deformB) {  // Plan sec:cylgrid, complex path
          const double L = opt.deformB - opt.deformA, x = PI * (u - opt.deformA) / L;
          k = cplx(u, -opt.delta * std::sin(x));
          dk = du * cplx(1.0, -opt.delta * PI / L * std::cos(x));
        }
        ks.push_back(k);
        ws.push_back(wt * dk);
      }
    }
  }
  SpectralGrid g;
  g.k = Eigen::Map<VecC>(ks.data(), static_cast<Eigen::Index>(ks.size()));
  g.w = Eigen::Map<VecC>(ws.data(), static_cast<Eigen::Index>(ws.size()));
  g.edges = std::move(edges);
  return g;
}

}  // namespace fg::cyl
