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
  // delta < 0 deforms into Im k > 0, past the poles of the guided modes: the wrong side (V21f shows the jump;
  // the configuration only accepts delta > 0)
  const bool deform = opt.delta != 0;
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
  const double hmax = PI / (2 * opt.Rs * opt.dkFactor);
  // graded panels around poles, between the neighbouring fixed edges: resonances on the real section (width
  // Im p) and bound modes on the deformed section (width = distance of the pole to the path)
  {
    std::vector<double> fixed = edges;
    std::sort(fixed.begin(), fixed.end());
    auto grade = [&](double c, double width) {
      const auto hiIt = std::upper_bound(fixed.begin(), fixed.end(), c);
      const double lo = *(hiIt - 1), hi = *hiIt, w0 = std::max(width, 1e-3 * hmax);
      for (double w = w0; w < hmax; w *= 2)
        for (const double e : {c - w, c + w})
          if (e > lo + 0.5 * w0 && e < hi - 0.5 * w0) edges.push_back(e);
    };
    for (const cplx p : opt.resonances) {
      const double c = p.real();
      if (c <= tol || c >= opt.kmax - tol || (deform && c > opt.deformA - tol && c < opt.deformB + tol)) continue;
      grade(c, std::abs(p.imag()));
    }
    if (deform)
      for (const double c : opt.boundModes) {
        if (c <= opt.deformA + tol || c >= opt.deformB - tol) continue;
        const double L = opt.deformB - opt.deformA;
        grade(c, std::abs(opt.delta) * std::sin(PI * (c - opt.deformA) / L));
      }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(
      std::unique(edges.begin(), edges.end(), [&](double a, double b) { return b - a <= tol; }),
      edges.end());

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

namespace {

/// Characteristic function of the bound modes at real k (Plan sec:cylgrid, v2.38): (u, w) = (E_z, E_z' / p)
/// resp. (H_z, H_z' / (p = eps)) from the decaying substrate solution through the layer matrices; zero when
/// the field also decays into the superstrate. Real for lossless materials; scaled by positive factors only.
double boundModeFunction(double k, bool tm, double epsSub, double epsSup, const std::vector<double>& eps,
                         const std::vector<double>& h) {
  const double gs = std::sqrt(std::max(0.0, k * k - epsSub)), gc = std::sqrt(std::max(0.0, k * k - epsSup));
  double u = 1, w = gs / (tm ? epsSub : 1.0);
  for (std::size_t j = 0; j < eps.size(); ++j) {
    const double p = tm ? eps[j] : 1.0, b2 = eps[j] - k * k;
    double c, s, bs;  // cos(beta h), sin(beta h) / beta, beta sin(beta h), entire in beta^2
    if (b2 > 0) {
      const double b = std::sqrt(b2);
      c = std::cos(b * h[j]);
      s = std::sin(b * h[j]) / b;
      bs = b * std::sin(b * h[j]);
    } else if (b2 < 0) {
      const double q = std::sqrt(-b2);
      c = std::cosh(q * h[j]);
      s = std::sinh(q * h[j]) / q;
      bs = -q * std::sinh(q * h[j]);
    } else {
      c = 1, s = h[j], bs = 0;
    }
    const double un = c * u + p * s * w, wn = -bs / p * u + c * w;
    const double sc = std::max(std::abs(un), std::abs(wn));
    u = un / sc;
    w = wn / sc;
  }
  return w + gc / (tm ? epsSup : 1.0) * u;
}

}  // namespace

std::vector<double> planarBoundModes(cplx epsSub, cplx epsSup, const std::vector<cplx>& epsLayers,
                                     const std::vector<double>& thickness, double ka, double kb) {
  FG_REQUIRE(epsLayers.size() == thickness.size(), "planarBoundModes: one thickness per layer");
  auto lossless = [](cplx e) { return std::abs(e.imag()) <= 1e-14 * std::abs(e) && e.real() != 0; };
  if (!lossless(epsSub) || !lossless(epsSup) || epsSub.real() < 0 || epsSup.real() < 0 || !(ka < kb)) return {};
  std::vector<double> eps;
  for (const cplx e : epsLayers) {
    if (!lossless(e)) return {};
    eps.push_back(e.real());
  }
  FG_REQUIRE(ka * ka >= std::max(epsSub.real(), epsSup.real()) * (1 - 1e-12),
             "planarBoundModes: k_a below a branch point of the half spaces");
  std::vector<double> roots;
  const int n = 4096;  // sign changes on a fine grid, then bisection
  for (const bool tm : {false, true}) {
    auto f = [&](double k) { return boundModeFunction(k, tm, epsSub.real(), epsSup.real(), eps, thickness); };
    double k0 = ka + 1e-9 * (kb - ka), f0 = f(k0);
    for (int i = 1; i <= n; ++i) {
      const double k1 = ka + (kb - ka) * (i == n ? 1 - 1e-9 : double(i) / n), f1 = f(k1);
      if (f0 == 0) roots.push_back(k0);
      if (f0 * f1 < 0) {
        double a = k0, b = k1, fa = f0;
        for (int it = 0; it < 200 && b - a > 1e-15 * kb; ++it) {
          const double m = 0.5 * (a + b), fm = f(m);
          if (fa * fm <= 0)
            b = m;
          else
            a = m, fa = fm;
        }
        roots.push_back(0.5 * (a + b));
      }
      k0 = k1, f0 = f1;
    }
  }
  std::sort(roots.begin(), roots.end());
  return roots;
}

}  // namespace fg::cyl
