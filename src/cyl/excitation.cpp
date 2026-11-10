#include "fg/cyl/excitation.hpp"

#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/cyl/bessel.hpp"

namespace fg::cyl {

VecC besselProjection(const SpectralGrid& grid, int nu, int nup, const RadialSteps& g, cplx kp) {
  FG_REQUIRE(g.radii.size() == g.value.size() && !g.radii.empty(), "besselProjection: radii and values");
  for (std::size_t i = 0; i < g.radii.size(); ++i)
    FG_REQUIRE(g.radii[i] > (i ? g.radii[i - 1] : 0.0), "besselProjection: radii increasing and > 0");
  // J_{-n} = (-1)^n J_n: Lommel for equal moduli with the sign of both orders
  auto sgn = [](int v) { return (v < 0 && (-v) % 2 == 1) ? -1.0 : 1.0; };
  const int n = std::abs(nu), np = std::abs(nup);
  const double sg = sgn(nu) * sgn(nup);
  const Eigen::Index nk = grid.k.size();
  // int_{r0}^{r1} J_n(a r) J_np(kp r) r dr by composite Gauss-Legendre (panels <= pi / (2 max|k|))
  auto quad = [&](cplx a, double r0, double r1) {
    const auto& gl = quadrature::gaussLegendre(16);
    const double hmax = PI / (2 * std::max({std::abs(a), std::abs(kp), 1e-300}));
    const int panels = std::max(1, static_cast<int>(std::ceil((r1 - r0) / hmax - 1e-12)));
    const double h = (r1 - r0) / panels;
    cplx q = 0.0;
    for (int p = 0; p < panels; ++p)
      for (std::size_t s = 0; s < gl.x.size(); ++s) {
        const double r = r0 + h * (p + 0.5 * (gl.x[s] + 1));
        q += besselJ(n, a * r) * besselJ(np, kp * r) * r * 0.5 * h * gl.w[s];
      }
    return q;
  };
  VecC out = VecC::Zero(nk);
  if (n != np) {
    for (Eigen::Index j = 0; j < nk; ++j)
      for (std::size_t i = 0; i < g.radii.size(); ++i)
        out(j) += g.value[i] * quad(grid.k(j), i ? g.radii[i - 1] : 0.0, g.radii[i]);
    return sg * out;
  }
  auto der = [n](const std::vector<cplx>& v) {
    return 0.5 * ((n == 0 ? -v[1] : v[static_cast<std::size_t>(n - 1)]) - v[static_cast<std::size_t>(n + 1)]);
  };
  for (Eigen::Index j = 0; j < nk; ++j) {
    const cplx a = grid.k(j);
    cplx sum = 0.0, F0 = 0.0;
    for (std::size_t i = 0; i < g.radii.size(); ++i) {
      const double r1 = g.radii[i];
      cplx F1;
      if (std::abs(a - kp) * r1 < 1e-2) {  // near-coincident pair: the closed form cancels
        F1 = quad(a, 0.0, r1);
      } else {
        const auto va = besselJAll(n + 1, a * r1), vb = besselJAll(n + 1, kp * r1);
        const cplx Ja = va[static_cast<std::size_t>(n)], Jb = vb[static_cast<std::size_t>(n)];
        F1 = r1 * (kp * Ja * der(vb) - a * der(va) * Jb) / (a * a - kp * kp);  // Plan eq:lommel
      }
      sum += g.value[i] * (F1 - F0);
      F0 = F1;
    }
    out(j) = sg * sum;
  }
  return out;
}

PlaneWaveSource planeWaveSource(const SpectralGrid& grid, int m, const Rings& rings, double kpar, double Psi,
                                CylFactorization f, const CylQ& Q) {
  const Eigen::Index n = grid.k.size();
  const int p = m + 1, q = m - 1;
  // i^nu e^{-i nu Psi}: E_bg,pm = c_pm phase(m +- 1) J_{m +- 1}(kpar r) (Plan eq:cyl_jacobi)
  auto phase = [Psi](int nu) { return std::polar(1.0, nu * (0.5 * PI - Psi)); };
  RadialSteps d{rings.radii, {}};  // eps - eps_inf
  for (cplx e : rings.eps) d.value.push_back(e - rings.epsInf);
  PlaneWaveSource B;
  B.By = phase(m) * besselProjection(grid, m, m, d, kpar);
  B.Bpm = MatC::Zero(2 * n, 2);
  if (f == CylFactorization::Laurent) {
    B.Bpm.col(0).head(n) = phase(p) * besselProjection(grid, p, p, d, kpar);
    B.Bpm.col(1).tail(n) = phase(q) * besselProjection(grid, q, q, d, kpar);
    return B;
  }
  // Li in r: P = b^ - Q a^ with the columns of A - 1 = [[a]] and B - eps_inf = [[b]] (signs of Plan
  // eq:cyl_li as in cylQ) applied to E_bg
  RadialSteps ra{rings.radii, {}}, rb{rings.radii, {}};
  for (cplx e : rings.eps) {
    ra.value.push_back(0.5 * (rings.epsInf / e - 1.0));
    rb.value.push_back(0.5 * (e - rings.epsInf));
  }
  MatC ah(2 * n, 2), bh(2 * n, 2);
  ah << besselProjection(grid, p, p, ra, kpar), besselProjection(grid, p, q, ra, kpar),
      besselProjection(grid, q, p, ra, kpar), besselProjection(grid, q, q, ra, kpar);
  bh << besselProjection(grid, p, p, rb, kpar), -besselProjection(grid, p, q, rb, kpar),
      -besselProjection(grid, q, p, rb, kpar), besselProjection(grid, q, q, rb, kpar);
  const Eigen::Vector2cd ph(phase(p), phase(q));
  ah = ah * ph.asDiagonal();
  bh = bh * ph.asDiagonal();
  MatC Qf(2 * n, 2 * n);
  Qf << Q.Qpp, Q.Qpm, Q.Qmp, Q.Qmm;
  B.Bpm = bh - Qf * ah;
  return B;
}

VecC scatteredSource(const VecC& k, const Eigen::PartialPivLU<MatC>& Qyy, const PlaneWaveSource& B, cplx cp,
                     cplx cm, cplx ey) {
  const Eigen::Index n = k.size();
  const VecC kEpy = k.cwiseProduct(Qyy.solve(VecC(ey * B.By)));
  VecC s(4 * n);
  const VecC ppm = B.Bpm * Eigen::Vector2cd(cp, cm);
  s << kEpy, -kEpy, ppm.head(n), -ppm.tail(n);  // Plan eq:cyl_scat
  return s;
}

int harmonicCutoff(double x) {
  return static_cast<int>(std::ceil(x + 4.0 * std::cbrt(x) + 2.0));
}

}  // namespace fg::cyl
