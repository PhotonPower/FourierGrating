// FourierGrating - shape derivatives of indicator coefficients (Plan eq:shape2d).
#include "fg/adjoint/shape.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"

namespace fg::adjoint {

namespace {

void contourTerm(const geometry::Contour2D& c, const geometry::Contour2D& cp, const geometry::Contour2D& cm, double delta,
                 double dx, double dz, int Px, int Pz, fourier::CoefficientGrid& out) {
  FG_REQUIRE(c.segments.size() == cp.segments.size() && c.segments.size() == cm.segments.size(),
             "shape derivative: the perturbed contour has another segment structure");
  const auto& g = quadrature::gaussLegendre(16);
  const double kx = 2 * PI / dx, kz = 2 * PI / dz;
  const double kmax = std::hypot(kx * Px, kz * Pz);
  std::vector<cplx> ex(static_cast<std::size_t>(2 * Px + 1)), ez(static_cast<std::size_t>(2 * Pz + 1));
  for (std::size_t s = 0; s < c.segments.size(); ++s) {
    const auto& seg = c.segments[s];
    const int m = std::max(4, static_cast<int>(std::ceil(seg.controlLength() * kmax)));  // <= 1 rad per panel
    for (int panel = 0; panel < m; ++panel)
      for (std::size_t q = 0; q < g.x.size(); ++q) {
        const double t = (panel + 0.5 * (1 + g.x[q])) / m, w = 0.5 * g.w[q] / m;
        geometry::Point2 r, dr;
        seg.eval(t, r, dr);
        const geometry::Point2 rp = cp.segments[s].point(t), rm = cm.segments[s].point(t);
        const double Vx = (rp[0] - rm[0]) / (2 * delta), Vz = (rp[1] - rm[1]) / (2 * delta);
        const double val = w * (Vx * dr[1] - Vz * dr[0]) / (dx * dz);
        if (val == 0.0) continue;
        for (int p = -Px; p <= Px; ++p) ex[static_cast<std::size_t>(p + Px)] = std::exp(-I * (p * kx * r[0]));
        for (int k = -Pz; k <= Pz; ++k) ez[static_cast<std::size_t>(k + Pz)] = std::exp(-I * (k * kz * r[1]));
        for (int p = -Px; p <= Px; ++p) {
          const cplx a = val * ex[static_cast<std::size_t>(p + Px)];
          for (int k = -Pz; k <= Pz; ++k) out(p, k) += a * ez[static_cast<std::size_t>(k + Pz)];
        }
      }
  }
}

}  // namespace

fourier::CoefficientGrid shapeDerivative(const geometry::Region& base, const geometry::Region& plus,
                                         const geometry::Region& minus, double delta, double dx, double dz, int Px,
                                         int Pz) {
  FG_REQUIRE(base.holes.size() == plus.holes.size() && base.holes.size() == minus.holes.size(),
             "shape derivative: the perturbed region has another number of holes");
  fourier::CoefficientGrid out(Px, Pz);
  contourTerm(base.outer, plus.outer, minus.outer, delta, dx, dz, Px, Pz, out);
  for (std::size_t h = 0; h < base.holes.size(); ++h)
    contourTerm(base.holes[h], plus.holes[h], minus.holes[h], delta, dx, dz, Px, Pz, out);
  return out;
}

}  // namespace fg::adjoint
