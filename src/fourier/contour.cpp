// FourierGrating - Fourier coefficients from contour integrals (Plan eq:contour, eq:contour0).
#include "fg/fourier/contour.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"

namespace fg::fourier {

namespace {

/// Contour sums of one piece: G1_pq = sum w z' e^{-i(kx x + kz z)}, G2 with x', and the area term.
struct PieceSums {
  MatC G1, G2;
  double area = 0;
};

PieceSums pieceSums(const geometry::BezierSegment& s, int n, double Kx, double Kz, int Px, int Pz) {
  const auto& g = quadrature::gaussLegendre(n);
  const int np = 2 * Px + 1, nq = 2 * Pz + 1;
  MatC Ex(n, np), Ez(n, nq);
  VecC wz(n), wx(n);
  PieceSums out;
  for (int l = 0; l < n; ++l) {
    geometry::Point2 r, dr;
    s.eval(0.5 * (g.x[static_cast<std::size_t>(l)] + 1), r, dr);
    const double w = 0.5 * g.w[static_cast<std::size_t>(l)];
    wz(l) = w * dr[1];
    wx(l) = w * dr[0];
    out.area += 0.5 * w * (r[0] * dr[1] - r[1] * dr[0]);  // Plan eq:contour0
    // E^x_lp = e^{-i p Kx x}, powers from p = 0 outwards (rounding ~ |p| eps)
    const cplx ex = std::exp(-I * (Kx * r[0])), ez = std::exp(-I * (Kz * r[1]));
    Ex(l, Px) = 1.0;
    for (int p = 1; p <= Px; ++p) {
      Ex(l, Px + p) = Ex(l, Px + p - 1) * ex;
      Ex(l, Px - p) = Ex(l, Px - p + 1) * std::conj(ex);
    }
    Ez(l, Pz) = 1.0;
    for (int q = 1; q <= Pz; ++q) {
      Ez(l, Pz + q) = Ez(l, Pz + q - 1) * ez;
      Ez(l, Pz - q) = Ez(l, Pz - q + 1) * std::conj(ez);
    }
  }
  out.G1 = Ex.transpose() * (wz.asDiagonal() * Ez);
  out.G2 = Ex.transpose() * (wx.asDiagonal() * Ez);
  return out;
}

struct Accumulator {
  double Kx, Kz, kmax, tolAbs;
  int Px, Pz, n0;
  MatC G1, G2;
  double area = 0;

  /// max_pq |i (kx dG1 - kz dG2)| / |k|^2 (contribution to dx dz f_pq), plus the area difference
  double difference(const PieceSums& a, const PieceSums& b) const {
    double d = std::abs(a.area - b.area);
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q) {
        if (p == 0 && q == 0) continue;
        const double kx = p * Kx, kz = q * Kz;
        const cplx v = kx * (a.G1(p + Px, q + Pz) - b.G1(p + Px, q + Pz)) - kz * (a.G2(p + Px, q + Pz) - b.G2(p + Px, q + Pz));
        d = std::max(d, std::abs(v) / (kx * kx + kz * kz));
      }
    return d;
  }

  void add(const geometry::BezierSegment& s, int depth) {
    const int n = n0 > 0 ? n0 : std::max(16, static_cast<int>(std::ceil(0.5 * kmax * s.controlLength())) + 8);
    const PieceSums a = pieceSums(s, n, Kx, Kz, Px, Pz), b = pieceSums(s, 2 * n, Kx, Kz, Px, Pz);
    if (difference(a, b) > tolAbs && depth < 30) {
      geometry::BezierSegment l, r;
      s.split(0.5, l, r);
      add(l, depth + 1);
      add(r, depth + 1);
      return;
    }
    G1 += b.G1;
    G2 += b.G2;
    area += b.area;
  }
};

}  // namespace

CoefficientGrid regionCoefficients(const geometry::Region& region, double dx, double dz, int Px, int Pz, int n0) {
  FG_REQUIRE(dx > 0 && dz > 0 && Px >= 0 && Pz >= 0, "regionCoefficients: invalid period or orders");
  const double Kx = 2 * PI / dx, Kz = 2 * PI / dz;
  Accumulator acc{Kx, Kz, std::hypot(Px * Kx, Pz * Kz), 1e-14 * std::abs(region.area()), Px, Pz, n0,
                  MatC::Zero(2 * Px + 1, 2 * Pz + 1), MatC::Zero(2 * Px + 1, 2 * Pz + 1), 0.0};
  for (const auto& s : region.outer.segments) acc.add(s, 0);
  for (const auto& h : region.holes)
    for (const auto& s : h.segments) acc.add(s, 0);
  CoefficientGrid g(Px, Pz);
  const double cell = dx * dz;
  for (int p = -Px; p <= Px; ++p)
    for (int q = -Pz; q <= Pz; ++q) {
      if (p == 0 && q == 0) {
        g(0, 0) = acc.area / cell;  // Plan eq:contour0
        continue;
      }
      const double kx = p * Kx, kz = q * Kz;
      // Plan eq:contour: f_pq = i (kx G1 - kz G2) / (|k|^2 dx dz)
      g(p, q) = I * (kx * acc.G1(p + Px, q + Pz) - kz * acc.G2(p + Px, q + Pz)) / ((kx * kx + kz * kz) * cell);
    }
  return g;
}

}  // namespace fg::fourier
