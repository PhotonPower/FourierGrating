// FourierGrating - Fourier coefficients from contour integrals (Plan eq:contour, eq:contour0).
#include "fg/fourier/contour.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "fg/core/error.hpp"
#include "fg/core/parallel.hpp"
#include "fg/core/prof.hpp"
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

  /// Accepted (converged) piece sums of segment s in contour order (adaptive bisection).
  void pieces(const geometry::BezierSegment& s, int depth, std::vector<PieceSums>& out) const {
    const int n = n0 > 0 ? n0 : std::max(16, static_cast<int>(std::ceil(0.5 * kmax * s.controlLength())) + 8);
    PieceSums a = pieceSums(s, n, Kx, Kz, Px, Pz), b = pieceSums(s, 2 * n, Kx, Kz, Px, Pz);
    if (difference(a, b) > tolAbs && depth < 30) {
      geometry::BezierSegment l, r;
      s.split(0.5, l, r);
      pieces(l, depth + 1, out);
      pieces(r, depth + 1, out);
      return;
    }
    out.push_back(std::move(b));
  }
  /// The segments are independent (Plan section "M15", step 4); the sums are added in contour order.
  void add(const std::vector<const geometry::BezierSegment*>& segs) {
    std::vector<std::vector<PieceSums>> parts(segs.size());
    parallelFor(static_cast<int>(segs.size()), innerThreads(),
                [&](int k) { pieces(*segs[static_cast<std::size_t>(k)], 0, parts[static_cast<std::size_t>(k)]); });
    for (const auto& part : parts)
      for (const PieceSums& b : part) {
        G1 += b.G1;
        G2 += b.G2;
        area += b.area;
      }
  }
};

/// Piece [t0, t1] of a Bezier segment mapped to transformed coordinates.
struct MappedPiece {
  const geometry::BezierSegment* s;
  double t0, t1;
};

struct MappedAccumulator {
  const CoordinateMap *xm, *zm;
  double Kx, Kz, kmax, tolAbs, stretch;
  int Px, Pz;
  MatC G1, G2;
  double area = 0;

  PieceSums sums(const MappedPiece& pc, int n) const {
    const auto& g = quadrature::gaussLegendre(n);
    const int np = 2 * Px + 1, nq = 2 * Pz + 1;
    MatC Ex(n, np), Ez(n, nq);
    VecC wz(n), wx(n);
    PieceSums out;
    const double h = pc.t1 - pc.t0;
    for (int l = 0; l < n; ++l) {
      geometry::Point2 r, dr;
      pc.s->eval(pc.t0 + 0.5 * h * (g.x[static_cast<std::size_t>(l)] + 1), r, dr);
      // mapped point and derivative d(X, Z)/dt
      const double X = xm ? xm->inverse(r[0]) : r[0], Z = zm ? zm->inverse(r[1]) : r[1];
      const double dX = (xm ? dr[0] / xm->df(X) : dr[0]) * h, dZ = (zm ? dr[1] / zm->df(Z) : dr[1]) * h;
      const double w = 0.5 * g.w[static_cast<std::size_t>(l)];
      wz(l) = w * dZ;
      wx(l) = w * dX;
      out.area += 0.5 * w * (X * dZ - Z * dX);  // Plan eq:contour0
      const cplx ex = std::exp(-I * (Kx * X)), ez = std::exp(-I * (Kz * Z));
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

  void pieces(const MappedPiece& pc, int depth, std::vector<PieceSums>& out) const {
    const double len = pc.s->controlLength() * (pc.t1 - pc.t0) * stretch;
    const int n = std::max(16, static_cast<int>(std::ceil(0.5 * kmax * len)) + 8);
    PieceSums a = sums(pc, n), b = sums(pc, 2 * n);
    if (difference(a, b) > tolAbs && depth < 40) {
      const double tm = 0.5 * (pc.t0 + pc.t1);
      pieces({pc.s, pc.t0, tm}, depth + 1, out);
      pieces({pc.s, tm, pc.t1}, depth + 1, out);
      return;
    }
    out.push_back(std::move(b));
  }
};

}  // namespace

CoefficientGrid regionCoefficients(const geometry::Region& region, double dx, double dz, int Px, int Pz,
                                   const CoordinateMap* xmap, const CoordinateMap* zmap) {
  FG_PROF_SCOPE("contour.coefficients");
  FG_REQUIRE(dx > 0 && dz > 0 && Px >= 0 && Pz >= 0, "regionCoefficients: invalid period or orders");
  const double Kx = 2 * PI / dx, Kz = 2 * PI / dz;
  // the mapped boundary is longer by at most 1 / min f' = 1 / (1 - eta)
  const double stretch = 1.0 / std::min(xmap ? 1.0 - xmap->eta() : 1.0, zmap ? 1.0 - zmap->eta() : 1.0);
  MappedAccumulator acc{xmap, zmap, Kx, Kz, std::hypot(Px * Kx, Pz * Kz), 1e-14 * std::abs(region.area()), stretch,
                        Px, Pz, MatC::Zero(2 * Px + 1, 2 * Pz + 1), MatC::Zero(2 * Px + 1, 2 * Pz + 1), 0.0};
  std::vector<const geometry::BezierSegment*> segs;
  for (const auto& s : region.outer.segments) segs.push_back(&s);
  for (const auto& h : region.holes)
    for (const auto& s : h.segments) segs.push_back(&s);
  std::vector<std::vector<PieceSums>> parts(segs.size());
  parallelFor(static_cast<int>(segs.size()), innerThreads(), [&](int k) {
    acc.pieces({segs[static_cast<std::size_t>(k)], 0.0, 1.0}, 0, parts[static_cast<std::size_t>(k)]);
  });
  for (const auto& part : parts)
    for (const PieceSums& b : part) {
      acc.G1 += b.G1;
      acc.G2 += b.G2;
      acc.area += b.area;
    }
  CoefficientGrid g(Px, Pz);
  const double cell = dx * dz;
  for (int p = -Px; p <= Px; ++p)
    for (int q = -Pz; q <= Pz; ++q) {
      if (p == 0 && q == 0) {
        g(0, 0) = acc.area / cell;  // Plan eq:contour0
        continue;
      }
      const double kx = p * Kx, kz = q * Kz;
      g(p, q) = I * (kx * acc.G1(p + Px, q + Pz) - kz * acc.G2(p + Px, q + Pz)) / ((kx * kx + kz * kz) * cell);
    }
  return g;
}

CoefficientGrid regionCoefficients(const geometry::Region& region, double dx, double dz, int Px, int Pz, int n0) {
  FG_PROF_SCOPE("contour.coefficients");
  FG_REQUIRE(dx > 0 && dz > 0 && Px >= 0 && Pz >= 0, "regionCoefficients: invalid period or orders");
  const double Kx = 2 * PI / dx, Kz = 2 * PI / dz;
  Accumulator acc{Kx, Kz, std::hypot(Px * Kx, Pz * Kz), 1e-14 * std::abs(region.area()), Px, Pz, n0,
                  MatC::Zero(2 * Px + 1, 2 * Pz + 1), MatC::Zero(2 * Px + 1, 2 * Pz + 1), 0.0};
  std::vector<const geometry::BezierSegment*> segs;
  for (const auto& s : region.outer.segments) segs.push_back(&s);
  for (const auto& h : region.holes)
    for (const auto& s : h.segments) segs.push_back(&s);
  acc.add(segs);
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
