// FourierGrating - angular coefficients of cross sections (M21, Plan section "Winkelkoeffizienten aus
// Querschnitten", eq:cylP_chi).
#include "fg/cyl/angular.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "fg/core/error.hpp"
#include "fg/geometry/bernstein.hpp"

namespace fg::cyl {

using geometry::BezierSegment;
using geometry::Contour2D;
using geometry::Point2;
using geometry::Region;

namespace {

/// Homogeneous coordinates (w x, w z, w) of a segment as Bernstein coefficients.
struct Homogeneous {
  std::vector<double> X, Z, W;
};

Homogeneous homogeneous(const BezierSegment& s) {
  Homogeneous h;
  for (std::size_t j = 0; j < s.P.size(); ++j) {
    const double w = s.weight(j);
    h.X.push_back(w * s.P[j][0]);
    h.Z.push_back(w * s.P[j][1]);
    h.W.push_back(w);
  }
  return h;
}

std::vector<double> add(std::vector<double> a, const std::vector<double>& b, double f = 1.0) {
  for (std::size_t j = 0; j < a.size(); ++j) a[j] += f * b[j];
  return a;
}

std::vector<const Contour2D*> contours(const Region& region) {
  std::vector<const Contour2D*> c{&region.outer};
  for (const auto& h : region.holes) c.push_back(&h);
  return c;
}

double wrap(double phi) {
  phi = std::fmod(phi, 2 * PI);
  return phi < 0 ? phi + 2 * PI : phi;
}

/// Distance of q to a segment: nearest of 64 samples, then golden-section search on the neighbouring interval.
double distanceToSegment(const BezierSegment& s, const Point2& q) {
  auto d2 = [&](double t) {
    const Point2 p = s.point(t);
    return (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
  };
  const int n = 64;
  int best = 0;
  double bestD = d2(0.0);
  for (int k = 1; k <= n; ++k) {
    const double d = d2(double(k) / n);
    if (d < bestD) bestD = d, best = k;
  }
  double a = std::max(0, best - 1) / double(n), b = std::min(n, best + 1) / double(n);
  const double g = 0.5 * (std::sqrt(5.0) - 1);
  double c = b - g * (b - a), e = a + g * (b - a), fc = d2(c), fe = d2(e);
  for (int it = 0; it < 100 && b - a > 1e-16; ++it) {
    if (fc < fe) {
      b = e, e = c, fe = fc, c = b - g * (b - a), fc = d2(c);
    } else {
      a = c, c = e, fc = fe, e = a + g * (b - a), fe = d2(e);
    }
  }
  return std::sqrt(std::min({bestD, fc, fe}));
}

}  // namespace

std::vector<Arc> arcsInside(const Region& region, double r) {
  FG_REQUIRE(r > 0, "arcsInside: r > 0");
  struct Crossing {
    double phi;
    bool enter;
  };
  std::vector<Crossing> xs;
  for (const Contour2D* c : contours(region))
    for (const BezierSegment& s : c->segments) {
      const Homogeneous h = homogeneous(s);
      const auto f = add(add(geometry::bernsteinProduct(h.X, h.X), geometry::bernsteinProduct(h.Z, h.Z)),
                         geometry::bernsteinProduct(h.W, h.W), -r * r);
      for (const double t : geometry::bernsteinRoots(f)) {
        Point2 p, d;
        s.eval(t, p, d);
        const double phi = wrap(std::atan2(p[1], p[0]));
        // the circle runs along e_phi = (-sin, cos); it enters the material if e_phi points to the left of the
        // boundary tangent d (material on the left after normalize)
        const double left = -std::sin(phi) * (-d[1]) + std::cos(phi) * d[0];
        xs.push_back({phi, left > 0});
      }
    }
  if (xs.empty()) {
    if (region.contains(r, 0.0)) return {{0.0, 2 * PI}};
    return {};
  }
  std::sort(xs.begin(), xs.end(), [](const Crossing& a, const Crossing& b) { return a.phi < b.phi; });
  std::vector<Arc> arcs;
  const std::size_t n = xs.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Crossing& a = xs[i];
    const Crossing& b = xs[(i + 1) % n];
    FG_REQUIRE(a.enter != b.enter, "arcsInside: crossings of the circle r = " + std::to_string(r) +
                                       " do not alternate (region not normalized or self-intersecting?)");
    if (a.enter) arcs.push_back({a.phi, i + 1 < n ? b.phi : b.phi + 2 * PI});
  }
  return arcs;
}

VecC angularCoefficients(const Region& region, double r, int P, int qmax) {
  FG_REQUIRE(P >= 1 && qmax >= 0, "angularCoefficients: P >= 1, qmax >= 0");
  VecC chi = VecC::Zero(qmax + 1);
  for (const Arc& a : arcsInside(region, r)) {
    chi(0) += (a.phiB - a.phiA) / (2 * PI);
    for (int q = 1; q <= qmax; ++q) {
      const double n = double(q) * P;
      // 1/(-i n) = i/n (a real division; MSVC flags the complex division by (0, -n) as C4723)
      chi(q) += cplx(0, 1) * (std::exp(cplx(0, -n * a.phiB)) - std::exp(cplx(0, -n * a.phiA))) / (2 * PI * n);
    }
  }
  return chi;
}

VecC angularCoefficientsShape(const Region& region, double r, int P, int qmax, const BoundaryVelocity& velocity) {
  FG_REQUIRE(r > 0 && P >= 1 && qmax >= 0, "angularCoefficientsShape: r > 0, P >= 1, qmax >= 0");
  VecC d = VecC::Zero(qmax + 1);
  const auto cs = contours(region);
  for (std::size_t c = 0; c < cs.size(); ++c)
    for (std::size_t s = 0; s < cs[c]->segments.size(); ++s) {
      const BezierSegment& seg = cs[c]->segments[s];
      const Homogeneous h = homogeneous(seg);
      const auto f = add(add(geometry::bernsteinProduct(h.X, h.X), geometry::bernsteinProduct(h.Z, h.Z)),
                         geometry::bernsteinProduct(h.W, h.W), -r * r);
      for (const double t : geometry::bernsteinRoots(f)) {
        Point2 p, dp;
        seg.eval(t, p, dp);
        const double phi = std::atan2(p[1], p[0]);
        const double left = -std::sin(phi) * (-dp[1]) + std::cos(phi) * dp[0];  // entering (as arcsInside)
        const Point2 V = velocity(c, s, t);
        const double pt = p[0] * dp[0] + p[1] * dp[1];
        FG_REQUIRE(pt != 0.0, "angularCoefficientsShape: r is a tangent radius");
        const double dt = -(p[0] * V[0] + p[1] * V[1]) / pt;
        const double dphi = (p[0] * (V[1] + dp[1] * dt) - p[1] * (V[0] + dp[0] * dt)) / (r * r);
        const double sg = left > 0 ? -1.0 : 1.0;
        for (int q = 0; q <= qmax; ++q)
          d(q) += sg * std::exp(cplx(0, -double(q) * P * phi)) * dphi / (2 * PI);
      }
    }
  return d;
}

std::vector<double> criticalRadii(const Region& region) {
  std::vector<double> radii{0.0};
  double rmax = 0;
  for (const Contour2D* c : contours(region))
    for (const BezierSegment& s : c->segments) {
      radii.push_back(std::hypot(s.P.front()[0], s.P.front()[1]));
      rmax = std::max(rmax, radii.back());
      if (s.degree() < 1) continue;
      // d/dt |C|^2 = 2 (X (X'W - XW') + Z (Z'W - ZW')) / W^3
      const Homogeneous h = homogeneous(s);
      const auto dX = geometry::bernsteinDerivative(h.X), dZ = geometry::bernsteinDerivative(h.Z),
                 dW = geometry::bernsteinDerivative(h.W);
      using geometry::bernsteinProduct;
      auto g = add(bernsteinProduct(bernsteinProduct(h.X, dX), h.W), bernsteinProduct(bernsteinProduct(h.X, h.X), dW), -1);
      g = add(g, bernsteinProduct(bernsteinProduct(h.Z, dZ), h.W));
      g = add(g, bernsteinProduct(bernsteinProduct(h.Z, h.Z), dW), -1);
      for (const double t : geometry::bernsteinRoots(g)) {
        const Point2 p = s.point(t);
        radii.push_back(std::hypot(p[0], p[1]));
      }
    }
  std::sort(radii.begin(), radii.end());
  const double tol = 1e-12 * std::max(rmax, 1e-300);
  radii.erase(std::unique(radii.begin(), radii.end(), [&](double a, double b) { return b - a <= tol; }),
              radii.end());
  return radii;
}

bool hasRotationSymmetry(const std::vector<Region>& regions, int P, double tol) {
  FG_REQUIRE(P >= 1 && tol > 0, "hasRotationSymmetry: P >= 1, tol > 0");
  if (P == 1) return true;
  const double c = std::cos(2 * PI / P), s = std::sin(2 * PI / P);
  for (const Region& reg : regions)
    for (const Contour2D* con : contours(reg))
      for (const BezierSegment& seg : con->segments)
        for (int k = 0; k <= 16; ++k) {
          const Point2 p = seg.point(k / 16.0), q{c * p[0] - s * p[1], s * p[0] + c * p[1]};
          bool onBoundary = false;
          for (const Region& other : regions) {
            if (other.material != reg.material) continue;
            for (const Contour2D* oc : contours(other)) {
              for (const BezierSegment& os : oc->segments)
                if (distanceToSegment(os, q) <= tol) {
                  onBoundary = true;
                  break;
                }
              if (onBoundary) break;
            }
            if (onBoundary) break;
          }
          if (!onBoundary) return false;
        }
  return true;
}

int rotationOrder(const std::vector<Region>& regions, int maxP, double tol) {
  FG_REQUIRE(maxP >= 1, "rotationOrder: maxP >= 1");
  for (int P = maxP; P > 1; --P)
    if (hasRotationSymmetry(regions, P, tol)) return P;
  return 1;
}

}  // namespace fg::cyl
