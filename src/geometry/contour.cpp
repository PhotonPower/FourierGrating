// FourierGrating - rational Bezier contours and regions (Plan section "M13: Konturdarstellung").
#include "fg/geometry/contour.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"

namespace fg::geometry {

namespace {

using H3 = std::array<double, 3>;  ///< homogeneous (w x, w z, w)

std::vector<H3> homogeneous(const BezierSegment& s) {
  std::vector<H3> h(s.P.size());
  for (std::size_t j = 0; j < s.P.size(); ++j) {
    const double w = s.weight(j);
    h[j] = {w * s.P[j][0], w * s.P[j][1], w};
  }
  return h;
}

BezierSegment fromHomogeneous(const std::vector<H3>& h) {
  BezierSegment s;
  bool unit = true;
  for (const H3& c : h) {
    s.P.push_back({c[0] / c[2], c[1] / c[2]});
    s.w.push_back(c[2]);
    unit = unit && c[2] == 1.0;
  }
  if (unit) s.w.clear();
  return s;
}

H3 casteljau(std::vector<H3> c, double t) {
  for (std::size_t r = c.size() - 1; r > 0; --r)
    for (std::size_t j = 0; j < r; ++j)
      for (int e = 0; e < 3; ++e) c[j][e] = (1 - t) * c[j][e] + t * c[j + 1][e];
  return c[0];
}

double pointChordDistance(const Point2& p, const Point2& a, const Point2& b) {
  const double tx = b[0] - a[0], tz = b[1] - a[1], l2 = tx * tx + tz * tz;
  if (l2 == 0.0) return std::hypot(p[0] - a[0], p[1] - a[1]);
  return std::abs((p[0] - a[0]) * tz - (p[1] - a[1]) * tx) / std::sqrt(l2);
}

void flattenSegment(const BezierSegment& s, double tol, int depth, std::vector<Point2>& out) {
  // convex hull property (positive weights): control polygon distance bounds the sagitta
  double dev = 0;
  for (std::size_t j = 1; j + 1 < s.P.size(); ++j) dev = std::max(dev, pointChordDistance(s.P[j], s.P.front(), s.P.back()));
  if (dev <= tol || depth >= 40) {
    out.push_back(s.P.front());
    return;
  }
  BezierSegment l, r;
  s.split(0.5, l, r);
  flattenSegment(l, tol, depth + 1, out);
  flattenSegment(r, tol, depth + 1, out);
}

bool insidePoly(const std::vector<Point2>& v, double x, double z) {
  bool in = false;
  for (std::size_t j = 0, k = v.size() - 1; j < v.size(); k = j++) {
    const double xj = v[j][0], zj = v[j][1], xk = v[k][0], zk = v[k][1];
    if ((zj > z) != (zk > z) && x < (xk - xj) * (z - zj) / (zk - zj) + xj) in = !in;
  }
  return in;
}

}  // namespace

// ------------------------------------------------------------------ BezierSegment
Point2 BezierSegment::point(double t) const {
  const H3 h = casteljau(homogeneous(*this), t);
  return {h[0] / h[2], h[1] / h[2]};
}

void BezierSegment::eval(double t, Point2& r, Point2& dr) const {
  const std::vector<H3> h = homogeneous(*this);
  const H3 c = casteljau(h, t);
  const int p = degree();
  H3 d{0, 0, 0};
  if (p > 0) {
    std::vector<H3> diff(h.size() - 1);
    for (std::size_t j = 0; j + 1 < h.size(); ++j)
      for (int e = 0; e < 3; ++e) diff[j][e] = p * (h[j + 1][e] - h[j][e]);
    d = casteljau(diff, t);
  }
  r = {c[0] / c[2], c[1] / c[2]};
  dr = {(d[0] * c[2] - c[0] * d[2]) / (c[2] * c[2]), (d[1] * c[2] - c[1] * d[2]) / (c[2] * c[2])};
}

void BezierSegment::split(double t, BezierSegment& left, BezierSegment& right) const {
  std::vector<H3> c = homogeneous(*this), l, r;
  const std::size_t n = c.size();
  l.push_back(c.front());
  r.push_back(c.back());
  for (std::size_t k = n - 1; k > 0; --k) {
    for (std::size_t j = 0; j < k; ++j)
      for (int e = 0; e < 3; ++e) c[j][e] = (1 - t) * c[j][e] + t * c[j + 1][e];
    l.push_back(c.front());
    r.push_back(c[k - 1]);
  }
  std::reverse(r.begin(), r.end());
  left = fromHomogeneous(l);
  right = fromHomogeneous(r);
}

double BezierSegment::controlLength() const {
  double L = 0;
  for (std::size_t j = 1; j < P.size(); ++j) L += std::hypot(P[j][0] - P[j - 1][0], P[j][1] - P[j - 1][1]);
  return L;
}

BezierSegment BezierSegment::reversed() const {
  BezierSegment s = *this;
  std::reverse(s.P.begin(), s.P.end());
  std::reverse(s.w.begin(), s.w.end());
  return s;
}

// ------------------------------------------------------------------ Contour2D
double Contour2D::signedArea() const {
  // Plan eq:contour0: A = 1/2 sum_s int (x z' - z x') dt, 48-point Gauss per half segment
  const auto& g = quadrature::gaussLegendre(48);
  double a = 0;
  for (const BezierSegment& s : segments) {
    BezierSegment halves[2];
    s.split(0.5, halves[0], halves[1]);
    for (const BezierSegment& h : halves)
      for (std::size_t l = 0; l < g.x.size(); ++l) {
        Point2 r, dr;
        h.eval(0.5 * (g.x[l] + 1), r, dr);
        a += 0.25 * g.w[l] * (r[0] * dr[1] - r[1] * dr[0]);
      }
  }
  return a;
}

Contour2D Contour2D::reversed() const {
  Contour2D c;
  for (auto it = segments.rbegin(); it != segments.rend(); ++it) c.segments.push_back(it->reversed());
  return c;
}

std::array<double, 4> Contour2D::bounds() const {
  std::array<double, 4> b{1e300, 1e300, -1e300, -1e300};
  for (const auto& s : segments)
    for (const auto& p : s.P) {
      b[0] = std::min(b[0], p[0]);
      b[1] = std::min(b[1], p[1]);
      b[2] = std::max(b[2], p[0]);
      b[3] = std::max(b[3], p[1]);
    }
  return b;
}

void Contour2D::validate() const {
  FG_REQUIRE(!segments.empty(), "contour: no segments");
  const auto b = bounds();
  const double tol = 1e-12 * std::max({b[2] - b[0], b[3] - b[1], 1e-300});
  for (std::size_t k = 0; k < segments.size(); ++k) {
    const BezierSegment& s = segments[k];
    FG_REQUIRE(s.degree() >= 1 && s.degree() <= 9, "contour: segment degree must be 1..9");
    FG_REQUIRE(s.w.empty() || s.w.size() == s.P.size(), "contour: one weight per control point");
    for (double w : s.w) FG_REQUIRE(w > 0, "contour: weights must be positive");
    const Point2& e = s.P.back();
    const Point2& a = segments[(k + 1) % segments.size()].P.front();
    FG_REQUIRE(std::hypot(e[0] - a[0], e[1] - a[1]) <= tol, "contour: segments do not form a closed chain");
  }
}

std::vector<Point2> Contour2D::flatten(double tol) const {
  std::vector<Point2> out;
  for (const auto& s : segments) flattenSegment(s, tol, 0, out);
  return out;
}

Contour2D Contour2D::polygon(const std::vector<Point2>& v) {
  FG_REQUIRE(v.size() >= 3, "contour: polygon needs at least 3 vertices");
  Contour2D c;
  for (std::size_t j = 0; j < v.size(); ++j) c.segments.push_back({{v[j], v[(j + 1) % v.size()]}, {}});
  return c;
}

Contour2D Contour2D::circle(Point2 center, double radius) { return ellipse(center, radius, radius, 0.0); }

Contour2D Contour2D::ellipse(Point2 c, double a, double b, double angle) {
  FG_REQUIRE(a > 0 && b > 0, "contour: radii must be positive");
  const double ca = std::cos(angle), sa = std::sin(angle), h = std::sqrt(0.5);
  auto map = [&](double u, double v) { return Point2{c[0] + a * ca * u - b * sa * v, c[1] + a * sa * u + b * ca * v}; };
  // quarter arcs (1,0)-(1,1)-(0,1) etc. of the unit circle, exact as rational quadratics
  const double U[5] = {1, 0, -1, 0, 1}, V[5] = {0, 1, 0, -1, 0};
  Contour2D out;
  for (int k = 0; k < 4; ++k)
    out.segments.push_back(
        {{map(U[k], V[k]), map(U[k] + U[k + 1], V[k] + V[k + 1]), map(U[k + 1], V[k + 1])}, {1.0, h, 1.0}});
  return out;
}

Contour2D Contour2D::nurbs(int p, const std::vector<double>& U, const std::vector<Point2>& pts,
                           const std::vector<double>& wts) {
  // Piegl-Tiller A5.6 (DecomposeCurve) in homogeneous coordinates
  const int n = static_cast<int>(pts.size()) - 1, m = n + p + 1;
  FG_REQUIRE(p >= 1 && p <= 9, "nurbs: degree must be 1..9");
  FG_REQUIRE(n >= p, "nurbs: need at least degree + 1 control points");
  FG_REQUIRE(static_cast<int>(U.size()) == m + 1, "nurbs: knot vector size must be points + degree + 1");
  FG_REQUIRE(wts.empty() || wts.size() == pts.size(), "nurbs: one weight per control point");
  for (int j = 1; j <= m; ++j) FG_REQUIRE(U[static_cast<std::size_t>(j)] >= U[static_cast<std::size_t>(j - 1)], "nurbs: knots must be non-decreasing");
  for (int j = 1; j <= p; ++j)
    FG_REQUIRE(U[static_cast<std::size_t>(j)] == U[0] && U[static_cast<std::size_t>(m - j)] == U[static_cast<std::size_t>(m)],
               "nurbs: knot vector must be clamped (first and last knot with multiplicity degree + 1)");
  std::vector<H3> Pw(pts.size());
  for (std::size_t j = 0; j < pts.size(); ++j) {
    const double w = wts.empty() ? 1.0 : wts[j];
    FG_REQUIRE(w > 0, "nurbs: weights must be positive");
    Pw[j] = {w * pts[j][0], w * pts[j][1], w};
  }
  auto u = [&](int j) { return U[static_cast<std::size_t>(j)]; };
  std::vector<std::vector<H3>> Q(1, std::vector<H3>(static_cast<std::size_t>(p + 1)));
  int a = p, b = p + 1, nb = 0;
  for (int i = 0; i <= p; ++i) Q[0][static_cast<std::size_t>(i)] = Pw[static_cast<std::size_t>(i)];
  std::vector<double> alphas(static_cast<std::size_t>(p));
  while (b < m) {
    const int i0 = b;
    while (b < m && u(b + 1) == u(b)) ++b;
    const int mult = b - i0 + 1;
    Q.resize(static_cast<std::size_t>(nb + 2), std::vector<H3>(static_cast<std::size_t>(p + 1)));
    if (mult < p) {
      const double numer = u(b) - u(a);
      for (int j = p; j > mult; --j) alphas[static_cast<std::size_t>(j - mult - 1)] = numer / (u(a + j) - u(a));
      const int r = p - mult;
      for (int j = 1; j <= r; ++j) {
        const int save = r - j, s = mult + j;
        auto& q = Q[static_cast<std::size_t>(nb)];
        for (int k = p; k >= s; --k) {
          const double al = alphas[static_cast<std::size_t>(k - s)];
          for (int e = 0; e < 3; ++e)
            q[static_cast<std::size_t>(k)][e] = al * q[static_cast<std::size_t>(k)][e] + (1 - al) * q[static_cast<std::size_t>(k - 1)][e];
        }
        if (b < m) Q[static_cast<std::size_t>(nb + 1)][static_cast<std::size_t>(save)] = q[static_cast<std::size_t>(p)];
      }
    }
    ++nb;
    if (b < m) {
      for (int i = p - mult; i <= p; ++i)
        Q[static_cast<std::size_t>(nb)][static_cast<std::size_t>(i)] = Pw[static_cast<std::size_t>(b - p + i)];
      a = b;
      ++b;
    }
  }
  Contour2D c;
  for (int k = 0; k < nb; ++k) c.segments.push_back(fromHomogeneous(Q[static_cast<std::size_t>(k)]));
  return c;
}

// ------------------------------------------------------------------ Region
void Region::normalize(double flattenTol) {
  outer.validate();
  if (outer.signedArea() < 0) outer = outer.reversed();
  for (auto& h : holes) {
    h.validate();
    if (h.signedArea() > 0) h = h.reversed();
  }
  poly_.clear();
  poly_.push_back(outer.flatten(flattenTol));
  for (const auto& h : holes) poly_.push_back(h.flatten(flattenTol));
}

bool Region::contains(double x, double z) const {
  FG_REQUIRE(!poly_.empty(), "region: normalize() not called");
  if (!insidePoly(poly_[0], x, z)) return false;
  for (std::size_t k = 1; k < poly_.size(); ++k)
    if (insidePoly(poly_[k], x, z)) return false;
  return true;
}

double Region::area() const {
  double a = outer.signedArea();
  for (const auto& h : holes) a += h.signedArea();
  return a;
}

}  // namespace fg::geometry
