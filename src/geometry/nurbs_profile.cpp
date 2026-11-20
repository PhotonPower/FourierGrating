// FourierGrating - open periodic NURBS profile of a 1D grating (Plan section "Erweiterung ab Version 2.37:
// NURBS-Profil in 1D-Gittern (M24)").
#include "fg/geometry/nurbs_profile.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "fg/core/error.hpp"

namespace fg {
namespace geometry {

namespace {

double binom(int n, int k) {
  double r = 1;
  for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
  return r;
}

/// Value of a Bernstein polynomial at t (de Casteljau).
double bernstein(std::vector<double> c, double t) {
  for (std::size_t r = 1; r < c.size(); ++r)
    for (std::size_t j = 0; j + r < c.size(); ++j) c[j] = (1 - t) * c[j] + t * c[j + 1];
  return c.front();
}

/// Coefficients on [0, 1/2] and [1/2, 1] of a Bernstein polynomial (de Casteljau at 1/2).
void split(const std::vector<double>& c, std::vector<double>& l, std::vector<double>& r) {
  std::vector<double> w = c;
  const std::size_t n = c.size();
  l.assign(n, 0.0);
  r.assign(n, 0.0);
  l[0] = w[0];
  r[n - 1] = w[n - 1];
  for (std::size_t k = 1; k < n; ++k) {
    for (std::size_t j = 0; j + k < n; ++j) w[j] = 0.5 * (w[j] + w[j + 1]);
    l[k] = w[0];
    r[n - 1 - k] = w[n - 1 - k];
  }
}

int signChanges(const std::vector<double>& c) {
  int v = 0, last = 0;
  for (double x : c) {
    const int s = x > 0 ? 1 : x < 0 ? -1 : 0;
    if (s == 0) continue;
    if (last != 0 && s != last) ++v;
    last = s;
  }
  return v;
}

/// Simple roots of the Bernstein polynomial `full` inside (a, b), whose coefficients on [a, b] are c: subdivision
/// until one sign change (exactly one root, variation diminishing), then bisection on the full polynomial.
void isolate(const std::vector<double>& full, const std::vector<double>& c, double a, double b, int depth,
             std::vector<double>& roots) {
  const int v = signChanges(c);
  if (v == 0) return;
  if (v == 1 && c.front() != 0 && c.back() != 0) {
    double lo = a, hi = b, flo = c.front();
    for (int it = 0; it < 200 && hi - lo > 1e-17; ++it) {
      const double m = 0.5 * (lo + hi), fm = bernstein(full, m);
      if (fm == 0) {
        lo = hi = m;
        break;
      }
      if ((fm > 0) == (flo > 0)) {
        lo = m;
        flo = fm;
      } else {
        hi = m;
      }
    }
    roots.push_back(0.5 * (lo + hi));
    return;
  }
  if (depth > 60) {  // a multiple root (no sign change of y'), not an extremum we can separate
    return;
  }
  std::vector<double> l, r;
  split(c, l, r);
  const double m = 0.5 * (a + b);
  isolate(full, l, a, m, depth + 1, roots);
  if (l.back() == 0) roots.push_back(m);  // root exactly at the split point
  isolate(full, r, m, b, depth + 1, roots);
}

/// Message of an error without the source location FG_REQUIRE puts in front ("file.cpp:123: ").
std::string withoutLocation(const std::string& w) {
  const auto p = w.find(".cpp:");
  if (p == std::string::npos) return w;
  const auto q = w.find(": ", p + 5);
  return q == std::string::npos ? w : w.substr(q + 2);
}

std::string fmt(double v) {
  std::ostringstream o;
  o.precision(10);
  o << v;
  return o.str();
}

}  // namespace

std::vector<double> yExtremaParameters(const BezierSegment& s) {
  // Plan eq:nurbs1d_crit: N'D - ND' as a Bernstein polynomial of degree 2p - 1
  const int p = s.degree();
  std::vector<double> N(static_cast<std::size_t>(p + 1)), D(static_cast<std::size_t>(p + 1));
  for (int j = 0; j <= p; ++j) {
    D[static_cast<std::size_t>(j)] = s.weight(static_cast<std::size_t>(j));
    N[static_cast<std::size_t>(j)] = D[static_cast<std::size_t>(j)] * s.P[static_cast<std::size_t>(j)][1];
  }
  const int n = 2 * p - 1;
  std::vector<double> c(static_cast<std::size_t>(n + 1), 0.0);
  double scale = 0;
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= p; ++j) {
      const auto ui = static_cast<std::size_t>(i), uj = static_cast<std::size_t>(j);
      const double f = binom(p - 1, i) * binom(p, j) / binom(n, i + j) * p;
      const double dN = N[ui + 1] - N[ui], dD = D[ui + 1] - D[ui];
      c[static_cast<std::size_t>(i + j)] += f * (dN * D[uj] - dD * N[uj]);
      scale = std::max(scale, std::abs(f * dN * D[uj]) + std::abs(f * dD * N[uj]));
    }
  // constant y (horizontal piece): the numerator vanishes up to rounding
  double mx = 0;
  for (double x : c) mx = std::max(mx, std::abs(x));
  if (mx <= 1e-13 * scale || mx == 0) return {};
  std::vector<double> t;
  isolate(c, c, 0.0, 1.0, 0, t);
  std::vector<double> out;
  for (double x : t)
    if (x > 1e-12 && x < 1 - 1e-12) out.push_back(x);
  std::sort(out.begin(), out.end());
  return out;
}

NurbsProfile::NurbsProfile(const NurbsCurve& c, double dx_um, std::string below, std::string above)
    : dx_(dx_um), below_(std::move(below)), above_(std::move(above)) {
  if (!(dx_ > 0)) throw Error("[NURBS_INPUT] nurbs: the period must be positive");
  try {
    FG_REQUIRE(c.points.size() >= 2, "nurbs: at least two control points");
    pieces_ = Contour2D::nurbs(c.degree, c.knots, c.points, c.weights).segments;
  } catch (const Error& e) {
    throw Error("[NURBS_INPUT] " + withoutLocation(e.what()));
  }
  if (pieces_.empty()) throw Error("[NURBS_INPUT] nurbs: the curve has no non-empty knot span");
  const Point2 a = pieces_.front().P.front(), b = pieces_.back().P.back();
  double ext = dx_;
  for (const auto& s : pieces_)
    for (const auto& q : s.P) ext = std::max({ext, std::abs(q[0] - a[0]), std::abs(q[1] - a[1])});
  const double tol = 1e-12 * ext;
  if (std::abs(b[1] - a[1]) > tol)
    throw Error("[NURBS_NOT_PERIODIC] nurbs: the curve ends at y = " + fmt(b[1]) + " um but starts at y = " + fmt(a[1]) +
                " um (end = start + (period, 0) required; no wall is inserted)");
  if (std::abs(b[0] - a[0] - dx_) > tol)
    throw Error("[NURBS_PERIOD] nurbs: the curve runs from x = " + fmt(a[0]) + " to x = " + fmt(b[0]) +
                " um, not over one period " + fmt(dx_) + " um in +x" +
                (std::abs(b[0] - a[0] + dx_) <= tol ? " (it runs in -x: reverse the control points)" : ""));
  pieces_.back().P.back() = {a[0] + dx_, a[1]};  // the seam exactly
  // range from the piece ends and the interior extrema (exact, not the control points)
  ymin_ = ymax_ = a[1];
  for (const auto& s : pieces_) {
    for (const double y : {s.P.front()[1], s.P.back()[1]}) {
      ymin_ = std::min(ymin_, y);
      ymax_ = std::max(ymax_, y);
    }
    for (double t : yExtremaParameters(s)) {
      const double y = s.point(t)[1];
      ymin_ = std::min(ymin_, y);
      ymax_ = std::max(ymax_, y);
    }
  }
  // self-intersection: the flattened curve against itself and its images at +-dx (Plan "M24: Prüfungen")
  Contour2D open;
  open.segments = pieces_;
  std::vector<Point2> poly = open.flatten(1e-7 * dx_);
  if (poly.empty() || std::hypot(poly.back()[0] - b[0], poly.back()[1] - b[1]) > tol) poly.push_back(pieces_.back().P.back());
  const double ptol = 1e-9 * dx_;
  for (const double shift : {0.0, dx_, -dx_})
    for (std::size_t i = 0; i + 1 < poly.size(); ++i)
      for (std::size_t j = (shift == 0 ? i + 2 : 0); j + 1 < poly.size(); ++j) {
        const Point2 c0{poly[j][0] + shift, poly[j][1]}, c1{poly[j + 1][0] + shift, poly[j + 1][1]};
        if (properCrossing(poly[i], poly[i + 1], c0, c1, ptol))
          throw Error("[NURBS_SELF_INTERSECTION] nurbs: the curve crosses " +
                      std::string(shift == 0 ? "itself" : "its periodic image") + " near (" + fmt(poly[i][0]) + ", " +
                      fmt(poly[i][1]) + ") um");
      }
}

std::vector<double> NurbsProfile::criticalHeights(double tol) const {
  std::vector<double> h;
  for (const auto& s : pieces_) {
    h.push_back(s.P.front()[1]);
    h.push_back(s.P.back()[1]);
    for (double t : yExtremaParameters(s)) h.push_back(s.point(t)[1]);
  }
  std::sort(h.begin(), h.end());
  std::vector<double> out;
  for (double y : h)
    if (out.empty() || y > out.back() + tol) out.push_back(y);
  return out;
}

LineSection NurbsProfile::section(double y) const {
  std::vector<Crossing> cr;
  for (const auto& s : pieces_) lineCrossings(s, y, cr);  // t in [0, 1): the seam once
  LineSection out;
  if (cr.empty()) {  // under or over the whole curve
    out.breaks = {0.0, dx_};
    out.materials = {y < ymin_ ? below_ : above_};
    out.normals = {Vec3D(1, 0, 0)};
    return out;
  }
  // crossings within 1e-10 dx of the seam (t = 0 of the first piece, t -> 1 of the last within rounding) lie on it
  const double tol = 1e-10 * dx_;
  struct K {
    double x;
    Vec3D n;
    bool below;  // material right of the crossing (Plan eq:nurbs1d_rule)
  };
  std::vector<K> ks;
  for (const auto& c : cr) {
    double x = c.x - dx_ * std::floor(c.x / dx_);
    if (x >= dx_ - tol || x < 0) x = 0.0;  // the seam
    ks.push_back({x, c.n, c.n(0) > 0});
  }
  std::sort(ks.begin(), ks.end(), [](const K& p, const K& q) { return p.x < q.x; });
  // crossings closer than tol: the same direction is one crossing found twice (the seam: t = 0 of the first piece
  // and t -> 1 of the last), opposite directions bound a sliver narrower than tol next to an extremum (e.g. a valley
  // at the layer bottom, sampled 1e-14 above it), which is dropped
  {
    std::vector<K> m;
    for (const K& k : ks) {
      if (!m.empty() && k.x - m.back().x <= tol) {
        if (k.below != m.back().below) m.pop_back();
        continue;
      }
      m.push_back(k);
    }
    ks = std::move(m);
  }
  if (ks.empty()) {
    out.breaks = {0.0, dx_};
    out.materials = {y < 0.5 * (ymin_ + ymax_) ? below_ : above_};
    out.normals = {Vec3D(1, 0, 0)};
    return out;
  }
  for (std::size_t k = 0; k < ks.size(); ++k)
    if (ks[k].below == ks[(k + 1) % ks.size()].below || (k + 1 < ks.size() && ks[k + 1].x - ks[k].x <= tol))
      throw Error("[NURBS_SELF_INTERSECTION] nurbs: at y = " + fmt(y) +
                  " um the materials below/above do not alternate along x (self-intersecting curve)");
  const std::size_t K0 = ks.size();
  const bool atZero = ks.front().x <= tol;
  out.breaks.push_back(0.0);
  if (atZero) {
    out.normals.push_back(ks.front().n);
    out.materials.push_back(ks.front().below ? below_ : above_);
  } else {
    out.normals.push_back(Vec3D(1, 0, 0));  // no boundary at 0: the material of the last crossing continues
    out.materials.push_back(ks.back().below ? below_ : above_);
  }
  for (std::size_t k = atZero ? 1 : 0; k < K0; ++k) {
    out.breaks.push_back(ks[k].x);
    out.normals.push_back(ks[k].n);
    out.materials.push_back(ks[k].below ? below_ : above_);
  }
  out.breaks.push_back(dx_);
  return out;
}

}  // namespace geometry

std::vector<double> NurbsProfileLineProvider::bandBreaks() const {
  const double tol = 1e-12 * h_;
  std::vector<double> out;
  for (double y : p_->criticalHeights(tol))
    if (y > tol && y < h_ - tol) out.push_back(y);
  return out;
}

}  // namespace fg
