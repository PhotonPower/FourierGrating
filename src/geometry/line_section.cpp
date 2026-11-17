// FourierGrating - line sections of 1D profiles (Plan section "M23").
#include "fg/geometry/line_section.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

#include "fg/core/error.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/geometry/coordinate_transform.hpp"
#include "fg/geometry/normal_field.hpp"

namespace fg {
namespace geometry {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Point (x, y) and derivative of the segment at t.
void evalXY(const BezierSegment& s, double t, double& x, double& y, double& dx, double& dy) {
  Point2 r, d;
  s.eval(t, r, d);
  x = r[0];
  y = r[1];
  dx = d[0];
  dy = d[1];
}

/// Root of y(t) - c in [a, b] with a sign change: safeguarded Newton (bisection fallback).
double refine(const BezierSegment& s, double c, double a, double b) {
  double x, ya, yb, dx, dy;
  evalXY(s, a, x, ya, dx, dy);
  evalXY(s, b, x, yb, dx, dy);
  ya -= c;
  yb -= c;
  if (ya == 0) return a;
  if (yb == 0) return b;
  double lo = a, hi = b;  // g(lo) has the sign of ya
  double t = 0.5 * (a + b);
  for (int it = 0; it < 200; ++it) {
    double y;
    evalXY(s, t, x, y, dx, dy);
    const double g = y - c;
    if (g == 0) return t;
    if ((g < 0) == (ya < 0))
      lo = t;
    else
      hi = t;
    double tn = dy != 0 ? t - g / dy : 0.5 * (lo + hi);
    if (!(tn > std::min(lo, hi) && tn < std::max(lo, hi))) tn = 0.5 * (lo + hi);
    if (std::abs(tn - t) <= 1e-16 * std::max(1.0, std::abs(t)) || std::abs(hi - lo) <= 4e-16) return tn;
    t = tn;
  }
  return t;
}

int sgn(double v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

/// Roots of the Bernstein polynomial with coefficients c on [t0, t1] (of the whole segment).
void isolate(const BezierSegment& s, double level, std::vector<double> c, double t0, double t1, int depth,
             std::vector<double>& roots) {
  const std::size_t n = c.size() - 1;
  if (c.front() == 0) {  // a root at the start belongs to this interval (t = 1 to the next segment)
    if (t0 < 1) roots.push_back(t0);
    std::size_t k = 1;
    while (k <= n && c[k] == 0) ++k;
    if (k > n) return;  // identically zero: the segment lies on the line (critical height)
    c.front() = c[k];
  }
  if (c.back() == 0) {
    std::size_t k = n;
    while (k > 0 && c[k - 1] == 0) --k;
    c.back() = c[k - 1];
  }
  int changes = 0, last = 0;
  for (double v : c) {
    const int sv = sgn(v);
    if (sv == 0) continue;
    if (last != 0 && sv != last) ++changes;
    last = sv;
  }
  if (changes == 0) return;
  if (changes == 1) {  // exactly one root (Descartes' rule in Bernstein form), ends of opposite sign
    roots.push_back(refine(s, level, t0, t1));
    return;
  }
  if (depth > 60 || t1 - t0 < 1e-15) {
    roots.push_back(0.5 * (t0 + t1));
    return;
  }
  // de Casteljau at 1/2
  std::vector<double> L(n + 1), R(n + 1), b = c;
  for (std::size_t r = 0; r <= n; ++r) {
    L[r] = b[0];
    R[n - r] = b[n - r];
    for (std::size_t j = 0; j + r < n; ++j) b[j] = 0.5 * (b[j] + b[j + 1]);
  }
  const double tm = 0.5 * (t0 + t1);
  isolate(s, level, L, t0, tm, depth + 1, roots);
  isolate(s, level, R, tm, t1, depth + 1, roots);
}

struct Piece {
  double a = 0, b = 0;
  bool hasA = false, hasB = false;  ///< false: seam (end of a periodic image, no interface)
  Vec3D nA = Vec3D(1, 0, 0), nB = Vec3D(1, 0, 0);
  int body = -1;
};

/// Intervals of one body folded into [0, dx): parity rule, split at the cell boundary, merge of
/// touching images (seams), overlaps of images united (and reported).
std::vector<Piece> bodyPieces(const ProfileBody& body, int index, double y, double dx, double tol, bool unite,
                              LineSectionFindings* findings) {
  std::vector<Crossing> cr;
  for (const auto& loop : body.loops)
    for (const auto& seg : loop.segments) lineCrossings(seg, y, cr);
  FG_REQUIRE(cr.size() % 2 == 0, "line section: body '" + body.name + "' has an odd number of crossings at y = " +
                                     std::to_string(y) + " um (open profile or a tangent outside a critical height)");
  std::sort(cr.begin(), cr.end(), [](const Crossing& p, const Crossing& q) { return p.x < q.x; });
  std::vector<Piece> raw;
  bool overlapping = false;
  for (std::size_t k = 0; k + 1 < cr.size(); k += 2) {
    const double a = cr[k].x, b = cr[k + 1].x;
    if (b - a <= tol) continue;
    if (b - a > dx + tol) overlapping = true;
    if (b - a >= dx - tol) {  // covers the period: a seam onto itself
      raw.push_back({0.0, dx, false, false, Vec3D(1, 0, 0), Vec3D(1, 0, 0), index});
      continue;
    }
    const double s = std::floor(a / dx) * dx, a0 = a - s, b0 = b - s;
    if (b0 <= dx + tol) {
      raw.push_back({a0, std::min(b0, dx), true, true, cr[k].n, cr[k + 1].n, index});
    } else {
      raw.push_back({a0, dx, true, false, cr[k].n, Vec3D(1, 0, 0), index});
      raw.push_back({0.0, b0 - dx, false, true, Vec3D(1, 0, 0), cr[k + 1].n, index});
    }
  }
  std::sort(raw.begin(), raw.end(), [](const Piece& p, const Piece& q) { return p.a < q.a; });
  std::vector<Piece> out;
  for (const Piece& p : raw) {
    if (out.empty() || p.a > out.back().b + tol) {
      out.push_back(p);
      continue;
    }
    Piece& c = out.back();
    if (p.a < c.b - tol) overlapping = true;  // an image overlaps: periodic union or finding
    if (p.b > c.b) {
      c.b = p.b;
      c.hasB = p.hasB;
      c.nB = p.nB;
    }
    // touching (|p.a - c.b| <= tol): merged, the shared point is a seam
  }
  // across the cell boundary: [.., dx] and [0, ..] of the same body join without an interface
  if (out.size() > 1 && out.front().a <= tol && out.back().b >= dx - tol) {
    out.front().hasA = false;
    out.back().hasB = false;
  }
  if (out.size() == 1 && out.front().a <= tol && out.front().b >= dx - tol) {
    out.front() = {0.0, dx, false, false, Vec3D(1, 0, 0), Vec3D(1, 0, 0), index};
  }
  if (overlapping) {
    if (findings)
      findings->selfOverlap.push_back(index);
    else
      FG_REQUIRE(unite, "line section: body '" + body.name +
                            "' overlaps its periodic image (cad.periodic_union unites them)");
  }
  return out;
}

/// p minus [a, b] (the cut ends take the normals of the cutting piece).
std::vector<Piece> subtract(const Piece& p, const Piece& q) {
  if (q.b <= p.a || q.a >= p.b) return {p};
  std::vector<Piece> out;
  if (q.a > p.a) {
    Piece l = p;
    l.b = q.a;
    l.hasB = q.hasA;
    l.nB = q.nA;
    out.push_back(l);
  }
  if (q.b < p.b) {
    Piece r = p;
    r.a = q.b;
    r.hasA = q.hasB;
    r.nA = q.nB;
    out.push_back(r);
  }
  return out;
}

}  // namespace

void lineCrossings(const BezierSegment& s, double c, std::vector<Crossing>& out) {
  const std::size_t n = s.P.size();
  FG_REQUIRE(n >= 2, "line section: Bezier segment with fewer than two control points");
  std::vector<double> coef(n);
  bool pos = false, neg = false;
  for (std::size_t j = 0; j < n; ++j) {
    coef[j] = s.weight(j) * (s.P[j][1] - c);  // weights > 0: same sign as y(t) - c
    pos = pos || coef[j] > 0;
    neg = neg || coef[j] < 0;
  }
  if (!(pos && neg) && coef.front() != 0) return;  // entirely on one side
  std::vector<double> roots;
  isolate(s, c, coef, 0.0, 1.0, 0, roots);
  for (double t : roots) {
    double x, y, dx, dy;
    evalXY(s, t, x, y, dx, dy);
    const double l = std::hypot(dx, dy);
    Crossing k;
    k.x = x;
    if (l > 0) k.n = Vec3D(dy / l, -dx / l, 0.0);  // Plan eq:cad1d_n
    out.push_back(k);
  }
}

LineSection lineSection(const std::vector<ProfileBody>& bodies, double y, double dx, const std::string& background,
                        const LineSectionOptions& opt, LineSectionFindings* findings) {
  FG_REQUIRE(dx > 0, "line section: period must be positive");
  const double tol = opt.tol * dx;
  // bodies one after the other; earlier bodies win overlaps (priority)
  std::vector<Piece> all;
  for (std::size_t i = 0; i < bodies.size(); ++i) {
    std::vector<Piece> mine =
        bodyPieces(bodies[i], static_cast<int>(i), y, dx, tol, opt.periodicUnion, findings);
    std::vector<int> hit;
    for (const Piece& q : all) {
      std::vector<Piece> next;
      for (const Piece& p : mine) {
        if (std::min(p.b, q.b) - std::max(p.a, q.a) > tol &&
            std::find(hit.begin(), hit.end(), q.body) == hit.end())
          hit.push_back(q.body);
        for (Piece& r : subtract(p, q))
          if (r.b - r.a > tol) next.push_back(r);
      }
      mine = std::move(next);
    }
    for (int j : hit) {
      if (findings) {
        findings->bodyOverlap.push_back({static_cast<int>(i), j});
      } else {
        FG_REQUIRE(opt.priority, "line section: bodies '" + bodies[i].name + "' and '" +
                                     bodies[static_cast<std::size_t>(j)].name +
                                     "' overlap (cad.overlap = \"priority\" resolves it)");
      }
    }
    all.insert(all.end(), mine.begin(), mine.end());
  }
  std::sort(all.begin(), all.end(), [](const Piece& p, const Piece& q) { return p.a < q.a; });
  // cover [0, dx]: pieces and background gaps
  struct Seg {
    double a, b;
    std::string mat;
    bool hasA, hasB;
    Vec3D nA, nB;
  };
  std::vector<Seg> seg;
  double x = 0;
  for (const Piece& p : all) {
    const bool gap = p.a > x + tol;  // otherwise touching: snapped to the previous end
    if (gap) seg.push_back({x, p.a, background, false, false, Vec3D(1, 0, 0), Vec3D(1, 0, 0)});
    seg.push_back({gap ? p.a : x, p.b, bodies[static_cast<std::size_t>(p.body)].material, p.hasA, p.hasB, p.nA,
                   p.nB});
    x = p.b;
  }
  if (x < dx - tol) seg.push_back({x, dx, background, false, false, Vec3D(1, 0, 0), Vec3D(1, 0, 0)});
  if (seg.empty()) seg.push_back({0.0, dx, background, false, false, Vec3D(1, 0, 0), Vec3D(1, 0, 0)});
  seg.back().b = dx;
  seg.front().a = 0;
  // breaks with normals; boundaries without material change are removed
  const std::size_t S = seg.size();
  LineSection out;
  for (std::size_t k = 0; k < S; ++k) {
    const Seg& prev = seg[(k + S - 1) % S];
    const Seg& cur = seg[k];
    if (k > 0 && prev.mat == cur.mat) continue;
    out.breaks.push_back(cur.a);
    out.materials.push_back(cur.mat);
    out.normals.push_back(prev.hasB ? prev.nB : (cur.hasA ? cur.nA : Vec3D(1, 0, 0)));
  }
  out.breaks.push_back(dx);
  return out;
}

std::vector<std::size_t> materialBoundaries(const LineSection& s) {
  std::vector<std::size_t> b;
  const std::size_t K = s.materials.size();
  for (std::size_t k = 0; k < K; ++k)
    if (s.materials[k] != s.materials[(k + K - 1) % K]) b.push_back(k);
  return b;
}

double normalAngle(const LineSection& s, double dx, double x) {
  const auto b = materialBoundaries(s);
  if (b.empty()) return 0.0;
  const std::size_t M = b.size();
  auto theta = [&](std::size_t m) {
    const Vec3D& n = s.normals[b[m]];
    double t = std::atan2(n(1), n(0));
    t = std::fmod(t, kPi);
    if (t < 0) t += kPi;
    return t;
  };
  const double xr = x - std::floor(x / dx) * dx;
  // last boundary <= xr (cyclic: the last one, shifted by -dx)
  std::size_t m = M - 1;
  double xi = s.breaks[b[M - 1]] - dx;
  for (std::size_t k = 0; k < M; ++k)
    if (s.breaks[b[k]] <= xr) {
      m = k;
      xi = s.breaks[b[k]];
    }
  const std::size_t m1 = (m + 1) % M;
  double xn = s.breaks[b[m1]];
  if (xn <= xi) xn += dx;
  double d = theta(m1) - theta(m);
  while (d > 0.5 * kPi) d -= kPi;
  while (d <= -0.5 * kPi) d += kPi;
  return theta(m) + (xr - xi) / (xn - xi) * d;  // Plan eq:cad1d_N
}

}  // namespace geometry

NormalFieldCoefficients normalField(const LineSection& s, double dx_um, const IndexMap& map,
                                    const CoordinateMap* xmap) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  NormalFieldCoefficients nf;
  for (auto& row : nf.NN)
    for (auto& g : row) g = fourier::CoefficientGrid(Px, Pz);
  for (auto& g : nf.N) g = fourier::CoefficientGrid(Px, Pz);
  nf.hasSigned = false;  // projector form (Plan section "M23: Annahmen")
  // all boundaries parallel (e.g. vertical walls): N is constant, exact coefficients (Li rule)
  const auto bnd = geometry::materialBoundaries(s);
  bool constant = true;
  for (std::size_t k : bnd) {
    const Vec3D& a = s.normals[bnd.front()];
    const Vec3D& b = s.normals[k];
    constant = constant && std::abs(a(0) * b(1) - a(1) * b(0)) <= 1e-15;
  }
  if (constant) {
    const double t = bnd.empty() ? 0.0 : geometry::normalAngle(s, dx_um, s.breaks[bnd.front()]);
    const double cx = bnd.empty() ? 1.0 : std::cos(t), sy = bnd.empty() ? 0.0 : std::sin(t);
    nf.NN[0][0](0, 0) = cx * cx;
    nf.NN[0][1](0, 0) = cx * sy;
    nf.NN[1][0](0, 0) = cx * sy;
    nf.NN[1][1](0, 0) = sy * sy;
    nf.N[0](0, 0) = cx;
    nf.N[1](0, 0) = sy;
    return nf;
  }
  // quadrature breaks: the boundaries (and the breakpoints of the map) in the coordinate X
  std::vector<double> br{0.0, dx_um};
  for (std::size_t k : geometry::materialBoundaries(s)) br.push_back(xmap ? xmap->inverse(s.breaks[k]) : s.breaks[k]);
  if (xmap)
    for (double x : xmap->breakpoints()) br.push_back(x);
  for (double& x : br) x = std::clamp(x - std::floor(x / dx_um) * dx_um, 0.0, dx_um);
  br.push_back(0.0);
  std::sort(br.begin(), br.end());
  br.erase(std::unique(br.begin(), br.end(), [&](double a, double b) { return b - a <= 1e-14 * dx_um; }), br.end());
  if (br.back() < dx_um) br.push_back(dx_um);
  br.back() = dx_um;
  const auto c = fourier::piecewiseSmoothCoefficients(
      br, 5,
      [&](double X, int, cplx* out) {
        const double t = geometry::normalAngle(s, dx_um, xmap ? xmap->f(X) : X);
        const double cx = std::cos(t), sy = std::sin(t);
        out[0] = cx * cx;
        out[1] = cx * sy;
        out[2] = sy * sy;
        out[3] = cx;
        out[4] = sy;
      },
      dx_um, Px);
  for (int p = -Px; p <= Px; ++p) {
    const auto i = static_cast<std::size_t>(p + Px);
    nf.NN[0][0](p, 0) = c[0][i];
    nf.NN[0][1](p, 0) = c[1][i];
    nf.NN[1][0](p, 0) = c[1][i];
    nf.NN[1][1](p, 0) = c[2][i];
    nf.N[0](p, 0) = c[3][i];
    nf.N[1](p, 0) = c[4][i];
  }
  return nf;
}

ProfileLineProvider::ProfileLineProvider(std::shared_ptr<const std::vector<geometry::ProfileBody>> bodies, double y0,
                                         double h, double dx, std::string background,
                                         geometry::LineSectionOptions opt)
    : bodies_(std::move(bodies)), y0_(y0), h_(h), dx_(dx), bg_(std::move(background)), opt_(opt),
      margin_(1e-8 * std::max(dx, h)) {
  FG_REQUIRE(h_ > 0 && dx_ > 0, "profile line provider: thickness and period must be positive");
}

LineSection ProfileLineProvider::section(double y_um) const {
  return geometry::lineSection(*bodies_, y0_ + std::clamp(y_um, margin_, h_ - margin_), dx_, bg_, opt_);
}

}  // namespace fg
