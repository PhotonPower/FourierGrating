// V25a (Plan section "M24: Tests und Abnahme"): open periodic NURBS profile of a 1D grating - Bezier pieces against
// de Boor, critical heights (Plan eq:nurbs1d_crit) against closed forms, line sections with the material rule (Plan
// eq:nurbs1d_rule) against lineSection of an equivalent closed body, the seam, the checks NURBS_*.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numbers>

#include "fg/core/error.hpp"
#include "fg/geometry/nurbs_profile.hpp"

using namespace fg;
using geometry::BezierSegment;
using geometry::Contour2D;
using geometry::NurbsCurve;
using geometry::NurbsProfile;
using geometry::Point2;
using geometry::ProfileBody;

namespace {

constexpr double pi = std::numbers::pi;

/// Point of a NURBS curve at u (de Boor in homogeneous coordinates, Piegl-Tiller A4.1 via A2.1/A3.1).
Point2 deBoor(const NurbsCurve& c, double u) {
  const int p = c.degree, n = static_cast<int>(c.points.size()) - 1;
  const auto& U = c.knots;
  int k = p;  // span with U[k] <= u < U[k + 1] (u = end: last non-empty span)
  for (int i = p; i <= n; ++i)
    if (U[i] <= u && (u < U[i + 1] || (i == n && u == U[i + 1]))) {
      if (U[i] < U[i + 1]) k = i;
    }
  std::vector<std::array<double, 3>> d(static_cast<std::size_t>(p + 1));
  for (int j = 0; j <= p; ++j) {
    const auto i = static_cast<std::size_t>(j + k - p);
    const double w = c.weights.empty() ? 1.0 : c.weights[i];
    d[static_cast<std::size_t>(j)] = {w * c.points[i][0], w * c.points[i][1], w};
  }
  for (int r = 1; r <= p; ++r)
    for (int j = p; j >= r; --j) {
      const double a = (u - U[j + k - p]) / (U[j + 1 + k - r] - U[j + k - p]);
      for (int e = 0; e < 3; ++e)
        d[static_cast<std::size_t>(j)][e] = (1 - a) * d[static_cast<std::size_t>(j - 1)][e] + a * d[static_cast<std::size_t>(j)][e];
    }
  return {d[static_cast<std::size_t>(p)][0] / d[static_cast<std::size_t>(p)][2],
          d[static_cast<std::size_t>(p)][1] / d[static_cast<std::size_t>(p)][2]};
}

/// Clamped uniform knot vector for n + 1 points of degree p.
std::vector<double> clamped(int n, int p) {
  std::vector<double> U(static_cast<std::size_t>(p + 1), 0.0);
  const int inner = n - p;
  for (int i = 1; i <= inner; ++i) U.push_back(static_cast<double>(i) / (inner + 1));
  U.insert(U.end(), static_cast<std::size_t>(p + 1), 1.0);
  return U;
}

NurbsCurve polyline(const std::vector<Point2>& pts) {
  return {1, clamped(static_cast<int>(pts.size()) - 1, 1), pts, {}};
}

/// Cubic bump over one period dx = 0.5: up to 0.3, down to 0, back to the start height 0.05.
NurbsCurve bump() {
  return {3, {0, 0, 0, 0, 0.5, 1, 1, 1, 1}, {{0.0, 0.05}, {0.1, 0.25}, {0.25, 0.35}, {0.4, -0.05}, {0.5, 0.05}}, {}};
}

/// Mushroom (stem x in [0.2, 0.3] up to 0.2, cap x in [0.1, 0.4] between 0.2 and 0.3) on a base at y = 0.05.
const std::vector<Point2> kMushroom{{0.0, 0.05}, {0.2, 0.05}, {0.2, 0.2}, {0.1, 0.2}, {0.1, 0.3},
                                    {0.4, 0.3},  {0.4, 0.2},  {0.3, 0.2}, {0.3, 0.05}, {0.5, 0.05}};
/// The same geometry, starting at (0.4, 0.25) on the right side of the cap.
const std::vector<Point2> kMushroomShifted{{0.4, 0.25}, {0.4, 0.2},  {0.3, 0.2}, {0.3, 0.05}, {0.5, 0.05}, {0.7, 0.05},
                                           {0.7, 0.2},  {0.6, 0.2},  {0.6, 0.3}, {0.9, 0.3},  {0.9, 0.25}};

/// The equivalent closed body: the curve, then down to -0.1 at the end, back to the start below, up.
ProfileBody closedBody(const NurbsProfile& p, const std::string& mat) {
  Contour2D c;
  c.segments = p.pieces();
  const Point2 a = p.pieces().front().P.front(), b = p.pieces().back().P.back();
  for (const auto& seg : Contour2D::polygon({b, {b[0], -0.1}, {a[0], -0.1}, a}).segments) c.segments.push_back(seg);
  c.segments.pop_back();  // polygon closes back to b; replace the closing edge a -> b by nothing
  return {{c}, mat, mat};
}

double projDist(const Vec3D& a, const Vec3D& b) { return (a * a.transpose() - b * b.transpose()).cwiseAbs().maxCoeff(); }

void expectSameSection(const LineSection& a, const LineSection& b, double y) {
  const auto ia = geometry::materialBoundaries(a), ib = geometry::materialBoundaries(b);
  ASSERT_EQ(ia.size(), ib.size()) << "y = " << y;
  for (std::size_t k = 0; k < ia.size(); ++k) {
    EXPECT_NEAR(a.breaks[ia[k]], b.breaks[ib[k]], 1e-13) << "y = " << y;
    EXPECT_EQ(a.materials[ia[k]], b.materials[ib[k]]) << "y = " << y;
    EXPECT_LT(projDist(a.normals[ia[k]], b.normals[ib[k]]), 1e-13) << "y = " << y;
  }
  if (ia.empty()) EXPECT_EQ(a.materials.front(), b.materials.front()) << "y = " << y;
}

template <class F>
void expectCode(F f, const std::string& code) {
  try {
    f();
    ADD_FAILURE() << "no error, expected " << code;
  } catch (const Error& e) {
    EXPECT_NE(std::string(e.what()).find("[" + code + "]"), std::string::npos) << e.what();
  }
}

}  // namespace

TEST(NurbsProfile, PiecesReproduceTheCurve) {
  // rational cubic with a double interior knot (C1) and weights
  const NurbsCurve c{3, {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
                     {{0.0, 0.1}, {0.1, 0.3}, {0.2, 0.25}, {0.25, 0.0}, {0.35, 0.05}, {0.45, 0.2}, {0.5, 0.1}},
                     {1.0, 0.7, 1.3, 1.0, 0.8, 1.1, 1.0}};
  const NurbsProfile p(c, 0.5, "si", "air");
  ASSERT_EQ(p.pieces().size(), 3u);  // spans [0, 0.3], [0.3, 0.7], [0.7, 1]
  const double spans[4] = {0.0, 0.3, 0.7, 1.0};
  double err = 0;
  for (std::size_t s = 0; s < 3; ++s)
    for (int i = 0; i <= 20; ++i) {
      const double t = i / 20.0, u = spans[s] + t * (spans[s + 1] - spans[s]);
      const Point2 a = p.pieces()[s].point(t), b = deBoor(c, u);
      err = std::max({err, std::abs(a[0] - b[0]), std::abs(a[1] - b[1])});
    }
  EXPECT_LT(err, 1e-14);
}

TEST(NurbsProfile, YExtremaOfPolynomialAndRationalPieces) {
  // cubic: dy/dt = 3 sum dy_j B_j^2 is a quadratic, roots in closed form
  const BezierSegment s{{{0.0, 0.0}, {0.1, 0.9}, {0.2, -0.6}, {0.3, 0.4}}, {}};
  const double d0 = 0.9, d1 = -1.5, d2 = 1.0;  // differences of y
  const double A = d0 - 2 * d1 + d2, B = 2 * (d1 - d0), C = d0;  // quadratic in t
  const double disc = std::sqrt(B * B - 4 * A * C);
  std::vector<double> ref{(-B - disc) / (2 * A), (-B + disc) / (2 * A)};
  std::sort(ref.begin(), ref.end());
  const std::vector<double> t = geometry::yExtremaParameters(s);
  ASSERT_EQ(t.size(), 2u);
  for (std::size_t i = 0; i < 2; ++i) EXPECT_NEAR(t[i], ref[i], 1e-14);
  // circular arcs (rational quadratics, w1 = cos of half the opening): the top y = r at angle pi / 2
  for (const auto& [a0, a1] : {std::pair{30.0, 150.0}, std::pair{20.0, 130.0}, std::pair{95.0, 60.0}}) {
    const double t0 = a0 * pi / 180, t1 = a1 * pi / 180, h = 0.5 * (t1 - t0), r = 0.7;
    const Point2 P0{r * std::cos(t0), r * std::sin(t0)}, P2{r * std::cos(t1), r * std::sin(t1)};
    const double m = 0.5 * (t0 + t1), l = r / std::cos(h);
    const BezierSegment arc{{P0, {l * std::cos(m), l * std::sin(m)}, P2}, {1.0, std::cos(h), 1.0}};
    const std::vector<double> e = geometry::yExtremaParameters(arc);
    ASSERT_EQ(e.size(), 1u) << a0;
    EXPECT_NEAR(arc.point(e[0])[1], r, 1e-14) << a0;
    EXPECT_NEAR(arc.point(e[0])[0], 0.0, 1e-12) << a0;
  }
  // monotonic and constant pieces have none
  EXPECT_TRUE(geometry::yExtremaParameters(BezierSegment{{{0, 0}, {0.1, 0.2}, {0.3, 0.3}}, {}}).empty());
  EXPECT_TRUE(geometry::yExtremaParameters(BezierSegment{{{0, 0.2}, {0.1, 0.2}, {0.3, 0.2}}, {}}).empty());
}

TEST(NurbsProfile, CriticalHeightsAndRange) {
  const NurbsProfile m(polyline(kMushroom), 0.5, "si", "air");
  EXPECT_EQ(m.criticalHeights(1e-12), (std::vector<double>{0.05, 0.2, 0.3}));
  EXPECT_DOUBLE_EQ(m.yMin(), 0.05);
  EXPECT_DOUBLE_EQ(m.yMax(), 0.3);
  // the bump: interior extrema of the two cubic pieces, exact against a fine scan
  const NurbsProfile b(bump(), 0.5, "si", "air");
  double lo = 1e300, hi = -1e300;
  for (const auto& s : b.pieces())
    for (int i = 0; i <= 200000; ++i) {
      const double y = s.point(i / 200000.0)[1];
      lo = std::min(lo, y);
      hi = std::max(hi, y);
    }
  EXPECT_NEAR(b.yMin(), lo, 1e-11);
  EXPECT_NEAR(b.yMax(), hi, 1e-11);
  const std::vector<double> h = b.criticalHeights(1e-12);
  EXPECT_NE(std::find(h.begin(), h.end(), b.yMin()), h.end());
  EXPECT_NE(std::find(h.begin(), h.end(), b.yMax()), h.end());
  EXPECT_TRUE(std::is_sorted(h.begin(), h.end()));
  // provider: interior band breaks only
  const NurbsProfileLineProvider pr(std::make_shared<const NurbsProfile>(polyline(kMushroom), 0.5, "si", "air"), 0.4);
  EXPECT_EQ(pr.bandBreaks(), (std::vector<double>{0.05, 0.2, 0.3}));
  const NurbsProfileLineProvider tight(std::make_shared<const NurbsProfile>(polyline(kMushroom), 0.5, "si", "air"), 0.3);
  EXPECT_EQ(tight.bandBreaks(), (std::vector<double>{0.05, 0.2}));
}

TEST(NurbsProfile, MaterialRuleOnTheMushroom) {
  const NurbsProfile m(polyline(kMushroom), 0.5, "si", "air");
  auto check = [&](double y, const std::vector<double>& br, const std::vector<std::string>& mat) {
    const LineSection s = m.section(y);
    ASSERT_EQ(s.breaks.size(), br.size()) << y;
    for (std::size_t k = 0; k < br.size(); ++k) EXPECT_NEAR(s.breaks[k], br[k], 1e-15) << y;
    EXPECT_EQ(s.materials, mat) << y;
  };
  check(0.01, {0.0, 0.5}, {"si"});                              // under the curve
  check(0.1, {0.0, 0.2, 0.3, 0.5}, {"air", "si", "air"});        // stem
  check(0.25, {0.0, 0.1, 0.4, 0.5}, {"air", "si", "air"});       // cap (overhang)
  check(0.35, {0.0, 0.5}, {"air"});                             // above
  // normals: vertical walls have n = (+-1, 0)
  const LineSection s = m.section(0.25);
  EXPECT_NEAR(std::abs(s.normals[1](0)), 1.0, 1e-15);
}

TEST(NurbsProfile, EqualsTheClosedBody) {
  geometry::LineSectionOptions opt;
  // (the closed body needs a start without overhang across the vertical closing edges: not the shifted mushroom,
  // whose cap reaches left of its start x = 0.4; that one is compared with the unshifted description below)
  for (const NurbsCurve& c : {bump(), polyline(kMushroom)}) {
    const NurbsProfile p(c, 0.5, "si", "air");
    const std::vector<ProfileBody> body{closedBody(p, "si")};
    const std::vector<double> crit = p.criticalHeights(1e-12);
    for (int i = 1; i < 60; ++i) {
      const double y = p.yMin() - 0.02 + (p.yMax() - p.yMin() + 0.04) * i / 60.0;
      if (y <= 0 || std::any_of(crit.begin(), crit.end(), [&](double h) { return std::abs(h - y) < 1e-9; })) continue;
      expectSameSection(p.section(y), geometry::lineSection(body, y, 0.5, "air", opt), y);
    }
  }
  // the shifted description is the same geometry
  const NurbsProfile a(polyline(kMushroom), 0.5, "si", "air"), b(polyline(kMushroomShifted), 0.5, "si", "air");
  for (double y : {0.07, 0.1, 0.21, 0.24, 0.26, 0.29}) expectSameSection(a.section(y), b.section(y), y);
}

TEST(NurbsProfile, CrossingAtTheSeamCountsOnce) {
  // the bump starts at y0 = 0.05 going up; the cut y = y0 crosses at x = 0 (seam) and once going down
  const NurbsProfile b(bump(), 0.5, "si", "air");
  const LineSection s = b.section(0.05);
  const auto bd = geometry::materialBoundaries(s);
  ASSERT_EQ(bd.size(), 2u);
  EXPECT_EQ(bd[0], 0u);  // the boundary at 0 = dx
  EXPECT_EQ(s.materials.front(), "si");
  EXPECT_EQ(s.materials.back(), "air");
  // the same with a kink at the seam: a triangle wave starting at its valley
  const NurbsProfile t(polyline({{0.0, 0.0}, {0.25, 0.2}, {0.5, 0.0}}), 0.5, "si", "air");
  EXPECT_EQ(geometry::materialBoundaries(t.section(0.1)).size(), 2u);
}

TEST(NurbsProfile, Checks) {
  auto make = [](const NurbsCurve& c, double dx = 0.5) { return NurbsProfile(c, dx, "si", "air"); };
  NurbsCurve c = bump();
  c.degree = 0;
  expectCode([&] { make(c); }, "NURBS_INPUT");
  c = bump();
  c.knots.pop_back();
  expectCode([&] { make(c); }, "NURBS_INPUT");
  c = bump();
  c.knots[1] = 0.1;  // not clamped
  expectCode([&] { make(c); }, "NURBS_INPUT");
  c = bump();
  c.weights = {1, 1, 0, 1, 1};
  expectCode([&] { make(c); }, "NURBS_INPUT");
  c = bump();
  c.points.back()[1] = 0.06;
  expectCode([&] { make(c); }, "NURBS_NOT_PERIODIC");
  expectCode([&] { make(bump(), 0.4); }, "NURBS_PERIOD");
  // reversed orientation: from (0.5, 0.05) to (0, 0.05)
  c = bump();
  std::reverse(c.points.begin(), c.points.end());
  expectCode([&] { make(c); }, "NURBS_PERIOD");
  // a loop crossing itself
  expectCode([&] { make(polyline({{0.0, 0.1}, {0.3, 0.1}, {0.3, 0.3}, {0.2, 0.3}, {0.2, 0.0}, {0.5, 0.1}})); },
             "NURBS_SELF_INTERSECTION");
  // a cap reaching into the periodic image
  expectCode([&] {
    make(polyline({{0.0, 0.1}, {0.2, 0.1}, {0.2, 0.3}, {0.8, 0.3}, {0.8, 0.2}, {0.5, 0.2}, {0.5, 0.1}}));
  }, "NURBS_SELF_INTERSECTION");
  // touching the image only at the seam is fine (the mushroom touches nothing)
  EXPECT_NO_THROW(make(polyline(kMushroom)));
}
