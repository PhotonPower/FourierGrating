// Line sections of 1D profiles (Plan section "M23: Linienschnitt", "Normalenfeld im Linienschnitt"):
// crossings of rational Bezier segments, parity rule, periodic folding and seams, overlaps,
// the normal field eq:cad1d_N against its closed form.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "fg/core/error.hpp"
#include "fg/geometry/line_section.hpp"
#include "fg/geometry/normal_field.hpp"

using namespace fg;
using geometry::Contour2D;
using geometry::LineSectionOptions;
using geometry::ProfileBody;

namespace {

constexpr double pi = std::numbers::pi;

ProfileBody box(double x0, double x1, double y0, double y1, const std::string& mat, bool cw = false) {
  std::vector<geometry::Point2> v{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  if (cw) std::reverse(v.begin(), v.end());
  return {{Contour2D::polygon(v)}, mat, mat};
}

double thetaOf(const Vec3D& n) {
  double t = std::fmod(std::atan2(n(1), n(0)), pi);
  return t < 0 ? t + pi : t;
}

/// Projector distance of two in-plane normals.
double projDist(const Vec3D& a, const Vec3D& b) {
  return (a * a.transpose() - b * b.transpose()).cwiseAbs().maxCoeff();
}

}  // namespace

TEST(LineSection, CircleCrossingsAreExact) {
  const double r = 0.3;
  const Contour2D c = Contour2D::circle({0.1, 0.2}, r);
  for (double y : {0.2, 0.0, 0.45, 0.2 + r - 1e-9, 0.2 - r + 1e-7}) {
    std::vector<geometry::Crossing> cr;
    for (const auto& s : c.segments) geometry::lineCrossings(s, y, cr);
    ASSERT_EQ(cr.size(), 2u) << y;
    std::sort(cr.begin(), cr.end(), [](auto& a, auto& b) { return a.x < b.x; });
    const double h = std::sqrt(r * r - (y - 0.2) * (y - 0.2));
    // rounding of y (~1e-16) divided by the slope dy/dx = h / r of the circle at the crossing
    const double tol = 2e-15 + 4e-16 * r / h;
    EXPECT_NEAR(cr[0].x, 0.1 - h, tol) << y;
    EXPECT_NEAR(cr[1].x, 0.1 + h, tol) << y;
    for (const auto& k : cr) {  // radial normal (Plan eq:cad1d_n), sign arbitrary
      const Vec3D rad((k.x - 0.1) / r, (y - 0.2) / r, 0.0);
      EXPECT_LT(projDist(k.n, rad), 1e-12) << y;
    }
  }
}

TEST(LineSection, ParityRuleWithHoleAndOrientation) {
  // a frame: outer box and an inner box (any orientation): two intervals at mid height
  ProfileBody b = box(0.1, 0.7, 0.0, 1.0, "si");
  b.loops.push_back(box(0.3, 0.5, 0.2, 0.8, "", true).loops.front());
  b.loops.front() = b.loops.front().reversed();
  const LineSection s = geometry::lineSection({b}, 0.5, 1.0, "air", {});
  const std::vector<double> br{0.0, 0.1, 0.3, 0.5, 0.7, 1.0};
  ASSERT_EQ(s.breaks.size(), br.size());
  for (std::size_t k = 0; k < br.size(); ++k) EXPECT_NEAR(s.breaks[k], br[k], 1e-15);
  EXPECT_EQ(s.materials, (std::vector<std::string>{"air", "si", "air", "si", "air"}));
  // outside the hole's height: one interval
  const LineSection t = geometry::lineSection({b}, 0.1, 1.0, "air", {});
  EXPECT_EQ(t.materials, (std::vector<std::string>{"air", "si", "air"}));
}

TEST(LineSection, FoldingAndSeams) {
  const double d = 0.5;
  // a line over the cell boundary: [0.4, 0.6] -> [0, 0.1] u [0.4, 0.5]; break 0 is no boundary
  const LineSection s = geometry::lineSection({box(0.4, 0.6, 0, 1, "si")}, 0.5, d, "air", {});
  ASSERT_EQ(s.breaks.size(), 4u);
  EXPECT_NEAR(s.breaks[1], 0.1, 1e-15);
  EXPECT_NEAR(s.breaks[2], 0.4, 1e-15);
  EXPECT_EQ(s.materials, (std::vector<std::string>{"si", "air", "si"}));
  EXPECT_EQ(geometry::materialBoundaries(s), (std::vector<std::size_t>{1, 2}));
  // shifting by a period changes nothing
  const LineSection s2 = geometry::lineSection({box(-0.1, 0.1, 0, 1, "si")}, 0.5, d, "air", {});
  for (std::size_t k = 0; k < s.breaks.size(); ++k) EXPECT_NEAR(s.breaks[k], s2.breaks[k], 1e-15);
  // exactly one period wide (residual layer): homogeneous, no boundary
  const LineSection h = geometry::lineSection({box(0.2, 0.7, 0, 1, "si")}, 0.5, d, "air", {});
  EXPECT_EQ(h.materials, (std::vector<std::string>{"si"}));
  EXPECT_TRUE(geometry::materialBoundaries(h).empty());
  // two images touching (two boxes of one body, each half a period)
  ProfileBody two = box(0.0, 0.25, 0, 1, "si");
  two.loops.push_back(box(0.25, 0.5, 0, 1, "").loops.front());
  EXPECT_EQ(geometry::lineSection({two}, 0.5, d, "air", {}).materials, (std::vector<std::string>{"si"}));
}

TEST(LineSection, Overlaps) {
  const double d = 0.5;
  // wider than the period: overlaps its own image
  const ProfileBody wide = box(0.0, 0.7, 0, 1, "si");
  EXPECT_THROW((void)geometry::lineSection({wide}, 0.5, d, "air", {}), Error);
  LineSectionOptions u;
  u.periodicUnion = true;
  EXPECT_EQ(geometry::lineSection({wide}, 0.5, d, "air", u).materials, (std::vector<std::string>{"si"}));
  geometry::LineSectionFindings f;
  (void)geometry::lineSection({wide}, 0.5, d, "air", {}, &f);
  EXPECT_EQ(f.selfOverlap, (std::vector<int>{0}));
  // two bodies
  const std::vector<ProfileBody> two{box(0.1, 0.3, 0, 1, "si"), box(0.2, 0.4, 0, 1, "sio2")};
  EXPECT_THROW((void)geometry::lineSection(two, 0.5, d, "air", {}), Error);
  LineSectionOptions p;
  p.priority = true;
  const LineSection s = geometry::lineSection(two, 0.5, d, "air", p);
  EXPECT_EQ(s.materials, (std::vector<std::string>{"air", "si", "sio2", "air"}));
  EXPECT_NEAR(s.breaks[2], 0.3, 1e-15);
  geometry::LineSectionFindings g;
  (void)geometry::lineSection(two, 0.5, d, "air", {}, &g);
  ASSERT_EQ(g.bodyOverlap.size(), 1u);
  EXPECT_EQ(g.bodyOverlap[0], (std::array<int, 2>{1, 0}));
  // touching bodies keep their interface
  const LineSection t = geometry::lineSection({box(0.1, 0.2, 0, 1, "si"), box(0.2, 0.3, 0, 1, "sio2")}, 0.5, d,
                                              "air", {});
  EXPECT_EQ(t.materials, (std::vector<std::string>{"air", "si", "sio2", "air"}));
}

TEST(LineSection, NormalAngleInterpolates) {
  // trapezoid with 15 degree walls: boundary normals exact, angle linear in between (Plan eq:cad1d_N)
  const double d = 0.3, w = std::tan(15 * pi / 180);
  ProfileBody t{{Contour2D::polygon({{0.1, 0.0}, {0.2, 0.0}, {0.2 - 0.1 * w, 0.1}, {0.1 + 0.1 * w, 0.1}})}, "si",
                "T"};
  const LineSection s = geometry::lineSection({t}, 0.05, d, "air", {});
  const auto b = geometry::materialBoundaries(s);
  ASSERT_EQ(b.size(), 2u);
  for (std::size_t k : b) {
    const Vec3D N = geometry::lineNormal(s, d, s.breaks[k]);
    EXPECT_LT(projDist(N, s.normals[k]), 1e-15);
  }
  const double tl = thetaOf(s.normals[b[0]]), tr = thetaOf(s.normals[b[1]]);
  EXPECT_NEAR(std::abs(tl - tr), pi - 2 * 15 * pi / 180, 1e-14);  // 165 and 15 degrees
  // inside the line the angle passes through 0 / pi (N = x), the short way round
  const double xm = 0.5 * (s.breaks[b[0]] + s.breaks[b[1]]);
  EXPECT_LT(projDist(geometry::lineNormal(s, d, xm), Vec3D(1, 0, 0)), 1e-14);
}

TEST(LineSection, NormalFieldCoefficientsClosedForm) {
  // boundaries at a, b with angles ta, tb; theta linear on [a, b] and on [b, a + d]
  const double d = 0.4, a = 0.07, b = 0.23, ta = 0.3, tb = 1.2;
  LineSection s;
  s.breaks = {0.0, a, b, d};
  s.materials = {"air", "si", "air"};
  s.normals = {Vec3D(1, 0, 0), Vec3D(std::cos(ta), std::sin(ta), 0), Vec3D(std::cos(tb), std::sin(tb), 0)};
  const IndexMap map(6, 0);
  const NormalFieldCoefficients nf = normalField(s, d, map);
  EXPECT_FALSE(nf.hasSigned);
  // closed form: (1/d) int_x0^x1 e^{i(s 2 theta(x) - K x)} dx with theta = t0 + k (x - x0)
  auto piece = [&](double x0, double x1, double t0, double t1, int sgn2, double K) {
    const double k = (t1 - t0) / (x1 - x0), al = sgn2 * 2 * k - K;
    const cplx ph = std::exp(cplx(0, sgn2 * 2 * t0 - K * x0));
    const double L = x1 - x0;
    const cplx I = std::abs(al) < 1e-300 ? cplx(L, 0) : (std::exp(cplx(0, al * L)) - 1.0) / cplx(0, al);
    return ph * I / d;
  };
  auto interval = [&](int sgn2, double K) {
    double dl = ta - tb;  // from b to a + d: reduced difference
    while (dl > pi / 2) dl -= pi;
    while (dl <= -pi / 2) dl += pi;
    double du = tb - ta;
    while (du > pi / 2) du -= pi;
    while (du <= -pi / 2) du += pi;
    return piece(a, b, ta, ta + du, sgn2, K) + piece(b, a + d, tb, tb + dl, sgn2, K);
  };
  for (int p = -12; p <= 12; ++p) {
    const double K = 2 * pi * p / d;
    const cplx c0 = (p == 0) ? cplx(1, 0) : cplx(0, 0);
    const cplx ep = interval(1, K), em = interval(-1, K);
    const cplx xx = 0.5 * c0 + 0.25 * (ep + em), yy = 0.5 * c0 - 0.25 * (ep + em), xy = (ep - em) / cplx(0, 4);
    EXPECT_LT(std::abs(nf.NN[0][0](p, 0) - xx), 1e-14) << p;
    EXPECT_LT(std::abs(nf.NN[1][1](p, 0) - yy), 1e-14) << p;
    EXPECT_LT(std::abs(nf.NN[0][1](p, 0) - xy), 1e-14) << p;
    EXPECT_LT(std::abs(nf.NN[1][0](p, 0) - xy), 1e-14) << p;
  }
}

TEST(LineSection, VerticalWallsGiveConstantNormal) {
  const LineSection s = geometry::lineSection({box(0.1, 0.2, 0, 1, "si")}, 0.5, 0.4, "air", {});
  const NormalFieldCoefficients nf = normalField(s, 0.4, IndexMap(5, 0));
  for (int p = -10; p <= 10; ++p) {
    EXPECT_LT(std::abs(nf.NN[0][0](p, 0) - (p == 0 ? 1.0 : 0.0)), 1e-15) << p;
    EXPECT_LT(std::abs(nf.NN[1][1](p, 0)), 1e-15) << p;
    EXPECT_LT(std::abs(nf.NN[0][1](p, 0)), 1e-15) << p;
  }
}

TEST(LineSection, ProviderStaysOffTheBandEnds) {
  auto bodies = std::make_shared<const std::vector<ProfileBody>>(std::vector<ProfileBody>{box(0.1, 0.2, 0.3, 0.5, "si")});
  ProfileLineProvider pr(bodies, 0.3, 0.2, 0.4, "air", {});
  EXPECT_DOUBLE_EQ(pr.thickness_um(), 0.2);
  for (double y : {0.0, 0.1, 0.2}) EXPECT_EQ(pr.section(y).materials, (std::vector<std::string>{"air", "si", "air"}));
}
