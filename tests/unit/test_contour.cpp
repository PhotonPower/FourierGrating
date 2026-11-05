// Contours and contour-integral coefficients (Plan section "M13", eq:contour, eq:circle).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/core/error.hpp"
#include "fg/fourier/contour.hpp"
#include "fg/io/config.hpp"

using fg::cplx;
using fg::geometry::Contour2D;
using fg::geometry::Point2;
using fg::geometry::Region;

namespace {

Region region(Contour2D c, std::vector<Contour2D> holes = {}) {
  Region r;
  r.outer = std::move(c);
  r.holes = std::move(holes);
  r.normalize(1e-9);
  return r;
}

/// Plan eq:circle for an ellipse x = c + A u, |u| <= 1: f = det(A) 2 pi J1(|A^T k|) / (|A^T k| dx dz) e^{-i k c}
cplx ellipseExact(Point2 c, double a, double b, double angle, double dx, double dz, int p, int q) {
  const double kx = 2 * fg::PI * p / dx, kz = 2 * fg::PI * q / dz;
  const double ca = std::cos(angle), sa = std::sin(angle);
  if (p == 0 && q == 0) return fg::PI * a * b / (dx * dz);
  const double u = a * (ca * kx + sa * kz), v = b * (-sa * kx + ca * kz), kap = std::hypot(u, v);
  return 2 * fg::PI * a * b * std::cyl_bessel_j(1.0, kap) / (kap * dx * dz) * std::exp(-fg::I * (kx * c[0] + kz * c[1]));
}

double maxDiff(const fg::fourier::CoefficientGrid& a, const fg::fourier::CoefficientGrid& b) {
  double d = 0;
  for (int p = -a.Px(); p <= a.Px(); ++p)
    for (int q = -a.Pz(); q <= a.Pz(); ++q) d = std::max(d, std::abs(a(p, q) - b(p, q)));
  return d;
}

}  // namespace

TEST(Contour, DegreeOneEqualsPolygonFormula) {
  const std::vector<std::array<double, 2>> v{{0.1, 0.05}, {0.62, 0.12}, {0.55, 0.4}, {0.33, 0.27}, {0.12, 0.58}};
  const Region r = region(Contour2D::polygon(v));
  const double dx = 0.7, dz = 0.65;
  const auto g = fg::fourier::regionCoefficients(r, dx, dz, 14, 12);
  double d = 0;
  for (int p = -14; p <= 14; ++p)
    for (int q = -12; q <= 12; ++q) d = std::max(d, std::abs(g(p, q) - fg::fourier::polygonCoefficient(v, dx, dz, p, q)));
  EXPECT_LT(d, 1e-14);
}

TEST(Contour, CircleEqualsBessel) {
  const Point2 c{0.4, 0.45};
  const Region r = region(Contour2D::circle(c, 0.3));
  const auto g = fg::fourier::regionCoefficients(r, 1.0, 1.0, 60, 60);
  double d = 0;
  for (int p = -60; p <= 60; ++p)
    for (int q = -60; q <= 60; ++q) d = std::max(d, std::abs(g(p, q) - ellipseExact(c, 0.3, 0.3, 0, 1.0, 1.0, p, q)));
  EXPECT_LT(d, 1e-13);
  EXPECT_NEAR(r.area(), fg::PI * 0.09, 1e-15);
}

TEST(Contour, EllipseIsStretchedCircle) {
  const Point2 c{0.3, 0.2};
  const double dx = 0.8, dz = 0.6, a = 0.2, b = 0.11, ang = 0.5236;
  const auto g = fg::fourier::regionCoefficients(region(Contour2D::ellipse(c, a, b, ang)), dx, dz, 20, 20);
  double d = 0;
  for (int p = -20; p <= 20; ++p)
    for (int q = -20; q <= 20; ++q) d = std::max(d, std::abs(g(p, q) - ellipseExact(c, a, b, ang, dx, dz, p, q)));
  EXPECT_LT(d, 1e-13);
}

TEST(Contour, PeriodShiftAndHole) {
  const double dx = 0.5, dz = 0.4;
  const auto g0 = fg::fourier::regionCoefficients(region(Contour2D::circle({0.2, 0.15}, 0.12)), dx, dz, 16, 16);
  const auto g1 = fg::fourier::regionCoefficients(region(Contour2D::circle({0.2 + dx, 0.15 - dz}, 0.12)), dx, dz, 16, 16);
  EXPECT_LT(maxDiff(g0, g1), 1e-13);
  // annulus = outer disc - inner disc (hole = negative area)
  const auto ring = fg::fourier::regionCoefficients(
      region(Contour2D::circle({0.2, 0.15}, 0.12), {Contour2D::circle({0.22, 0.16}, 0.05)}), dx, dz, 16, 16);
  const auto inner = fg::fourier::regionCoefficients(region(Contour2D::circle({0.22, 0.16}, 0.05)), dx, dz, 16, 16);
  double d = 0;
  for (int p = -16; p <= 16; ++p)
    for (int q = -16; q <= 16; ++q) d = std::max(d, std::abs(ring(p, q) - (g0(p, q) - inner(p, q))));
  EXPECT_LT(d, 1e-14);
}

TEST(Contour, AdaptivityIndependentOfN0) {
  const Region r = region(Contour2D::ellipse({0.31, 0.27}, 0.2, 0.13, 0.3));
  const auto a = fg::fourier::regionCoefficients(r, 0.7, 0.6, 24, 24, 16);
  const auto b = fg::fourier::regionCoefficients(r, 0.7, 0.6, 24, 24, 80);
  const auto c = fg::fourier::regionCoefficients(r, 0.7, 0.6, 24, 24);
  EXPECT_LT(maxDiff(a, b), 1e-14);
  EXPECT_LT(maxDiff(a, c), 1e-14);
}

TEST(Contour, NurbsCircleEqualsExactCircle) {
  // standard 9-point rational quadratic NURBS circle (Piegl-Tiller, Ex. 7.2), radius 0.2 at (0.3, 0.3)
  const double h = std::sqrt(0.5), R = 0.2, cx = 0.3, cz = 0.3;
  const std::vector<Point2> P{{cx + R, cz},     {cx + R, cz + R}, {cx, cz + R}, {cx - R, cz + R}, {cx - R, cz},
                              {cx - R, cz - R}, {cx, cz - R},     {cx + R, cz - R}, {cx + R, cz}};
  const std::vector<double> w{1, h, 1, h, 1, h, 1, h, 1};
  const std::vector<double> U{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
  const Contour2D n = Contour2D::nurbs(2, U, P, w);
  ASSERT_EQ(n.segments.size(), 4u);
  const auto a = fg::fourier::regionCoefficients(region(n), 1.0, 1.0, 20, 20);
  const auto b = fg::fourier::regionCoefficients(region(Contour2D::circle({cx, cz}, R)), 1.0, 1.0, 20, 20);
  EXPECT_LT(maxDiff(a, b), 1e-14);
  // cubic B-spline (non-rational, interior knots of multiplicity 1): pieces join, area consistent
  const std::vector<Point2> Q{{0.1, 0.1}, {0.5, 0.05}, {0.6, 0.4}, {0.3, 0.6}, {0.05, 0.4}, {0.1, 0.1}};
  const Contour2D cub = Contour2D::nurbs(3, {0, 0, 0, 0, 0.4, 0.7, 1, 1, 1, 1}, Q, {});
  EXPECT_EQ(cub.segments.size(), 3u);
  EXPECT_NO_THROW(cub.validate());
}

TEST(Contour, SplitAndContains) {
  const Contour2D c = Contour2D::ellipse({0.5, 0.5}, 0.3, 0.1, 0.2);
  fg::geometry::BezierSegment l, r;
  c.segments[1].split(0.3, l, r);
  for (double t : {0.0, 0.4, 1.0}) {
    const Point2 a = l.point(t), b = c.segments[1].point(0.3 * t);
    EXPECT_NEAR(a[0], b[0], 1e-15);
    EXPECT_NEAR(a[1], b[1], 1e-15);
  }
  const Region reg = region(c.reversed(), {Contour2D::circle({0.5, 0.5}, 0.05)});  // clockwise input
  EXPECT_GT(reg.outer.signedArea(), 0);
  EXPECT_TRUE(reg.contains(0.7, 0.55));
  EXPECT_FALSE(reg.contains(0.5, 0.5));  // in the hole
  EXPECT_FALSE(reg.contains(0.5, 0.65));
  EXPECT_NEAR(reg.area(), fg::PI * (0.03 - 0.0025), 1e-15);
  Contour2D broken = c;
  broken.segments[2].P.back()[0] += 1e-6;
  EXPECT_THROW(broken.validate(), fg::Error);
}

namespace {
fg::Result solveRegions(const std::string& regions) {
  const auto c = fg::io::parseConfig(R"({"dimension": "2D", "truncation": {"Nx": 2, "Nz": 2},
    "period_um": {"x": 0.5, "z": 0.5}, "incidence": {"lambda_um": 0.6, "phi_deg": 0, "psi_deg": 0, "delta_deg": 0},
    "materials": {"air": {"n": 1}, "si": {"n": 3.5}}, "superstrate": "air", "substrate": "air",
    "layers": [{"type": "contour2d", "thickness_um": 0.1, "background": "air", "regions": [)" + regions + "]}]}");
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
}
}  // namespace

TEST(Contour, Contour2DDisjointRegions) {
  const auto ok = solveRegions(R"({"material": "si", "circle": {"center_um": [0.1, 0.1], "radius_um": 0.08}},
                                  {"material": "si", "circle": {"center_um": [0.3, 0.3], "radius_um": 0.08}})");
  EXPECT_NEAR(ok.diag.energyBalance, 1.0, 1e-12);
  // crossing boundaries
  EXPECT_THROW(solveRegions(R"({"material": "si", "circle": {"center_um": [0.2, 0.2], "radius_um": 0.1}},
                               {"material": "si", "circle": {"center_um": [0.3, 0.2], "radius_um": 0.1}})"),
               fg::Error);
  // one region inside another
  EXPECT_THROW(solveRegions(R"({"material": "si", "circle": {"center_um": [0.2, 0.2], "radius_um": 0.15}},
                               {"material": "si", "circle": {"center_um": [0.2, 0.2], "radius_um": 0.05}})"),
               fg::Error);
  // overlap with its own periodic image (diameter > period)
  EXPECT_THROW(solveRegions(R"({"material": "si", "circle": {"center_um": [0.25, 0.25], "radius_um": 0.3}})"), fg::Error);
  EXPECT_THROW(solveRegions(R"({"material": "si", "square": {"size_um": 0.1}})"), fg::Error);
}
