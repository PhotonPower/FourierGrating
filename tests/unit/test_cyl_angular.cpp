// Unit tests M21 (Plan section "Winkelkoeffizienten aus Querschnitten", eq:cylP_chi; "Tests und Abnahme (M21)"):
// arcs, angular coefficients chi_q(r), critical radii, rotational symmetry; Bernstein helpers.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>

#include "fg/cyl/angular.hpp"
#include "fg/cyl/coupled.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/geometry/bernstein.hpp"

using fg::cplx;
using fg::PI;
using fg::geometry::Contour2D;
using fg::geometry::Point2;
using fg::geometry::Region;
namespace cyl = fg::cyl;

namespace {

Region region(Contour2D outer, std::vector<Contour2D> holes = {}, const std::string& material = "m") {
  Region r;
  r.outer = std::move(outer);
  r.holes = std::move(holes);
  r.material = material;
  r.normalize(1e-9);
  return r;
}

/// Regular n-gon about the origin with apothem a, one vertex at angle psi0.
Contour2D ngon(int n, double a, double psi0) {
  const double R = a / std::cos(PI / n);
  std::vector<Point2> v;
  for (int k = 0; k < n; ++k) v.push_back({R * std::cos(psi0 + 2 * PI * k / n), R * std::sin(psi0 + 2 * PI * k / n)});
  return Contour2D::polygon(v);
}

/// chi_q of the regular n-gon with P = n (exact): for a < r < R the circle is inside on arcs of half width
/// theta = pi/n - acos(a/r) about the vertex directions, chi_0 = n theta / pi, chi_q = e^{-iqn psi0} sin(qn theta)/(pi q).
cplx ngonChi(int n, double a, double psi0, double r, int q) {
  const double R = a / std::cos(PI / n);
  if (r < a) return q == 0 ? 1.0 : 0.0;
  if (r > R) return 0.0;
  const double theta = PI / n - std::acos(a / r);
  if (q == 0) return n * theta / PI;
  return std::exp(cplx(0, -q * n * psi0)) * std::sin(q * n * theta) / (PI * q);
}

/// chi_q (P = 1) of the disc of radius R about (d, 0), d > 0: for |d - R| < r < d + R the arc |phi| < alpha,
/// cos alpha = (r^2 + d^2 - R^2) / (2 r d); chi_q = sin(q alpha) / (pi q).
cplx discChi(double d, double R, double r, int q) {
  if (r < R - d) return q == 0 ? 1.0 : 0.0;
  if (r <= std::abs(d - R) || r >= d + R) return 0.0;
  const double alpha = std::acos((r * r + d * d - R * R) / (2 * r * d));
  return q == 0 ? alpha / PI : std::sin(q * alpha) / (PI * q);
}

double maxError(const Region& reg, int P, int qmax, const std::vector<double>& radii,
                const std::function<cplx(double, int)>& exact) {
  double e = 0;
  for (const double r : radii) {
    const fg::VecC chi = cyl::angularCoefficients(reg, r, P, qmax);
    for (int q = 0; q <= qmax; ++q) e = std::max(e, std::abs(chi(q) - exact(r, q)));
  }
  return e;
}

std::vector<double> samples(double a, double b, int n) {
  std::vector<double> r;
  for (int i = 0; i < n; ++i) r.push_back(a + (b - a) * (i + 0.37) / n);
  return r;
}

}  // namespace

TEST(Bernstein, ProductDerivativeRoots) {
  const std::vector<double> a{1.0, -2.0, 0.5}, b{0.3, 1.1, -0.7, 2.0};
  const auto ab = fg::geometry::bernsteinProduct(a, b);
  ASSERT_EQ(ab.size(), 6u);
  for (double t : {0.0, 0.21, 0.5, 0.83, 1.0})
    EXPECT_NEAR(fg::geometry::bernsteinValue(ab, t),
                fg::geometry::bernsteinValue(a, t) * fg::geometry::bernsteinValue(b, t), 1e-14);
  const auto db = fg::geometry::bernsteinDerivative(b);
  const double h = 1e-6, t = 0.37;
  EXPECT_NEAR(fg::geometry::bernsteinValue(db, t),
              (fg::geometry::bernsteinValue(b, t + h) - fg::geometry::bernsteinValue(b, t - h)) / (2 * h), 1e-8);
  // (t - 0.2)(t - 0.7) = t^2 - 0.9 t + 0.14 in Bernstein form: c0 = 0.14, c1 = 0.14 - 0.45, c2 = 0.24
  const auto roots = fg::geometry::bernsteinRoots({0.14, 0.14 - 0.45, 0.24});
  ASSERT_EQ(roots.size(), 2u);
  EXPECT_NEAR(roots[0], 0.2, 1e-15);
  EXPECT_NEAR(roots[1], 0.7, 1e-15);
}

TEST(CylAngular, SquareAndHexagonAgainstExactArcs) {
  // square (P = 4) and regular hexagon (P = 6), also rotated by an arbitrary angle (phase of chi_q): 1e-13
  for (const auto& [n, psi0] : std::vector<std::pair<int, double>>{{4, PI / 4}, {6, 0.0}, {6, 0.3}, {4, -0.7}}) {
    const double a = 0.35, R = a / std::cos(PI / n);
    const Region reg = region(ngon(n, a, psi0));
    auto rs = samples(0.0, a, 4);
    for (double r : samples(a, R, 23)) rs.push_back(r);
    rs.push_back(1.2 * R);
    EXPECT_LT(maxError(reg, n, 6, rs, [&](double r, int q) { return ngonChi(n, a, psi0, r, q); }), 1e-13)
        << n << " " << psi0;
    // the symmetry kills every harmonic that is no multiple of n (P = 1 sees them all)
    const fg::VecC all = cyl::angularCoefficients(reg, 0.5 * (a + R), 1, 2 * n);
    for (int k = 1; k <= 2 * n; ++k)
      if (k % n) EXPECT_LT(std::abs(all(k)), 1e-15) << n << " " << k;
  }
}

TEST(CylAngular, CircleAsRationalBezierIsExact) {
  // a centred circle (four rational quadratic segments): chi_q = delta_q0 inside, 0 outside, no crossings
  const Region disc = region(Contour2D::circle({0.0, 0.0}, 0.4));
  EXPECT_LT(maxError(disc, 1, 5, {0.1, 0.39, 0.41, 0.8},
                     [](double r, int q) { return cplx(r < 0.4 && q == 0 ? 1.0 : 0.0); }),
            1e-15);
  // off centre (P = 1): crossings of the rational segments with the circle, 1e-13 against the exact arc
  const double d = 0.15, R = 0.4;
  const Region off = region(Contour2D::circle({d, 0.0}, R));
  auto rs = samples(0.0, R - d, 3);
  for (double r : samples(R - d, R + d, 31)) rs.push_back(r);
  EXPECT_LT(maxError(off, 1, 8, rs, [&](double r, int q) { return discChi(d, R, r, q); }), 1e-13);
  // the same disc as a fine polygon converges to it (consistency of the two boundary types)
  std::vector<Point2> v;
  for (int k = 0; k < 2000; ++k) v.push_back({d + R * std::cos(2 * PI * k / 2000), R * std::sin(2 * PI * k / 2000)});
  const Region poly = region(Contour2D::polygon(v));
  EXPECT_LT(maxError(poly, 1, 4, rs, [&](double r, int q) { return discChi(d, R, r, q); }), 1e-4);
}

TEST(CylAngular, RegionWithAHoleAndClockwiseInput) {
  // disc (r = 0.6) with a hexagonal hole: chi = chi_disc - chi_hexagon; orientation of the input is irrelevant
  const double a = 0.3;
  const Region holed = region(Contour2D::circle({0.0, 0.0}, 0.6), {ngon(6, a, 0.2).reversed()});
  const Region cw = region(Contour2D::circle({0.0, 0.0}, 0.6).reversed(), {ngon(6, a, 0.2)});
  auto rs = samples(0.0, 0.6, 37);
  rs.push_back(0.7);
  auto exact = [&](double r, int q) {
    return cplx(r < 0.6 && q == 0 ? 1.0 : 0.0) - ngonChi(6, a, 0.2, r, q);
  };
  EXPECT_LT(maxError(holed, 6, 5, rs, exact), 1e-13);
  EXPECT_LT(maxError(cw, 6, 5, rs, exact), 1e-13);
  // the axis inside the material: chi_q -> delta_q0 for r -> 0
  const fg::VecC chi = cyl::angularCoefficients(region(ngon(6, a, 0.2)), 1e-9, 6, 3);
  EXPECT_NEAR(chi(0).real(), 1.0, 1e-15);
  EXPECT_LT(std::abs(chi(1)), 1e-15);
}

TEST(CylAngular, CriticalRadii) {
  auto has = [](const std::vector<double>& v, double x) {
    return std::any_of(v.begin(), v.end(), [&](double y) { return std::abs(x - y) < 1e-12; });
  };
  // square: apothem (tangent point of the edges) and corner radius
  const auto sq = cyl::criticalRadii(region(ngon(4, 0.35, PI / 4)));
  EXPECT_EQ(sq.size(), 3u);
  EXPECT_EQ(sq.front(), 0.0);
  EXPECT_TRUE(has(sq, 0.35));
  EXPECT_TRUE(has(sq, 0.35 * std::sqrt(2.0)));
  EXPECT_TRUE(std::is_sorted(sq.begin(), sq.end()));
  // off-centre disc: nearest and farthest point (tangent points of the circle |x| = r) plus segment ends
  const auto disc = cyl::criticalRadii(region(Contour2D::circle({0.15, 0.0}, 0.4)));
  EXPECT_TRUE(has(disc, 0.25));
  EXPECT_TRUE(has(disc, 0.55));
  EXPECT_TRUE(std::is_sorted(disc.begin(), disc.end()));
  // a centred circle: all radii equal (|C| constant), one critical radius besides 0
  const auto c = cyl::criticalRadii(region(Contour2D::circle({0.0, 0.0}, 0.4)));
  EXPECT_EQ(c.size(), 2u);
  EXPECT_NEAR(c.back(), 0.4, 1e-12);
}

TEST(CylAngular, RotationalSymmetry) {
  const double tol = 1e-9 * 0.4;
  const std::vector<Region> hex{region(ngon(6, 0.3, 0.1))};
  EXPECT_TRUE(cyl::hasRotationSymmetry(hex, 6, tol));
  EXPECT_TRUE(cyl::hasRotationSymmetry(hex, 3, tol));
  EXPECT_TRUE(cyl::hasRotationSymmetry(hex, 2, tol));
  EXPECT_FALSE(cyl::hasRotationSymmetry(hex, 4, tol));
  EXPECT_EQ(cyl::rotationOrder(hex, 12, tol), 6);
  // one vertex moved by 1e-6 um: no symmetry left
  const double R = 0.3 / std::cos(PI / 6);
  std::vector<Point2> v;
  for (int k = 0; k < 6; ++k) v.push_back({R * std::cos(0.1 + PI * k / 3), R * std::sin(0.1 + PI * k / 3)});
  v[2][0] += 1e-6;
  const std::vector<Region> bent{region(Contour2D::polygon(v))};
  EXPECT_FALSE(cyl::hasRotationSymmetry(bent, 6, tol));
  EXPECT_EQ(cyl::rotationOrder(bent, 12, tol), 1);
  // square: 4; centred circle: any order; off-centre circle: none
  EXPECT_EQ(cyl::rotationOrder({region(ngon(4, 0.3, 0.0))}, 12, tol), 4);
  EXPECT_EQ(cyl::rotationOrder({region(Contour2D::circle({0.0, 0.0}, 0.4))}, 12, tol), 12);
  EXPECT_EQ(cyl::rotationOrder({region(Contour2D::circle({0.05, 0.0}, 0.4))}, 12, tol), 1);
  // two discs on opposite sides (P = 2), materials must match
  const Region l = region(Contour2D::circle({-0.5, 0.0}, 0.1)), r = region(Contour2D::circle({0.5, 0.0}, 0.1));
  EXPECT_EQ(cyl::rotationOrder({l, r}, 12, tol), 2);
  const Region r2 = region(Contour2D::circle({0.5, 0.0}, 0.1), {}, "other");
  EXPECT_EQ(cyl::rotationOrder({l, r2}, 12, tol), 1);
}

// ---- coupled Laurent blocks (Plan eq:cylP_Q)
namespace {

fg::cyl::SpectralGrid fiberGrid(double kmax, double Rs) {
  cyl::SpectrumOptions o;
  o.kmax = kmax;
  o.Rs = Rs;
  o.branchPoints = {1.0, 2.0};
  return cyl::spectralGrid(o);
}

double maxAbs(const fg::MatC& m) { return m.size() ? m.cwiseAbs().maxCoeff() : 0.0; }

}  // namespace

TEST(CylCoupled, CentredDiscEqualsM20PerBlock) {
  // a centred disc has chi_q = delta_q0: the blocks off the diagonal vanish, the diagonal blocks are the
  // rotationally symmetric Laurent matrices of M20 (Lommel integrals) of the harmonics m_l = m0 + l P
  const auto grid = fiberGrid(3.0, 2.0);
  const Eigen::Index nk = grid.k.size();
  cyl::AngularSection s;
  s.regions = {region(Contour2D::circle({0.0, 0.0}, 2.0))};
  s.eps = {4.0};
  s.epsBg = 1.0;
  const cyl::HarmonicMap map{1, 3, 2};  // m = -5, -2, 1, 4, 7
  const auto Q = cyl::coupledQ(grid, map, s);
  ASSERT_EQ(Q.Qyy.rows(), 5 * nk);
  double off = 0, diag = 0;
  for (int b = 0; b < 5; ++b)
    for (int bp = 0; bp < 5; ++bp) {
      if (b != bp) {
        off = std::max({off, maxAbs(Q.Qpp.block(b * nk, bp * nk, nk, nk)), maxAbs(Q.Qyy.block(b * nk, bp * nk, nk, nk))});
        continue;
      }
      const auto ref = cyl::cylQ(grid, map.m(b), {{2.0}, {4.0}, 1.0}, cyl::CylFactorization::Laurent);
      diag = std::max({diag, maxAbs(Q.Qpp.block(b * nk, b * nk, nk, nk) - ref.Qpp),
                       maxAbs(Q.Qmm.block(b * nk, b * nk, nk, nk) - ref.Qmm),
                       maxAbs(Q.Qyy.block(b * nk, b * nk, nk, nk) - ref.Qyy)});
    }
  EXPECT_LT(off, 1e-15);
  EXPECT_LT(diag, 1e-12);
  EXPECT_EQ(maxAbs(Q.Qpm), 0.0);
}

TEST(CylCoupled, RadialQuadratureConvergedForAnEccentricDisc) {
  // the square-root substitution at the tangent radii |a - x0|, a + x0: 16 nodes per panel = 32 to 1e-12
  const double a = 2.0, x0 = 0.6;
  const auto grid = fiberGrid(3.0, a + x0);
  cyl::AngularSection s;
  s.regions = {region(Contour2D::circle({x0, 0.0}, a))};
  s.eps = {4.0};
  const cyl::HarmonicMap map{1, 1, 3};
  const auto Q16 = cyl::coupledQ(grid, map, s), Q32 = cyl::coupledQ(grid, map, s, 32);
  const double scale = maxAbs(Q32.Qyy);
  std::printf("eccentric disc: |Q16 - Q32| = %.1e (Qyy), %.1e (Qpp) of %.2f\n", maxAbs(Q16.Qyy - Q32.Qyy),
              maxAbs(Q16.Qpp - Q32.Qpp), scale);
  EXPECT_LT(maxAbs(Q16.Qyy - Q32.Qyy), 1e-12 * scale);
  EXPECT_LT(maxAbs(Q16.Qpp - Q32.Qpp), 1e-12 * scale);
  EXPECT_LT(maxAbs(Q16.Qmm - Q32.Qmm), 1e-12 * scale);
  // the eccentric disc couples neighbouring harmonics
  EXPECT_GT(maxAbs(Q16.Qyy.block(0, grid.k.size(), grid.k.size(), grid.k.size())), 1e-3);
}

// ---- shape derivatives (M21 gradients, V22g): d chi_q / d theta from the crossings of the circle with the moving
// boundary, and the reverse rule of coupledQ, both against central differences
namespace {

/// Velocity (C(theta + h) - C(theta - h)) / (2 h) of the boundary of a region given at theta +- h.
cyl::BoundaryVelocity velocity(const Region& plus, const Region& minus, double h) {
  return [plus, minus, h](std::size_t c, std::size_t s, double t) {
    const auto& cp = c == 0 ? plus.outer : plus.holes[c - 1];
    const auto& cm = c == 0 ? minus.outer : minus.holes[c - 1];
    const Point2 a = cp.segments[s].point(t), b = cm.segments[s].point(t);
    return Point2{(a[0] - b[0]) / (2 * h), (a[1] - b[1]) / (2 * h)};
  };
}

/// Irregular pentagon about the origin with vertex 2 moved by (dx, dz). Irregular: in a regular one all vertex radii
/// and all tangent radii coincide, a moved vertex splits them by dx, and the discrete map (whose square-root
/// substitution then sits dx away from the singularity) changes by about dx^{3/2}: its difference quotients are off
/// by 1e-3, the analytic derivative is not.
Contour2D pentagon(double dx, double dz) {
  const double rad[5] = {1.5, 1.3, 1.45, 1.2, 1.38}, ang[5] = {0.3, 1.5, 2.6, 3.9, 5.0};
  std::vector<Point2> v;
  for (int k = 0; k < 5; ++k) v.push_back({rad[k] * std::cos(ang[k]), rad[k] * std::sin(ang[k])});
  v[2][0] += dx;
  v[2][1] += dz;
  return Contour2D::polygon(v);
}

}  // namespace

TEST(CylAngular, ShapeDerivativeAgainstDifferences) {
  const double h = 1e-6;
  double worst = 0;
  // eccentric disc: centre x0 and radius as parameters; pentagon: one vertex (with the tangent radii of its edges)
  struct Case {
    std::function<Region(double)> at;
    std::vector<double> radii;
  };
  const std::vector<Case> cases = {
      {[](double d) { return region(Contour2D::circle({0.6 + d, 0.0}, 2.0)); }, samples(1.41, 2.59, 7)},
      {[](double d) { return region(Contour2D::circle({0.6, 0.0}, 2.0 + d)); }, samples(1.41, 2.59, 7)},
      {[](double d) { return region(pentagon(d, 0.0)); }, samples(0.8, 1.6, 11)},
      {[](double d) { return region(pentagon(0.0, d)); }, samples(0.8, 1.6, 11)}};
  for (const auto& c : cases) {
    const Region r0 = c.at(0), rp = c.at(h), rm = c.at(-h);
    const auto v = velocity(rp, rm, h);
    for (const double r : c.radii) {
      const fg::VecC d = cyl::angularCoefficientsShape(r0, r, 1, 4, v);
      const fg::VecC fd = (cyl::angularCoefficients(rp, r, 1, 4) - cyl::angularCoefficients(rm, r, 1, 4)) / (2 * h);
      worst = std::max(worst, (d - fd).cwiseAbs().maxCoeff());
    }
  }
  std::printf("shape derivative of chi_q: max |analytic - fd| = %.1e\n", worst);
  EXPECT_LT(worst, 1e-8);
}

TEST(CylCoupled, ReverseRuleAgainstDifferences) {
  // dL = Re sum conj(Qbar) dQ for a fixed cotangent Qbar: permittivities (real and imaginary parts, background) and a
  // vertex of a pentagon (P = 1, harmonics -1, 0, 1) against central differences of coupledQ
  const auto grid = fiberGrid(2.0, 1.6);
  const Eigen::Index nf = 3 * grid.k.size();
  const cyl::HarmonicMap map{0, 1, 1};
  auto section = [](double dz, cplx e, cplx bg) {
    cyl::AngularSection s;
    s.regions = {region(pentagon(0.0, dz))};
    s.eps = {e};
    s.epsBg = bg;
    return s;
  };
  const cplx e0(4.0, 0.3), bg0(1.2, 0.05);
  cyl::CylQ Qb;
  auto pattern = [nf](double a) {
    fg::MatC M(nf, nf);
    for (Eigen::Index i = 0; i < nf; ++i)
      for (Eigen::Index j = 0; j < nf; ++j) M(i, j) = cplx(std::sin(a * i + 0.3 * j), std::cos(0.7 * i - a * j));
    return M;
  };
  Qb.Qpp = pattern(0.11);
  Qb.Qmm = pattern(0.23);
  Qb.Qyy = pattern(0.37);
  Qb.Qpm = Qb.Qmp = fg::MatC::Zero(nf, nf);
  auto L = [&](const cyl::AngularSection& s) {
    const auto Q = cyl::coupledQ(grid, map, s);
    return (Qb.Qpp.conjugate().cwiseProduct(Q.Qpp).sum() + Qb.Qmm.conjugate().cwiseProduct(Q.Qmm).sum() +
            Qb.Qyy.conjugate().cwiseProduct(Q.Qyy).sum()).real();
  };
  const auto bar = cyl::coupledQBar(grid, map, section(0, e0, bg0), Qb);
  const double h = 1e-6;
  auto check = [](const char* what, double analytic, double fd) {
    std::printf("coupledQ reverse %-8s analytic %+.10e  fd %+.10e\n", what, analytic, fd);
    EXPECT_NEAR(analytic, fd, 1e-7 * std::max(1.0, std::abs(fd))) << what;
  };
  for (const cplx d : {cplx(1, 0), cplx(0, 1)}) {
    check(d.real() ? "Re eps" : "Im eps", std::real(std::conj(bar.eps[0]) * d),
          (L(section(0, e0 + h * d, bg0)) - L(section(0, e0 - h * d, bg0))) / (2 * h));
    check(d.real() ? "Re bg" : "Im bg", std::real(std::conj(bar.epsBg) * d),
          (L(section(0, e0, bg0 + h * d)) - L(section(0, e0, bg0 - h * d))) / (2 * h));
  }
  const auto sp = section(h, e0, bg0), sm = section(-h, e0, bg0);
  const double shape = cyl::coupledQShape(bar, map, section(0, e0, bg0), {velocity(sp.regions[0], sm.regions[0], h)});
  check("vertex", shape, (L(sp) - L(sm)) / (2 * h));
}

TEST(CylCoupled, SourceReverseRuleAgainstDifferences) {
  // dL = Re sum conj(Bbar) dB of the class source (coupledSource) for fixed cotangents of Bpm and By: permittivities
  // (real and imaginary parts, background) and a vertex of a pentagon (P = 1, harmonics -1, 0, 1, oblique kpar, Psi)
  // against central differences
  const auto grid = fiberGrid(2.0, 1.6);
  const Eigen::Index nf = 3 * grid.k.size();
  const cyl::HarmonicMap map{0, 1, 1};
  const double kpar = 0.7, Psi = 0.4;
  auto section = [](double dz, cplx e, cplx bg) {
    cyl::AngularSection s;
    s.regions = {region(pentagon(0.0, dz))};
    s.eps = {e};
    s.epsBg = bg;
    return s;
  };
  const cplx e0(4.0, 0.3), bg0(1.2, 0.05);
  fg::MatC Bpmb(2 * nf, 2);
  fg::VecC Byb(nf);
  for (Eigen::Index i = 0; i < 2 * nf; ++i)
    for (Eigen::Index c = 0; c < 2; ++c) Bpmb(i, c) = cplx(std::sin(0.13 * i + c), std::cos(0.29 * i - 0.5 * c));
  for (Eigen::Index i = 0; i < nf; ++i) Byb(i) = cplx(std::cos(0.17 * i), std::sin(0.41 * i + 0.2));
  auto L = [&](const cyl::AngularSection& s) {
    const auto B = cyl::coupledSource(grid, map, s, kpar, Psi);
    return (Bpmb.conjugate().cwiseProduct(B.Bpm).sum() + Byb.conjugate().cwiseProduct(B.By).sum()).real();
  };
  const auto bar = cyl::coupledSourceBar(grid, map, section(0, e0, bg0), kpar, Psi, Bpmb, Byb);
  const double h = 1e-6;
  auto check = [](const char* what, double analytic, double fd) {
    std::printf("coupledSource reverse %-8s analytic %+.10e  fd %+.10e\n", what, analytic, fd);
    EXPECT_NEAR(analytic, fd, 1e-7 * std::max(1.0, std::abs(fd))) << what;
  };
  for (const cplx d : {cplx(1, 0), cplx(0, 1)}) {
    check(d.real() ? "Re eps" : "Im eps", std::real(std::conj(bar.eps[0]) * d),
          (L(section(0, e0 + h * d, bg0)) - L(section(0, e0 - h * d, bg0))) / (2 * h));
    check(d.real() ? "Re bg" : "Im bg", std::real(std::conj(bar.epsBg) * d),
          (L(section(0, e0, bg0 + h * d)) - L(section(0, e0, bg0 - h * d))) / (2 * h));
  }
  const auto sp = section(h, e0, bg0), sm = section(-h, e0, bg0);
  const double shape = cyl::coupledQShape(bar, map, section(0, e0, bg0), {velocity(sp.regions[0], sm.regions[0], h)});
  check("vertex", shape, (L(sp) - L(sm)) / (2 * h));
}
