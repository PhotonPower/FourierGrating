// V16a (Plan section "M13: Tests und Abnahme"): tilted walls without CAD. SectionedLayer with analytic
// section providers (frustum, spherical cap) against Profile2D with the same geometry y = g(x, z), RK,
// same factorization: 1e-6 at rk.tol = 1e-9; energy conservation for lossless media 1e-10.
//
// Factorization: Laurent. With FFF the two layer types use different normal fields by construction
// (Profile2D: profile normal of g, (0, 1, 0) on the flat parts; SectionedLayer: Algorithm NF, normal of
// the nearest wall), so they agree only in the limit N -> infinity; Laurent compares the exact
// coefficients and the y dependence. FFF_sym with 3D normals: energy conservation and RK ~ RCWA.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <memory>

#include "fg/geometry/section.hpp"
#include "fg/io/config.hpp"

namespace {

const double kC = 0.25;  // centre (x, z) of the pillar in the 0.5 x 0.5 cell

/// Frustum: radius R(y) = R0 + (R1 - R0) y / h, outward normal (cos t, -(R1 - R0)/h, sin t)/|.|.
class Frustum final : public fg::SectionProvider {
 public:
  Frustum(double r0, double r1, double h, std::string mat) : r0_(r0), r1_(r1), h_(h), mat_(std::move(mat)) {}
  double radius(double y) const { return r0_ + (r1_ - r0_) * y / h_; }
  fg::CrossSection section(double y) const override {
    fg::CrossSection cs;
    fg::geometry::Region r;
    r.material = mat_;
    r.outer = fg::geometry::Contour2D::circle({kC, kC}, radius(y));
    cs.regions.push_back(std::move(r));
    fg::BoundaryChain c;
    const int S = 720;
    for (int s = 0; s < S; ++s) {
      const double t = 2 * fg::PI * s / S;
      const fg::Vec3D n = fg::Vec3D(std::cos(t), -(r1_ - r0_) / h_, std::sin(t)).normalized();
      c.samples.push_back({kC + radius(y) * std::cos(t), kC + radius(y) * std::sin(t), n, false});
    }
    cs.chains.push_back(std::move(c));
    return cs;
  }
  double thickness_um() const override { return h_; }
  /// The same frustum as a profile: material below y = g(r)
  double g(double x, double z) const {
    const double r = std::hypot(x - kC, z - kC);
    return std::clamp(h_ * (r0_ - r) / (r0_ - r1_), 0.0, h_);
  }

 private:
  double r0_, r1_, h_;
  std::string mat_;
};

/// Spherical cap of radius rho and height hc: section radius sqrt(rho^2 - (y - yc)^2), yc = hc - rho.
class Cap final : public fg::SectionProvider {
 public:
  Cap(double rho, double hc, std::string mat) : rho_(rho), hc_(hc), yc_(hc - rho), mat_(std::move(mat)) {}
  fg::CrossSection section(double y) const override {
    const double r = std::sqrt(rho_ * rho_ - (y - yc_) * (y - yc_));
    fg::CrossSection cs;
    fg::geometry::Region reg;
    reg.material = mat_;
    reg.outer = fg::geometry::Contour2D::circle({kC, kC}, r);
    cs.regions.push_back(std::move(reg));
    fg::BoundaryChain c;
    const int S = 720;
    for (int s = 0; s < S; ++s) {
      const double t = 2 * fg::PI * s / S;
      const fg::Vec3D n = fg::Vec3D(r * std::cos(t), y - yc_, r * std::sin(t)) / rho_;
      c.samples.push_back({kC + r * std::cos(t), kC + r * std::sin(t), n, false});
    }
    cs.chains.push_back(std::move(c));
    return cs;
  }
  double thickness_um() const override { return hc_; }
  double g(double x, double z) const {
    const double r2 = (x - kC) * (x - kC) + (z - kC) * (z - kC);
    return r2 < rho_ * rho_ - yc_ * yc_ ? std::clamp(yc_ + std::sqrt(rho_ * rho_ - r2), 0.0, hc_) : 0.0;
  }

 private:
  double rho_, hc_, yc_;
  std::string mat_;
};

fg::io::RunConfig config(const std::string& solver, int N) {
  return fg::io::parseConfig(R"({"dimension": "2D", "solver": ")" + solver + R"(", "truncation": {"Nx": )" +
                             std::to_string(N) + R"(, "Nz": )" + std::to_string(N) + R"(},
    "period_um": {"x": 0.5, "z": 0.5},
    "incidence": {"lambda_um": 0.633, "phi_deg": 15, "psi_deg": 30, "delta_deg": 35},
    "rk": {"tol": 1e-9}, "rcwa": {"slices": 60},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "pillar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass", "layers": []})");
}

fg::Result solve(fg::io::RunConfig c, fg::LayerSpec spec, fg::Factorization fac) {
  fg::LayerEntry le;
  le.spec = std::move(spec);
  le.factorization = fac;
  c.stack.layers.push_back(std::move(le));
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
}

fg::SectionedLayer sectioned(std::shared_ptr<const fg::SectionProvider> p) {
  fg::SectionedLayer l;
  l.provider = std::move(p);
  l.background = "air";
  return l;
}

template <typename P>
fg::Profile2DLayer profile(std::shared_ptr<const P> p) {
  fg::Profile2DLayer l;
  l.shape = fg::Profile2DLayer::Shape::Function;
  l.depth_um = p->thickness_um();
  l.g = [p](double x, double z) { return p->g(x, z); };
  l.below = "pillar";
  l.above = "air";
  return l;
}

double maxDiff(const fg::Result& a, const fg::Result& b) {
  EXPECT_EQ(a.reflected.size(), b.reflected.size());
  EXPECT_EQ(a.transmitted.size(), b.transmitted.size());
  double d = 0;
  for (std::size_t k = 0; k < std::min(a.reflected.size(), b.reflected.size()); ++k)
    d = std::max(d, std::abs(a.reflected[k].efficiency - b.reflected[k].efficiency));
  for (std::size_t k = 0; k < std::min(a.transmitted.size(), b.transmitted.size()); ++k)
    d = std::max(d, std::abs(a.transmitted[k].efficiency - b.transmitted[k].efficiency));
  return d;
}

}  // namespace

TEST(V16a_SectionedLayer, FrustumEqualsProfile2D) {
  const auto c = config("rk", 2);
  const auto fr = std::make_shared<const Frustum>(0.17, 0.07, 0.2, "pillar");
  const auto a = solve(c, sectioned(fr), fg::Factorization::Laurent);
  const auto b = solve(c, profile(fr), fg::Factorization::Laurent);
  const double d = maxDiff(a, b);
  std::printf("frustum: max |dEff| = %.3e, energy %.3e / %.3e, RK steps %ld / %ld\n", d, a.diag.energyBalance - 1,
              b.diag.energyBalance - 1, a.diag.rkSteps, b.diag.rkSteps);
  EXPECT_LT(d, 1e-6);
  EXPECT_NEAR(a.diag.energyBalance, 1.0, 1e-10);
  EXPECT_NEAR(b.diag.energyBalance, 1.0, 1e-10);
}

TEST(V16a_SectionedLayer, SphericalCapEqualsProfile2D) {
  const auto c = config("rk", 2);
  const auto cap = std::make_shared<const Cap>(0.25, 0.12, "pillar");
  const auto a = solve(c, sectioned(cap), fg::Factorization::Laurent);
  const auto b = solve(c, profile(cap), fg::Factorization::Laurent);
  const double d = maxDiff(a, b);
  std::printf("cap: max |dEff| = %.3e, energy %.3e / %.3e\n", d, a.diag.energyBalance - 1, b.diag.energyBalance - 1);
  EXPECT_LT(d, 1e-6);
  EXPECT_NEAR(a.diag.energyBalance, 1.0, 1e-10);
}

TEST(V16a_SectionedLayer, TiltedNormalsFFF) {
  // FFF_sym with the 3D normals of the frustum wall (n_y != 0): lossless, RK and the RCWA staircase agree
  const auto fr = std::make_shared<const Frustum>(0.17, 0.07, 0.2, "pillar");
  fg::SectionedLayer l = sectioned(fr);
  l.normalField.grid = 64;
  const auto rk = solve(config("rk", 2), l, fg::Factorization::FFF_sym);
  const auto rc = solve(config("rcwa", 2), l, fg::Factorization::FFF_sym);
  const auto la = solve(config("rk", 2), l, fg::Factorization::Laurent);
  std::printf("FFF_sym: RK - RCWA(60) %.3e, FFF - Laurent %.3e\n", maxDiff(rk, rc), maxDiff(rk, la));
  EXPECT_NEAR(rk.diag.energyBalance, 1.0, 1e-10);
  EXPECT_NEAR(rc.diag.energyBalance, 1.0, 1e-10);
  EXPECT_LT(maxDiff(rk, rc), 5e-3);  // staircase with 60 slices
  EXPECT_GT(maxDiff(rk, la), 1e-6);  // the normal field is in effect
  EXPECT_LT(maxDiff(rk, la), 5e-2);
}

TEST(V16a_SectionedLayer, NearFieldFollowsTheSection) {
  // vertical cut z = 0.25 through the frustum: the material width shrinks with y (R0 -> R1)
  auto c = config("rk", 2);
  fg::LayerEntry le;
  const auto fr = std::make_shared<const Frustum>(0.17, 0.07, 0.2, "pillar");
  fg::SectionedLayer l = sectioned(fr);
  l.normalField.grid = 64;
  le.spec = l;
  le.factorization = fg::Factorization::FFF_sym;
  c.stack.layers.push_back(le);
  fg::NearFieldRequest req;
  req.y_min_um = 0.01;
  req.y_max_um = 0.19;
  req.ny = 7;
  req.nx = 100;
  req.vertical_z_um = kC;
  req.horizontal_y_um = 0.1;
  const auto nf = fg::makeSolver(c.solver)->nearField(c.stack, c.scan.front(), req);
  const auto& v = nf.vertical;
  ASSERT_EQ(v.axis2_um.size(), 7u);
  for (std::size_t j = 0; j < v.axis2_um.size(); ++j) {
    int inside = 0;
    for (std::size_t i = 0; i < v.axis1_um.size(); ++i) inside += v.material[j * v.axis1_um.size() + i] > 2.0;
    const double width = 2 * fr->radius(v.axis2_um[j]);
    EXPECT_NEAR(inside * 0.5 / 100, width, 0.011) << "y = " << v.axis2_um[j];
    for (std::size_t i = 0; i < v.axis1_um.size(); ++i) EXPECT_TRUE(std::isfinite(std::abs(v.E[0][j * 100 + i])));
  }
}
