// Validation M19 (Plan section "M19: Tests und Abnahme"): gradients of the discrete adjoint against central
// difference quotients of the forward computation, per parameter kind (V3, V9, V13, V15 type structures).
// With Laurent / Li factorizations the normal field plays no role and the agreement is 1e-6; with FFF the
// adjoint keeps the normal field frozen and the deviation from the full difference quotient is tabulated.
#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <nlohmann/json.hpp>

#include "fg/adjoint/gradient.hpp"
#include "fg/io/config.hpp"
#include "fg/geometry/section.hpp"
#include "fg/solvers/simulation.hpp"

using json = nlohmann::json;

namespace {

fg::Objective objective(std::initializer_list<std::tuple<char, int, int, double>> w) {
  fg::Objective o;
  for (const auto& [s, n, m, v] : w) o.weights.push_back({fg::OrderRef{s, n, m}, v});
  return o;
}

/// Central difference of L at the JSON pointer (h relative to the value).
double fd(const json& cfg, const fg::Objective& obj, const std::string& ptr, double rel = 1e-6) {
  const json::json_pointer jp(ptr);
  const double v = cfg.at(jp).get<double>(), h = rel * std::max(1.0, std::abs(v));
  json a = cfg, b = cfg;
  a[jp] = v + h;
  b[jp] = v - h;
  return (fg::gradient(a.dump(), obj, {}).value - fg::gradient(b.dump(), obj, {}).value) / (2 * h);
}

void check(const json& cfg, const fg::Objective& obj, const std::vector<std::string>& params, double tol,
           const char* what) {
  const auto g = fg::gradient(cfg.dump(), obj, params);
  for (std::size_t k = 0; k < params.size(); ++k) {
    const double num = fd(cfg, obj, params[k]);
    std::printf("%s %-40s adjoint % .10e  fd % .10e  rel %.2e\n", what, params[k].c_str(), g.grad[k], num,
                std::abs(g.grad[k] - num) / std::max(1e-12, std::abs(num)));
    EXPECT_NEAR(g.grad[k], num, tol * std::max(1.0, std::abs(num))) << what << " " << params[k];
  }
}

json lamellar1D(const std::string& fac) {
  return json{{"dimension", "1D"},
              {"truncation", {{"Nx", 10}}},
              {"period_um", {{"x", 1.0}}},
              {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 15}, {"delta_deg", 30}}},
              {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.05}}}}, {"film", {{"eps", 2.2}}}}},
              {"superstrate", "air"},
              {"substrate", "glass"},
              {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.15}, {"material", "film"}},
                                      json{{"type", "lamellar1d"},
                                           {"thickness_um", 0.4},
                                           {"breakpoints_um", {0.0, 0.3, 0.75, 1.0}},
                                           {"materials", {"air", "bar", "air"}},
                                           {"factorization", fac}}})}};
}

}  // namespace

TEST(M19_Gradient, ValueEqualsSolve) {
  const json cfg = lamellar1D("li_lamellar");
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
  const auto g = fg::gradient(cfg.dump(), objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}}), {});
  ASSERT_EQ(r.reflected.size(), g.result.reflected.size());
  for (std::size_t k = 0; k < r.reflected.size(); ++k) EXPECT_NEAR(r.reflected[k].efficiency, g.result.reflected[k].efficiency, 1e-13);
  for (std::size_t k = 0; k < r.transmitted.size(); ++k)
    EXPECT_NEAR(r.transmitted[k].efficiency, g.result.transmitted[k].efficiency, 1e-13);
}

// V3 type: lamellae (Li and Laurent), conical; thickness, breakpoints, complex n, eps of a film.
TEST(M19_Gradient, LamellarAllKinds) {
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}, {'R', 0, 0, -0.3}});
  for (const char* fac : {"li_lamellar", "laurent"}) {
    check(lamellar1D(fac), obj,
          {"/layers/1/thickness_um", "/layers/0/thickness_um", "/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2",
           "/materials/bar/n/0", "/materials/bar/n/1", "/materials/film/eps"},
          1e-6, fac);
  }
}

// Exact degeneracy: normal incidence on a symmetric lamellar grating (Plan "Entartung").
TEST(M19_Gradient, SymmetricNormalIncidence) {
  json cfg = lamellar1D("li_lamellar");
  cfg["incidence"] = {{"lambda_um", 0.6328}, {"phi_deg", 0}, {"psi_deg", 0}, {"delta_deg", 90}};
  cfg["layers"][1]["breakpoints_um"] = {0.0, 0.25, 0.75, 1.0};
  const auto obj = objective({{'R', 1, 0, 1.0}, {'T', -1, 0, 1.0}, {'T', 0, 0, 0.2}});
  check(cfg, obj, {"/layers/1/thickness_um", "/materials/bar/n/0", "/layers/1/breakpoints_um/2"}, 1e-6, "symmetric");
}

namespace {
json pillars2D(const std::string& fac, bool contour) {
  json layer = contour ? json{{"type", "contour2d"},
                              {"thickness_um", 0.2},
                              {"background", "air"},
                              {"regions", json::array({json{{"material", "pillar"},
                                                            {"circle", {{"center_um", {0.25, 0.27}}, {"radius_um", 0.14}}}}})},
                              {"factorization", fac}}
                       : json{{"type", "polygon2d"},
                              {"thickness_um", 0.2},
                              {"background", "air"},
                              {"inclusions", json::array({json{{"material", "pillar"},
                                                               {"vertices_um", {{0.1, 0.1}, {0.38, 0.12}, {0.33, 0.4}, {0.12, 0.35}}}}})},
                              {"factorization", fac}};
  return json{{"dimension", "2D"},
              {"truncation", {{"Nx", 3}, {"Nz", 3}}},
              {"period_um", {{"x", 0.5}, {"z", 0.5}}},
              {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 15}, {"psi_deg", 30}, {"delta_deg", 35}}},
              {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"pillar", {{"n", 3.0}}}}},
              {"superstrate", "air"},
              {"substrate", "glass"},
              {"layers", json::array({layer})}};
}
}  // namespace

// V9 / V15 type: polygon2d vertices and contour2d circle (centre, radius), Laurent: exact.
TEST(M19_Gradient, ShapesLaurent) {
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}});
  check(pillars2D("laurent", false), obj,
        {"/layers/0/inclusions/0/vertices_um/1/0", "/layers/0/inclusions/0/vertices_um/2/1", "/layers/0/thickness_um",
         "/materials/pillar/n"},
        1e-6, "polygon");
  check(pillars2D("laurent", true), obj,
        {"/layers/0/regions/0/circle/radius_um", "/layers/0/regions/0/circle/center_um/0", "/materials/pillar/n"}, 1e-6,
        "circle");
}

// FFF_sym: the adjoint keeps the raster normal field of Algorithm NF frozen (Plan "Methode", 2.); the full
// difference quotient also moves it. Measured deviation: radius 1.8e-7 (N = 3), n 2e-9.
TEST(M19_Gradient, ShapesFFFFrozenNormalField) {
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}});
  check(pillars2D("FFF_sym", true), obj, {"/layers/0/regions/0/circle/radius_um", "/materials/pillar/n"}, 1e-5,
        "circle FFF_sym");
}

// V13 type: gold sinus profile, TM (staircase slices fixed), depth and shift. The analytic profile normal
// field is smooth and differentiated with the parameter (frozen, the FFF_rig gradient would be off by an
// order of magnitude: -7.0 instead of -83.0 for the depth); measured 7e-7 (depth), shift ~0 (translation).
TEST(M19_Gradient, SinusProfile) {
  json cfg{{"dimension", "1D"},
           {"truncation", {{"Nx", 8}}},
           {"period_um", {{"x", 1.0}}},
           {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"delta_deg", 90}}},
           {"rcwa", {{"slices", 12}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"gold", {{"n", {0.2, 3.4}}}}}},
           {"superstrate", "air"},
           {"substrate", "gold"},
           {"layers", json::array({json{{"type", "profile1d"},
                                        {"shape", "sinus"},
                                        {"depth_um", 0.25},
                                        {"x0_um", 0.1},
                                        {"below", "gold"},
                                        {"above", "air"},
                                        {"factorization", "laurent"}}})}};
  const auto obj = objective({{'R', -1, 0, 1.0}, {'R', 0, 0, 0.5}});
  check(cfg, obj, {"/layers/0/depth_um", "/layers/0/x0_um"}, 1e-6, "sinus laurent");
  cfg["layers"][0]["factorization"] = "FFF_rig";
  check(cfg, obj, {"/layers/0/depth_um", "/layers/0/x0_um"}, 1e-5, "sinus FFF_rig");
}

// Cost: gradient (all parameters) <= 5 x forward (Plan "Aufwand").
TEST(M19_Gradient, Cost) {
  const json cfg = pillars2D("laurent", true);
  const auto obj = objective({{'R', 0, 0, 1.0}});
  const std::vector<std::string> params{"/layers/0/regions/0/circle/radius_um", "/layers/0/regions/0/circle/center_um/0",
                                        "/layers/0/regions/0/circle/center_um/1", "/materials/pillar/n",
                                        "/layers/0/thickness_um"};
  auto time = [&](const std::vector<std::string>& p) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < 3; ++k) (void)fg::gradient(cfg.dump(), obj, p);
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / 3;
  };
  const double fwd = time({}), grad = time(params);
  std::printf("forward %.4f s, gradient (5 parameters) %.4f s, ratio %.2f\n", fwd, grad, grad / fwd);
  EXPECT_LT(grad / fwd, 5.0);
}

// RK: the adjoint of the recorded step sequence against the difference quotient (the quotient may see other
// adaptive steps: measured <= 9e-7), and RK gradient = RCWA gradient for lamellae (Plan: 1e-8; measured 2e-9).
TEST(M19_Gradient, RungeKutta) {
  json cfg = lamellar1D("li_lamellar");
  cfg["solver"] = "rk";
  cfg["rk"] = {{"tol", 1e-10}};
  cfg["truncation"] = {{"Nx", 6}};
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}});
  const std::vector<std::string> params{"/layers/1/thickness_um", "/layers/1/breakpoints_um/1", "/materials/bar/n/1"};
  check(cfg, obj, params, 1e-6, "rk lamellar");
  const auto grk = fg::gradient(cfg.dump(), obj, params);
  cfg["solver"] = "rcwa";
  const auto grc = fg::gradient(cfg.dump(), obj, params);
  for (std::size_t k = 0; k < params.size(); ++k) {
    std::printf("rk %.12e rcwa %.12e  %s\n", grk.grad[k], grc.grad[k], params[k].c_str());
    EXPECT_NEAR(grk.grad[k], grc.grad[k], 1e-8 * std::max(1.0, std::abs(grc.grad[k]))) << params[k];
  }
  // y-dependent layer (sinus profile, Laurent): stage matrices per node
  json prof{{"dimension", "1D"},
            {"solver", "rk"},
            {"rk", {{"tol", 1e-10}}},
            {"truncation", {{"Nx", 6}}},
            {"period_um", {{"x", 1.0}}},
            {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"delta_deg", 90}}},
            {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", 1.8}}}}},
            {"superstrate", "air"},
            {"substrate", "glass"},
            {"layers", json::array({json{{"type", "profile1d"},
                                         {"shape", "sinus"},
                                         {"depth_um", 0.25},
                                         {"x0_um", 0.1},
                                         {"below", "bar"},
                                         {"above", "air"},
                                         {"factorization", "laurent"}}})}};
  check(prof, objective({{'R', -1, 0, 1.0}, {'T', 1, 0, 0.3}}), {"/layers/0/depth_um", "/materials/bar/n"}, 1e-5,
        "rk profile");
}

// Simulation API (Plan "Schnittstellen"): valueAndGrad with the kept layer models = gradient(); vjp = grad.
TEST(M19_Gradient, SimulationApi) {
  const json cfg = lamellar1D("li_lamellar");
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}});
  const std::vector<std::string> params{"/layers/1/thickness_um", "/materials/bar/n/0"};
  fg::Simulation sim(cfg.dump());
  const auto a = sim.valueAndGrad(obj, params);
  const auto b = fg::gradient(cfg.dump(), obj, params);
  EXPECT_NEAR(a.value, b.value, 1e-14);
  for (std::size_t k = 0; k < params.size(); ++k) EXPECT_NEAR(a.grad[k], b.grad[k], 1e-12 * std::abs(b.grad[k]));
  const auto v = sim.vjp(obj, params);
  for (std::size_t k = 0; k < params.size(); ++k) EXPECT_EQ(v[k], a.grad[k]);
  sim.setParameter("/layers/1/thickness_um", 0.45);  // the models follow set_parameter
  json c2 = cfg;
  c2["layers"][1]["thickness_um"] = 0.45;
  EXPECT_NEAR(sim.valueAndGrad(obj, params).grad[0], fg::gradient(c2.dump(), obj, params).grad[0], 1e-12);
}

// V18 type: sectioned layer (frustum, analytic section provider), the top radius as a StackParameter
// (lateral move of the wall, Plan eq:shape3d); Laurent, so the adjoint equals the difference quotient.
namespace {
class FrustumProvider final : public fg::SectionProvider {
 public:
  FrustumProvider(double r0, double r1, double h) : r0_(r0), r1_(r1), h_(h) {}
  fg::CrossSection section(double y) const override {
    fg::CrossSection cs;
    fg::geometry::Region r;
    r.material = "pillar";
    r.outer = fg::geometry::Contour2D::circle({0.25, 0.25}, r0_ + (r1_ - r0_) * y / h_);
    r.normalize(1e-9);
    cs.regions.push_back(r);
    return cs;
  }
  double thickness_um() const override { return h_; }

 private:
  double r0_, r1_, h_;
};

fg::Stack frustumStack(double r1) {
  fg::Stack st;
  st.dx_um = 0.5;
  st.dz_um = 0.5;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["pillar"] = fg::Material::fromIndex("pillar", 2.0);
  st.superstrate = "air";
  st.substrate = "glass";
  fg::SectionedLayer l;
  l.provider = std::make_shared<FrustumProvider>(0.15, r1, 0.2);
  l.background = "air";
  l.slices = 6;
  st.layers.push_back(fg::LayerEntry{l, fg::Factorization::Laurent});
  return st;
}
}  // namespace

TEST(M19_Gradient, SectionedFrustum) {
  fg::SolverOptions o;
  o.Nx = 3;
  o.Nz = 3;
  const fg::Wavelength wl{0.6328, 15 * fg::PI / 180, 30 * fg::PI / 180, 35 * fg::PI / 180};
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}});
  auto L = [&](double r1) {
    const auto r = fg::makeSolver(o)->solve(frustumStack(r1), wl);
    double v = 0;
    for (const auto& e : r.reflected) v += (e.n == 0 && e.m == 0) ? e.efficiency : 0.0;
    for (const auto& e : r.transmitted) v += (e.n == -1 && e.m == 0) ? 0.7 * e.efficiency : 0.0;
    return v;
  };
  const double r1 = 0.1, h = 1e-6;
  const auto g = fg::gradient(o, frustumStack(r1), wl, obj, {fg::StackParameter{"r1", [&](double d) { return frustumStack(r1 + d); }}});
  const double num = (L(r1 + h) - L(r1 - h)) / (2 * h);
  std::printf("frustum r1: adjoint %.10e fd %.10e value %.10e %.10e\n", g.grad[0], num, g.value, L(r1));
  EXPECT_NEAR(g.value, L(r1), 1e-12);
  EXPECT_NEAR(g.grad[0], num, 1e-6 * std::max(1.0, std::abs(num)));
}

// NURBS region (closed rational quadratic circle, as examples/python/optimize_nurbs_pillar.py): control points
// (r linear in them: exact velocity) and weights (central difference of the curve, O(delta^2)).
TEST(M19_Gradient, NurbsControlPointsAndWeights) {
  const double s2 = std::sqrt(0.5), R = 0.15, cx = 0.25, cz = 0.25;
  json pts = json::array({{cx + R, cz}, {cx + R, cz + R}, {cx, cz + R}, {cx - R, cz + R}, {cx - R, cz},
                          {cx - R, cz - R}, {cx, cz - R}, {cx + R, cz - R}, {cx + R, cz}});
  json cfg = pillars2D("laurent", true);
  cfg["layers"][0]["regions"][0] = {{"material", "pillar"},
                                    {"nurbs", {{"degree", 2},
                                               {"knots", {0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4}},
                                               {"points_um", pts},
                                               {"weights", {1, s2, 1, s2, 1, s2, 1, s2, 1}}}}};
  check(cfg, objective({{'T', -1, 0, 1.0}, {'R', 0, 0, 0.5}}),
        {"/layers/0/regions/0/nurbs/points_um/2/1", "/layers/0/regions/0/nurbs/points_um/5/0",
         "/layers/0/regions/0/nurbs/weights/3"},
        1e-6, "nurbs");
}
