// Validation M19 (Plan section "M19: Tests und Abnahme"): gradients of the discrete adjoint against central
// difference quotients of the forward computation, per parameter kind (V3, V9, V13, V15 type structures).
// With Laurent / Li factorizations the normal field plays no role and the agreement is 1e-6; with FFF the
// raster normal field of Algorithm NF moves with the geometry (difference quotient over one raster cell; the
// forward map is not smooth below it, the trend is compared).
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

/// Against the difference quotient with the absolute step h: |adjoint - fd| <= rtol |fd| + atol (raster normal
/// field of Algorithm NF: the forward map is not smooth below the raster cell, the trend is compared).
void checkStep(const json& cfg, const fg::Objective& obj, const std::vector<std::string>& params, double h, double rtol,
               double atol, const char* what) {
  const auto g = fg::gradient(cfg.dump(), obj, params);
  for (std::size_t k = 0; k < params.size(); ++k) {
    const double v = cfg.at(json::json_pointer(params[k])).get<double>();
    const double num = fd(cfg, obj, params[k], h / std::max(1.0, std::abs(v)));
    std::printf("%s %-40s adjoint % .10e  fd(h=%.0e) % .10e\n", what, params[k].c_str(), g.grad[k], h, num);
    EXPECT_NEAR(g.grad[k], num, rtol * std::abs(num) + atol) << what << " " << params[k];
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

// FFF_sym: the raster normal field of Algorithm NF moves with the geometry (difference quotient over one raster
// cell). Frozen, a translation of the circle gave -0.44 against a trend of ~0 (the field stays behind); the
// forward map jumps by ~1e-7 below the raster cell, so translations and vertices are compared with the
// difference quotient at h = 1e-3..2e-3 um: vertices 4e-4..6e-4, centre |0 - (-7e-4)|; radius and n exact.
TEST(M19_Gradient, ShapesFFFNormalFieldMoves) {
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}});
  check(pillars2D("FFF_sym", true), obj, {"/layers/0/regions/0/circle/radius_um", "/materials/pillar/n"}, 1e-5,
        "circle FFF_sym");
  checkStep(pillars2D("FFF_sym", true), obj, {"/layers/0/regions/0/circle/center_um/0"}, 2e-3, 0.0, 2e-3,
            "circle FFF_sym centre");
  checkStep(pillars2D("FFF_sym", false), obj, {"/layers/0/inclusions/0/vertices_um/1/0", "/layers/0/inclusions/0/vertices_um/2/1"},
            1e-3, 2e-3, 0.0, "polygon FFF_sym");
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

// M19 Phase 2: wavelength, angles and the superstrate/substrate permittivity (through alpha, gamma, k0, the
// half-space bases, the incident vector and the efficiency normalization).
TEST(M19_Gradient, IncidenceAndHalfSpaces) {
  const std::vector<std::string> inc{"/incidence/lambda_um", "/incidence/phi_deg", "/incidence/psi_deg",
                                     "/incidence/delta_deg"};
  // conical lamellae (air in the lamellae and as superstrate), glass substrate
  auto params = inc;
  params.insert(params.end(), {"/materials/air/n", "/materials/glass/n"});
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}, {'T', 1, 0, -0.4}, {'R', 0, 0, -0.3}});
  check(lamellar1D("li_lamellar"), obj, params, 1e-6, "lamellar");
  // lossy substrate (gold sinus profile, TM): complex substrate index
  json gold{{"dimension", "1D"},
            {"truncation", {{"Nx", 8}}},
            {"period_um", {{"x", 1.0}}},
            {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"delta_deg", 90}}},
            {"rcwa", {{"slices", 8}}},
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
  check(gold, objective({{'R', -1, 0, 1.0}, {'R', 0, 0, 0.5}}),
        {"/incidence/lambda_um", "/incidence/phi_deg", "/materials/gold/n/0", "/materials/gold/n/1", "/materials/air/n"},
        1e-6, "gold");
  // 2D pillars (Laurent), conical
  check(pillars2D("laurent", true), objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}}), inc, 1e-6,
        "pillars");
  // RK: the recorded steps scale with k0
  json rk = lamellar1D("li_lamellar");
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-10}};
  rk["truncation"] = {{"Nx", 6}};
  check(rk, objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}}),
        {"/incidence/lambda_um", "/incidence/phi_deg", "/materials/glass/n", "/materials/air/n"}, 1e-6, "rk");
}

// checkpoint_every (Plan "Schnittstellen", Speicher): only every k-th prefix S-matrix is kept and the
// slices are recomputed per segment in the backward pass; the gradient is bit-identical.
TEST(M19_Gradient, CheckpointEveryIsBitIdentical) {
  json prof{{"dimension", "1D"},
            {"truncation", {{"Nx", 6}}},
            {"period_um", {{"x", 1.0}}},
            {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"delta_deg", 90}}},
            {"rcwa", {{"slices", 10}}},
            {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", 1.8}}}, {"film", {{"eps", 2.2}}}}},
            {"superstrate", "air"},
            {"substrate", "glass"},
            {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.15}, {"material", "film"}},
                                    json{{"type", "profile1d"},
                                         {"shape", "sinus"},
                                         {"depth_um", 0.25},
                                         {"x0_um", 0.1},
                                         {"below", "bar"},
                                         {"above", "air"},
                                         {"factorization", "laurent"}}})}};
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 1, 0, 0.3}});
  const std::vector<std::string> params{"/layers/1/depth_um", "/materials/bar/n", "/layers/0/thickness_um",
                                        "/incidence/lambda_um", "/materials/glass/n"};
  for (const char* solver : {"rcwa", "rk"}) {
    prof["solver"] = solver;
    prof["rk"] = {{"tol", 1e-9}};
    prof.erase("checkpoint_every");
    const auto ref = fg::gradient(prof.dump(), obj, params);
    for (int k : {2, 3, 7}) {
      prof["checkpoint_every"] = k;
      const auto g = fg::gradient(prof.dump(), obj, params);
      EXPECT_EQ(g.value, ref.value) << solver << " k = " << k;
      for (std::size_t q = 0; q < params.size(); ++q) EXPECT_EQ(g.grad[q], ref.grad[q]) << solver << " k = " << k << " " << params[q];
    }
  }
}

// M19 Phase 2: anisotropic media. Lamellae (Laurent, FFF_A, FFF_B): tensor entries, breakpoints, an isotropic
// material of the same layer, thickness, wavelength; homogeneous film: principal values and Euler angles;
// sinus profile (FFF_A, FFF_B): depth, shift, Euler angle; 2D pillar (Laurent): radius and tensor entry.
TEST(M19_Gradient, Anisotropic) {
  const json co = {{"eps_tensor",
                    {{{2.6, 0.05}, {0.0, 0.0}, {0.3, -0.02}}, {{0.0, 0.0}, {2.4, 0.04}, {0.0, 0.0}}, {{-0.3, 0.02}, {0.0, 0.0}, {2.9, 0.06}}}}};
  const json lnb = {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}};
  json lam{{"dimension", "1D"},
           {"truncation", {{"Nx", 8}}},
           {"period_um", {{"x", 1.0}}},
           {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 15}, {"delta_deg", 30}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.05}}}}, {"co", co}, {"lnb", lnb}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.12}, {"material", "lnb"}},
                                   json{{"type", "lamellar1d"},
                                        {"thickness_um", 0.3},
                                        {"breakpoints_um", {0.0, 0.3, 0.55, 1.0}},
                                        {"materials", {"bar", "co", "air"}}}})}};
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}, {'R', 0, 0, -0.3}});
  for (const char* fac : {"laurent", "FFF_A", "FFF_B"}) {
    lam["layers"][1]["factorization"] = fac;
    check(lam, obj,
          {"/materials/co/eps_tensor/0/0/0", "/materials/co/eps_tensor/0/2/1", "/materials/co/eps_tensor/2/0/0",
           "/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2", "/materials/bar/n/0", "/layers/1/thickness_um",
           "/materials/lnb/eps_principal/2", "/materials/lnb/euler_deg/1", "/incidence/lambda_um"},
          1e-6, fac);
  }
  json prof{{"dimension", "1D"},
            {"truncation", {{"Nx", 8}}},
            {"period_um", {{"x", 1.0}}},
            {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 10}, {"delta_deg", 45}}},
            {"rcwa", {{"slices", 8}}},
            {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"lnb", lnb}}},
            {"superstrate", "air"},
            {"substrate", "glass"},
            {"layers", json::array({json{{"type", "profile1d"},
                                         {"shape", "sinus"},
                                         {"depth_um", 0.2},
                                         {"x0_um", 0.1},
                                         {"below", "lnb"},
                                         {"above", "air"}}})}};
  for (const char* fac : {"FFF_A", "FFF_B"}) {
    prof["layers"][0]["factorization"] = fac;
    check(prof, objective({{'R', -1, 0, 1.0}, {'T', 1, 0, 0.3}}),
          {"/layers/0/depth_um", "/layers/0/x0_um", "/materials/lnb/euler_deg/0", "/materials/lnb/eps_principal/0"}, 1e-6,
          fac);
  }
  json pil = pillars2D("laurent", true);
  pil["materials"]["pillar"] = co;
  check(pil, objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}}),
        {"/layers/0/regions/0/circle/radius_um", "/layers/0/regions/0/circle/center_um/0", "/materials/pillar/eps_tensor/0/2/0"},
        1e-6, "pillar aniso");
}

// Phase 2: coordinate transformations (ASR, Plan section "M17"): M with A, G, numerical half-space bases in X,
// efficiencies from the physical Rayleigh amplitudes; lamellae, profiles (also anisotropic), 2D pillars, RK.
TEST(M19_Gradient, CoordinateTransform) {
  const auto obj = objective({{'R', -1, 0, 1.0}, {'T', 0, 0, 0.5}, {'R', 0, 0, -0.3}});
  json lam = lamellar1D("li_lamellar");
  lam["coordinate_transform"] = {{"x", {{"type", "asr"}, {"eta", 0.5}, {"breakpoints_um", {0.3, 0.75}}}}};
  {
    const auto c = fg::io::parseConfig(lam.dump());
    const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
    const auto g = fg::gradient(lam.dump(), obj, {});
    ASSERT_EQ(r.reflected.size(), g.result.reflected.size());
    for (std::size_t k = 0; k < r.reflected.size(); ++k)
      EXPECT_NEAR(r.reflected[k].efficiency, g.result.reflected[k].efficiency, 1e-12);
    for (std::size_t k = 0; k < r.transmitted.size(); ++k)
      EXPECT_NEAR(r.transmitted[k].efficiency, g.result.transmitted[k].efficiency, 1e-12);
  }
  for (const char* fac : {"li_lamellar", "laurent"}) {
    lam["layers"][1]["factorization"] = fac;
    check(lam, obj,
          {"/layers/1/thickness_um", "/layers/0/thickness_um", "/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2",
           "/materials/bar/n/0", "/materials/bar/n/1", "/materials/film/eps"},
          1e-6, fac);
  }
  // map breakpoints away from the jumps
  lam["coordinate_transform"]["x"]["breakpoints_um"] = {0.1, 0.6};
  check(lam, obj, {"/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2", "/materials/bar/n/0"}, 1e-6, "asr off-jump");
  // "auto": the map follows the jumps, only material and thickness parameters
  lam["coordinate_transform"]["x"]["breakpoints_um"] = "auto";
  check(lam, obj, {"/layers/1/thickness_um", "/materials/bar/n/1"}, 1e-6, "asr auto");
  EXPECT_THROW(fg::gradient(lam.dump(), obj, {"/layers/1/breakpoints_um/1"}), fg::Error);
  // global parameters: numerical half-space bases (eigenvector adjoint), A, G, conversions at k0 alpha0
  lam["coordinate_transform"]["x"]["breakpoints_um"] = {0.3, 0.75};
  lam["materials"]["glass"] = {{"n", {1.5, 0.01}}};
  check(lam, obj,
        {"/incidence/lambda_um", "/incidence/phi_deg", "/incidence/psi_deg", "/incidence/delta_deg", "/materials/air/n",
         "/materials/glass/n/0", "/materials/glass/n/1"},
        1e-6, "asr global");
  lam["materials"]["glass"] = {{"n", 1.5}};
  // RK
  json rk = lam;
  rk["coordinate_transform"]["x"]["breakpoints_um"] = {0.3, 0.75};
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-10}};
  rk["truncation"] = {{"Nx", 6}};
  check(rk, obj, {"/layers/1/thickness_um", "/layers/1/breakpoints_um/1", "/materials/bar/n/1", "/incidence/lambda_um"}, 1e-6,
        "asr rk");
  // sinus profile (Laurent exact, FFF_rig with the normal field at f(X))
  json prof{{"dimension", "1D"},
            {"truncation", {{"Nx", 8}}},
            {"period_um", {{"x", 1.0}}},
            {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 10}, {"delta_deg", 45}}},
            {"rcwa", {{"slices", 8}}},
            {"coordinate_transform", {{"x", {{"type", "asr"}, {"eta", 0.4}, {"breakpoints_um", {0.2, 0.7}}}}}},
            {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.02}}}}}},
            {"superstrate", "air"},
            {"substrate", "glass"},
            {"layers", json::array({json{{"type", "profile1d"},
                                         {"shape", "sinus"},
                                         {"depth_um", 0.25},
                                         {"x0_um", 0.1},
                                         {"below", "bar"},
                                         {"above", "air"},
                                         {"factorization", "laurent"}}})}};
  const auto pobj = objective({{'R', -1, 0, 1.0}, {'T', 1, 0, 0.3}});
  check(prof, pobj, {"/layers/0/depth_um", "/layers/0/x0_um", "/materials/bar/n/0", "/incidence/lambda_um", "/incidence/phi_deg"},
        1e-6, "asr sinus laurent");
  prof["layers"][0]["factorization"] = "FFF_rig";
  check(prof, pobj, {"/layers/0/depth_um", "/layers/0/x0_um", "/materials/bar/n/0"}, 1e-5, "asr sinus FFF_rig");
  // anisotropic lamella and profile
  const json co = {{"eps_tensor",
                    {{{2.6, 0.05}, {0.0, 0.0}, {0.3, -0.02}}, {{0.0, 0.0}, {2.4, 0.04}, {0.0, 0.0}}, {{-0.3, 0.02}, {0.0, 0.0}, {2.9, 0.06}}}}};
  const json lnb = {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}};
  json alam = lamellar1D("FFF_A");
  alam["truncation"] = {{"Nx", 8}};
  alam["materials"]["co"] = co;
  alam["layers"][1]["materials"] = {"bar", "co", "air"};
  alam["coordinate_transform"] = {{"x", {{"type", "asr"}, {"eta", 0.5}, {"breakpoints_um", {0.3, 0.75}}}}};
  for (const char* fac : {"laurent", "FFF_A", "FFF_B"}) {
    alam["layers"][1]["factorization"] = fac;
    check(alam, obj,
          {"/materials/co/eps_tensor/0/2/0", "/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2",
           "/layers/1/thickness_um"},
          1e-6, (std::string("asr aniso ") + fac).c_str());
  }
  prof["materials"]["lnb"] = lnb;
  prof["layers"][0]["below"] = "lnb";
  for (const char* fac : {"FFF_A", "FFF_B"}) {
    prof["layers"][0]["factorization"] = fac;
    check(prof, pobj, {"/layers/0/depth_um", "/layers/0/x0_um", "/materials/lnb/euler_deg/0"}, 1e-6,
          (std::string("asr aniso sinus ") + fac).c_str());
  }
  // 2D pillars, both axes transformed
  const json ct2 = {{"x", {{"type", "asr"}, {"eta", 0.4}, {"breakpoints_um", {0.11, 0.39}}}},
                    {"z", {{"type", "asr"}, {"eta", 0.3}, {"breakpoints_um", {0.13, 0.41}}}}};
  const auto obj2 = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}});
  json circ = pillars2D("laurent", true);
  circ["coordinate_transform"] = ct2;
  check(circ, obj2,
        {"/layers/0/regions/0/circle/radius_um", "/layers/0/regions/0/circle/center_um/0", "/materials/pillar/n",
         "/incidence/lambda_um", "/incidence/phi_deg", "/incidence/psi_deg", "/materials/air/n", "/materials/glass/n"},
        1e-6, "asr circle");
  json poly = pillars2D("laurent", false);
  poly["coordinate_transform"] = ct2;
  check(poly, obj2, {"/layers/0/inclusions/0/vertices_um/1/0", "/layers/0/inclusions/0/vertices_um/2/1", "/layers/0/thickness_um"},
        1e-6, "asr polygon");
  circ["materials"]["pillar"] = co;
  check(circ, obj2, {"/layers/0/regions/0/circle/radius_um", "/materials/pillar/eps_tensor/0/2/0"}, 1e-6, "asr circle aniso");
}

// Phase 2: anisotropic 2D cross-sections with FFF_B; shapes: indicator part exact, the raster normal field of
// Algorithm NF moves with the geometry (difference quotient over one raster cell), compared with the trend.
TEST(M19_Gradient, AnisotropicFFFB2D) {
  const json co = {{"eps_tensor",
                    {{{2.6, 0.05}, {0.0, 0.0}, {0.3, -0.02}}, {{0.0, 0.0}, {2.4, 0.04}, {0.0, 0.0}}, {{-0.3, 0.02}, {0.0, 0.0}, {2.9, 0.06}}}}};
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}});
  json circ = pillars2D("FFF_B", true);
  circ["materials"]["pillar"] = co;
  check(circ, obj, {"/materials/pillar/eps_tensor/0/2/0", "/materials/pillar/eps_tensor/1/1/1", "/layers/0/thickness_um"},
        1e-6, "FFF_B circle");
  check(circ, obj, {"/layers/0/regions/0/circle/radius_um"}, 1e-5, "FFF_B circle radius");
  checkStep(circ, obj, {"/layers/0/regions/0/circle/center_um/0"}, 2e-3, 0.0, 2e-4, "FFF_B circle centre");
  json poly = pillars2D("FFF_B", false);
  poly["materials"]["pillar"] = co;
  check(poly, obj, {"/materials/pillar/eps_tensor/2/2/0"}, 1e-6, "FFF_B polygon");
  checkStep(poly, obj, {"/layers/0/inclusions/0/vertices_um/1/0", "/layers/0/inclusions/0/vertices_um/2/1"}, 1e-3, 5e-3, 0.0,
            "FFF_B polygon vertices");
}

// Phase 2: profile2d (egg crate, sinus_x), isotropic (FFF_rig, Laurent) and anisotropic (FFF_A, FFF_B): materials,
// tensor entries, Euler angles and the depth (relative height kept, analytic normal field moves).
TEST(M19_Gradient, Profile2D) {
  const json lnb = {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}};
  json cfg{{"dimension", "2D"},
           {"truncation", {{"Nx", 3}, {"Nz", 3}}},
           {"period_um", {{"x", 0.5}, {"z", 0.6}}},
           {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 15}, {"psi_deg", 30}, {"delta_deg", 35}}},
           {"rcwa", {{"slices", 6}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.02}}}}, {"lnb", lnb}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({json{{"type", "profile2d"},
                                        {"shape", "egg_crate"},
                                        {"depth_um", 0.15},
                                        {"x0_um", 0.05},
                                        {"z0_um", 0.1},
                                        {"below", "bar"},
                                        {"above", "air"},
                                        {"factorization", "laurent"}}})}};
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}});
  check(cfg, obj, {"/layers/0/depth_um", "/materials/bar/n/0", "/materials/bar/n/1", "/incidence/lambda_um"}, 1e-6,
        "egg crate laurent");
  cfg["layers"][0]["factorization"] = "FFF_rig";
  check(cfg, obj, {"/layers/0/depth_um", "/materials/bar/n/0"}, 1e-6, "egg crate FFF_rig");
  cfg["layers"][0]["below"] = "lnb";
  for (const char* fac : {"FFF_A", "FFF_B"}) {
    cfg["layers"][0]["factorization"] = fac;
    check(cfg, obj, {"/layers/0/depth_um", "/materials/lnb/eps_principal/2", "/materials/lnb/euler_deg/1"}, 1e-6,
          (std::string("egg crate aniso ") + fac).c_str());
  }
  cfg["layers"][0]["shape"] = "sinus_x";
  cfg["layers"][0]["factorization"] = "FFF_A";
  check(cfg, obj, {"/layers/0/depth_um", "/materials/lnb/euler_deg/0"}, 1e-6, "sinus_x aniso FFF_A");
}

// Phase 2: shifts x0/z0 of profile2d (translation of the slice grids). A single layer only changes the phases, so a
// second, fixed egg crate below makes the efficiencies depend on the relative position.
TEST(M19_Gradient, Profile2DShift) {
  const json lnb = {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}};
  auto egg = [](double x0, double z0, const char* below, const char* above) {
    return json{{"type", "profile2d"}, {"shape", "egg_crate"}, {"depth_um", 0.12}, {"x0_um", x0}, {"z0_um", z0},
                {"below", below},      {"above", above},       {"factorization", "laurent"}};
  };
  json cfg{{"dimension", "2D"},
           {"truncation", {{"Nx", 3}, {"Nz", 3}}},
           {"period_um", {{"x", 0.5}, {"z", 0.6}}},
           {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 15}, {"psi_deg", 30}, {"delta_deg", 35}}},
           {"rcwa", {{"slices", 4}}},
           {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.02}}}}, {"lnb", lnb}}},
           {"superstrate", "air"},
           {"substrate", "glass"},
           {"layers", json::array({egg(0.07, 0.11, "bar", "air"), egg(0.0, 0.0, "glass", "bar")})}};
  const auto obj = objective({{'R', 0, 0, 1.0}, {'T', -1, 0, 0.7}, {'T', 0, -1, 0.4}, {'R', 1, 1, 0.3}});
  const std::vector<std::string> ps{"/layers/0/x0_um", "/layers/0/z0_um", "/layers/1/x0_um"};
  check(cfg, obj, ps, 1e-6, "laurent");
  for (auto& l : cfg["layers"]) l["factorization"] = "FFF_rig";
  check(cfg, obj, ps, 1e-6, "FFF_rig");
  cfg["layers"][0]["below"] = "lnb";
  cfg["layers"][1]["above"] = "lnb";
  for (const char* fac : {"FFF_A", "FFF_B"}) {
    for (auto& l : cfg["layers"]) l["factorization"] = fac;
    check(cfg, obj, ps, 1e-6, (std::string("aniso ") + fac).c_str());
  }
}

// Plan v2.13, Stufe A: aperiodic ridge with PML at normal incidence (total-field formulation): scattered powers up
// and down through the window; the stack and the background stack (film, materials, thickness).
namespace {
json pmlRidge() {
  return json{{"dimension", "1D"},
              {"truncation", {{"Nx", 24}}},
              {"period_um", {{"x", 3.0}}},
              {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 0}, {"psi_deg", 0}, {"delta_deg", 40}}},
              {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"bar", {{"n", {1.8, 0.05}}}}, {"film", {{"eps", 2.2}}}}},
              {"superstrate", "air"},
              {"substrate", "glass"},
              {"layers", json::array({json{{"type", "homogeneous"}, {"thickness_um", 0.12}, {"material", "film"}},
                                      json{{"type", "lamellar1d"},
                                           {"thickness_um", 0.3},
                                           {"breakpoints_um", {0.0, 1.2, 1.75, 3.0}},
                                           {"materials", {"air", "bar", "air"}}}})},
              {"aperiodic", {{"pml", {{"x", {{"width_um", 0.5}, {"beta", 1.5}}}}}, {"excitation", "plane_wave"}}}};
}
fg::Objective scatterObjective(double up, double down) {
  fg::Objective o;
  o.scattering = {{fg::ScatterRef{fg::ScatterRef::Kind::Up}, up}, {fg::ScatterRef{fg::ScatterRef::Kind::Down}, down}};
  return o;
}
}  // namespace

TEST(M19_Gradient, AperiodicNormalIncidence) {
  const json cfg = pmlRidge();
  const auto obj = scatterObjective(1.0, 0.6);
  {  // the value is the solver's
    const auto c = fg::io::parseConfig(cfg.dump());
    const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
    const auto g = fg::gradient(cfg.dump(), obj, {});
    EXPECT_NEAR(g.result.scattering.up, r.scattering.up, 1e-12 * std::abs(r.scattering.up));
    EXPECT_NEAR(g.result.scattering.down, r.scattering.down, 1e-12 * std::abs(r.scattering.down));
    EXPECT_NEAR(g.value, r.scattering.up + 0.6 * r.scattering.down, 1e-12);
  }
  const std::vector<std::string> ps{"/layers/1/thickness_um", "/layers/0/thickness_um", "/layers/1/breakpoints_um/1",
                                    "/layers/1/breakpoints_um/2", "/materials/bar/n/0", "/materials/bar/n/1",
                                    "/materials/film/eps"};
  check(cfg, obj, ps, 1e-6, "pml ridge laurent");
  json li = cfg;
  li["layers"][1]["factorization"] = "li_lamellar";
  check(li, scatterObjective(-0.4, 1.0), ps, 1e-6, "pml ridge li");
  // global parameters: wavelength, polarization, superstrate/substrate (numerical bases with PML, S_inc)
  check(cfg, obj, {"/incidence/lambda_um", "/incidence/delta_deg", "/incidence/psi_deg", "/materials/glass/n",
                   "/materials/air/n"},
        1e-6, "pml ridge global");
}

// Stufe A: RK in the same ridge.
TEST(M19_Gradient, AperiodicRK) {
  json rk = pmlRidge();
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-8}};
  rk["truncation"]["Nx"] = 12;
  check(rk, scatterObjective(1.0, 0.6), {"/layers/1/breakpoints_um/1", "/materials/film/eps"}, 2e-6, "pml ridge rk");
}

// Stufe A: 2D contact hole (round hole in a film, PML in x and z, Laurent): radius, centre, film index (stack and
// background) and thickness.
TEST(M19_Gradient, AperiodicContactHole) {
  const json hole{{"dimension", "2D"},
                  {"truncation", {{"Nx", 4}, {"Nz", 4}}},
                  {"period_um", {{"x", 2.0}, {"z", 2.0}}},
                  {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 0}, {"psi_deg", 0}, {"delta_deg", 30}}},
                  {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"film", {{"n", {2.1, 0.01}}}}}},
                  {"superstrate", "air"},
                  {"substrate", "glass"},
                  {"layers", json::array({json{{"type", "contour2d"},
                                               {"thickness_um", 0.3},
                                               {"background", "film"},
                                               {"regions", json::array({json{{"material", "air"},
                                                                             {"circle", {{"center_um", {1.05, 0.98}}, {"radius_um", 0.25}}}}})},
                                               {"factorization", "laurent"}}})},
                  {"aperiodic", {{"pml", {{"x", {{"width_um", 0.3}, {"beta", 1.0}}}, {"z", {{"width_um", 0.3}, {"beta", 1.0}}}}}}}};
  check(hole, scatterObjective(1.0, 0.5),
        {"/layers/0/regions/0/circle/radius_um", "/layers/0/regions/0/circle/center_um/0", "/materials/film/n/0",
         "/materials/film/n/1", "/layers/0/thickness_um"},
        1e-6, "contact hole");
}

// Stufe A: modal excitation. A slab waveguide (core n = 2) ends in a wider block of another material on a
// structured substrate (the same waveguide continues below): |r_0|^2, |t_0|^2 and a higher mode; the core index
// acts on both half-space bases (eigenvector reverse rule) and on the layer.
TEST(M19_Gradient, AperiodicModal) {
  const json slab = {{"type", "lamellar1d"}, {"breakpoints_um", {0.0, 1.3, 1.7, 3.0}}, {"materials", {"air", "core", "air"}}};
  const json cfg{{"dimension", "1D"},
                 {"truncation", {{"Nx", 20}}},
                 {"period_um", {{"x", 3.0}}},
                 {"incidence", {{"lambda_um", 0.6328}}},
                 {"materials", {{"air", {{"n", 1.0}}}, {"core", {{"n", 2.0}}}, {"block", {{"n", {1.7, 0.01}}}}}},
                 {"superstrate", "air"},
                 {"substrate", "air"},
                 {"superstrate_structure", slab},
                 {"substrate_structure", slab},
                 {"layers", json::array({json{{"type", "lamellar1d"},
                                              {"thickness_um", 0.35},
                                              {"breakpoints_um", {0.0, 1.1, 1.85, 3.0}},
                                              {"materials", {"air", "block", "air"}}}})},
                 {"aperiodic", {{"pml", {{"x", {{"width_um", 0.4}, {"beta", 1.0}}}}}, {"excitation", {{"mode", 0}}}}}};
  fg::Objective obj;
  obj.modes = {{fg::ModeRef{'R', 0}, 1.0}, {fg::ModeRef{'T', 0}, -0.5}, {fg::ModeRef{'R', 1}, 0.3}};
  {  // the powers are the solver's
    const auto c = fg::io::parseConfig(cfg.dump());
    const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
    const auto g = fg::gradient(cfg.dump(), obj, {});
    ASSERT_EQ(g.result.modal.reflected.size(), r.modal.reflected.size());
    ASSERT_EQ(g.result.modal.transmitted.size(), r.modal.transmitted.size());
    EXPECT_NEAR(g.value, r.modal.reflected[0].power - 0.5 * r.modal.transmitted[0].power + 0.3 * r.modal.reflected[1].power,
                1e-12);
  }
  check(cfg, obj, {"/layers/0/thickness_um", "/layers/0/breakpoints_um/1", "/materials/block/n/0", "/materials/block/n/1",
                   "/materials/core/n", "/incidence/lambda_um"},
        1e-6, "modal");
}

// Plan v2.13, Stufe B: oblique incidence on the PML axis (scattered-field formulation, Bloch-free expansion with the
// sources in the Bloch basis); stack parameters (edges, ridge material); the value is the solver's.
TEST(M19_Gradient, AperiodicScatteredField) {
  json cfg = pmlRidge();
  cfg["incidence"]["phi_deg"] = 25;
  const auto obj = scatterObjective(1.0, 0.6);
  {
    const auto c = fg::io::parseConfig(cfg.dump());
    const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
    const auto g = fg::gradient(cfg.dump(), obj, {});
    EXPECT_NEAR(g.result.scattering.up, r.scattering.up, 1e-12 * std::abs(r.scattering.up));
    EXPECT_NEAR(g.result.scattering.down, r.scattering.down, 1e-12 * std::abs(r.scattering.down));
  }
  const std::vector<std::string> ps{"/layers/1/breakpoints_um/1", "/layers/1/breakpoints_um/2", "/materials/bar/n/0",
                                    "/materials/bar/n/1"};
  check(cfg, obj, ps, 1e-6, "scattered ridge laurent");
  cfg["layers"][1]["factorization"] = "li_lamellar";
  check(cfg, scatterObjective(-0.3, 1.0), ps, 1e-6, "scattered ridge li");
  // Stufe B2: parameters acting on the background field (film, thicknesses, wavelength, angles, half spaces)
  check(cfg, obj, {"/materials/film/eps", "/layers/0/thickness_um", "/layers/1/thickness_um", "/incidence/lambda_um",
                   "/incidence/phi_deg", "/incidence/psi_deg", "/incidence/delta_deg", "/materials/glass/n",
                   "/materials/air/n"},
        1e-6, "scattered ridge background");
}

// Stufe B: 2D contact hole, oblique in x and z (Bloch-free on both PML axes); the scattered-field formulation at
// normal incidence ("scattered_field": sources without re-expansion) on the ridge.
TEST(M19_Gradient, AperiodicScatteredField2D) {
  const json hole{{"dimension", "2D"},
                  {"truncation", {{"Nx", 3}, {"Nz", 3}}},
                  {"period_um", {{"x", 2.0}, {"z", 2.0}}},
                  {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", 20}, {"psi_deg", 30}, {"delta_deg", 30}}},
                  {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"film", {{"n", {2.1, 0.01}}}}}},
                  {"superstrate", "air"},
                  {"substrate", "glass"},
                  {"layers", json::array({json{{"type", "contour2d"},
                                               {"thickness_um", 0.3},
                                               {"background", "film"},
                                               {"regions", json::array({json{{"material", "air"},
                                                                             {"circle", {{"center_um", {1.05, 0.98}}, {"radius_um", 0.25}}}}})},
                                               {"factorization", "laurent"}}})},
                  {"aperiodic", {{"pml", {{"x", {{"width_um", 0.3}, {"beta", 1.0}}}, {"z", {{"width_um", 0.3}, {"beta", 1.0}}}}}}}};
  check(hole, scatterObjective(1.0, 0.5),
        {"/layers/0/regions/0/circle/radius_um", "/materials/film/n/0", "/materials/film/n/1", "/layers/0/thickness_um",
         "/incidence/lambda_um", "/incidence/phi_deg"},
        1e-6, "scattered hole 2D");
  json ridge = pmlRidge();
  ridge["scattered_field"] = true;
  check(ridge, scatterObjective(1.0, 0.6),
        {"/layers/1/breakpoints_um/1", "/materials/bar/n/1", "/materials/film/eps", "/layers/0/thickness_um",
         "/incidence/lambda_um"},
        1e-6, "scattered ridge normal");
}
