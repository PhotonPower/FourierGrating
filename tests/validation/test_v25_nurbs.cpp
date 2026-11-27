// V25b-g (Plan section "Erweiterung ab Version 2.37: NURBS-Profil in 1D-Gittern (M24)", "Tests und Abnahme"):
// profile1d "nurbs" against profile1d points/triangle/Function, an overhang (energy, reciprocity, RCWA -> RK,
// truncation), the seam and the start of the curve, Chebyshev and ASR, near field, the checks NURBS_*, gradient
// and sweep behaviour.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <nlohmann/json.hpp>
#include <numbers>

#include "fg/adjoint/gradient.hpp"
#include "fg/core/error.hpp"
#include "fg/io/config.hpp"
#include "fg/solvers/simulation.hpp"
#include "fg/solvers/solver.hpp"

using json = nlohmann::json;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double kD = 0.3, kH = 0.1;  // period and line height (as V24)

fg::Result solve(const json& cfg) {
  const auto c = fg::io::parseConfig(cfg.dump());
  return fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
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

double energy(const fg::Result& r) {
  double v = 0;
  for (const auto& o : r.reflected) v += o.efficiency;
  for (const auto& o : r.transmitted) v += o.efficiency;
  return v;
}

/// 1D stack with the given layers (as V24: period 0.3, conical incidence).
json base(const json& layers, const std::string& solver = "rcwa", int N = 6) {
  return json{{"dimension", "1D"},
              {"solver", solver},
              {"truncation", {{"Nx", N}}},
              {"period_um", {{"x", kD}}},
              {"incidence", {{"lambda_um", 0.633}, {"phi_deg", 20}, {"psi_deg", 30}, {"delta_deg", 35}}},
              {"rk", {{"tol", 1e-10}}},
              {"rcwa", {{"slices", 40}}},
              {"materials",
               {{"air", {{"n", 1.0}}},
                {"glass", {{"n", 1.5}}},
                {"pillar", {{"n", 2.0}}},
                {"crystal", {{"eps_principal", {4.98, 4.98, 4.6}}, {"euler_deg", {20, 70, 10}}}}}},
              {"superstrate", "air"},
              {"substrate", "glass"},
              {"layers", layers}};
}

using Pts = std::vector<std::array<double, 2>>;

/// Clamped knot vector of a chain of Bezier pieces of degree p (interior knots with multiplicity p).
std::vector<double> chainKnots(int pieces, int p) {
  std::vector<double> U(static_cast<std::size_t>(p + 1), 0.0);
  for (int i = 1; i < pieces; ++i) U.insert(U.end(), static_cast<std::size_t>(p), static_cast<double>(i) / pieces);
  U.insert(U.end(), static_cast<std::size_t>(p + 1), 1.0);
  return U;
}

json nurbs(int degree, const Pts& pts, const std::vector<double>& weights = {}, const std::string& below = "pillar",
           const std::string& fac = "laurent") {
  const int pieces = static_cast<int>(pts.size() - 1) / degree;
  json n{{"degree", degree}, {"knots", chainKnots(pieces, degree)}, {"points_um", pts}};
  if (!weights.empty()) n["weights"] = weights;
  return json{{"type", "profile1d"}, {"shape", "nurbs"}, {"below", below}, {"above", "air"}, {"nurbs", n},
              {"factorization", fac}};
}

/// A chain of quadratic Bezier pieces, each {P0, P1, P2} (P2 = next P0), as points of a degree-2 NURBS.
Pts chain(const std::vector<Pts>& pieces) {
  Pts out{pieces.front().front()};
  for (const auto& p : pieces) out.insert(out.end(), p.begin() + 1, p.end());
  return out;
}

/// The same chain starting at piece k (the pieces before it moved by one period).
Pts rotated(const std::vector<Pts>& pieces, std::size_t k) {
  std::vector<Pts> r(pieces.begin() + static_cast<long>(k), pieces.end());
  for (std::size_t i = 0; i < k; ++i) {
    Pts p = pieces[i];
    for (auto& q : p) q[0] += kD;
    r.push_back(p);
  }
  return chain(r);
}

/// Smooth overhang ("mushroom"): stem x ~ 0.13 ... 0.17, cap overhanging to x ~ 0.08 / 0.22, quadratic pieces.
const std::vector<Pts> kMushroom{{{0.0, 0.02}, {0.065, 0.02}, {0.13, 0.02}},   {{0.13, 0.02}, {0.14, 0.06}, {0.12, 0.09}},
                                 {{0.12, 0.09}, {0.05, 0.10}, {0.08, 0.13}},  {{0.08, 0.13}, {0.15, 0.19}, {0.22, 0.13}},
                                 {{0.22, 0.13}, {0.25, 0.10}, {0.18, 0.09}},  {{0.18, 0.09}, {0.16, 0.06}, {0.17, 0.02}},
                                 {{0.17, 0.02}, {0.235, 0.02}, {0.3, 0.02}}};

template <class F>
void expectError(F f, const std::string& part, const std::string& part2 = "") {
  try {
    f();
    ADD_FAILURE() << "no error, expected " << part;
  } catch (const fg::Error& e) {
    const std::string w = e.what();
    EXPECT_NE(w.find(part), std::string::npos) << w;
    if (!part2.empty()) EXPECT_NE(w.find(part2), std::string::npos) << w;
  }
}

/// Rounded trapezoid of V24c (CD 100 nm at half height, h 100 nm, 15 degree walls, r = 10 nm).
struct RoundedLine {
  double xc = 0.15, cd = 0.1, h = kH, r = 0.01, a = 15 * pi / 180;
  double s() const { return std::sin(a); }
  double c() const { return std::cos(a); }
  double xf() const { return xc - 0.5 * cd - 0.5 * h * std::tan(a); }
  double cfx() const { return xf() - r * (1 - s()) / c(); }
  double ctx() const { return xf() + (s() * (h - r) + r) / c(); }
  double g(double x) const {  // y = g(x), as V24c
    const double u = xc - std::abs(x - xc);
    if (u <= cfx()) return 0.0;
    if (u <= cfx() + r * c()) return r - std::sqrt(std::max(0.0, r * r - (u - cfx()) * (u - cfx())));
    if (u <= ctx() - r * c()) return (u - xf()) / std::tan(a);
    if (u <= ctx()) return h - r + std::sqrt(std::max(0.0, r * r - (u - ctx()) * (u - ctx())));
    return h;
  }
  /// Exact degree-2 NURBS: flat, foot arc, wall, head arc, top, head arc, wall, foot arc, flat.
  void nurbs(Pts& pts, std::vector<double>& w) const {
    std::vector<Pts> pieces;
    std::vector<double> mid;
    auto line = [&](std::array<double, 2> p, std::array<double, 2> q) {
      pieces.push_back({p, {0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])}, q});
      mid.push_back(1.0);
    };
    auto arc = [&](double cx, double cy, double t0, double t1) {  // angles in degrees, |t1 - t0| < 180
      const double a0 = t0 * pi / 180, a1 = t1 * pi / 180, m = 0.5 * (a0 + a1), hh = 0.5 * (a1 - a0);
      const double l = r / std::cos(hh);
      pieces.push_back({{cx + r * std::cos(a0), cy + r * std::sin(a0)}, {cx + l * std::cos(m), cy + l * std::sin(m)},
                        {cx + r * std::cos(a1), cy + r * std::sin(a1)}});
      mid.push_back(std::cos(hh));
    };
    const double m2 = 2 * xc, deg = a * 180 / pi;
    auto at = [&](double cx, double cy, double t) {
      return std::array<double, 2>{cx + r * std::cos(t * pi / 180), cy + r * std::sin(t * pi / 180)};
    };
    line({0.0, 0.0}, {cfx(), 0.0});
    arc(cfx(), r, -90, -deg);                                   // foot (concave)
    line(at(cfx(), r, -deg), at(ctx(), h - r, 180 - deg));      // left wall
    arc(ctx(), h - r, 180 - deg, 90);                           // head (convex)
    line({ctx(), h}, {m2 - ctx(), h});                          // top
    arc(m2 - ctx(), h - r, 90, deg);
    line(at(m2 - ctx(), h - r, deg), at(m2 - cfx(), r, 180 + deg));  // right wall
    arc(m2 - cfx(), r, 180 + deg, 270);
    line({m2 - cfx(), 0.0}, {kD, 0.0});
    pts = chain(pieces);
    w.assign(pts.size(), 1.0);
    for (std::size_t i = 0; i < mid.size(); ++i) w[2 * i + 1] = mid[i];
  }
};

}  // namespace

TEST(V25b, Degree1EqualsProfile1DPointsAndTriangle) {
  // trapezoid of V24b: profile1d points (periodic polygon) and the same polyline as a degree-1 NURBS; Laurent: Q is
  // the same, the efficiencies too (RCWA 1e-12, RK 1e-10); isotropic and anisotropic line
  const double t = kH * std::tan(15 * pi / 180);
  const Pts trap{{0.0, 0.0}, {0.0866, 0.0}, {0.0866 + t, kH}, {0.2134 - t, kH}, {0.2134, 0.0}, {kD, 0.0}};
  for (const std::string mat : {"pillar", "crystal"}) {
    const json prof = {{"type", "profile1d"}, {"shape", "points"}, {"depth_um", kH},
                       {"x_um", {0.0, 0.0866, 0.0866 + t, 0.2134 - t, 0.2134}}, {"y_um", {0.0, 0.0, kH, kH, 0.0}},
                       {"below", mat}, {"above", "air"}, {"factorization", "laurent"}};
    const json nl = nurbs(1, trap, {}, mat);
    const double rc = maxDiff(solve(base(json::array({nl}))), solve(base(json::array({prof})))),
                 rk = maxDiff(solve(base(json::array({nl}), "rk")), solve(base(json::array({prof}), "rk")));
    std::printf("V25b trapezoid %s: RCWA %.2e, RK %.2e\n", mat.c_str(), rc, rk);
    EXPECT_LT(rc, 1e-12) << mat;
    EXPECT_LT(rk, 1e-10) << mat;
  }
  // triangle (apex at 0.12): profile1d "triangle"
  const json tri = {{"type", "profile1d"}, {"shape", "triangle"}, {"depth_um", kH}, {"apex_x_um", 0.12},
                    {"below", "pillar"}, {"above", "air"}, {"factorization", "laurent"}};
  const json ntri = nurbs(1, {{0.0, 0.0}, {0.12, kH}, {kD, 0.0}});
  const double rc = maxDiff(solve(base(json::array({ntri}))), solve(base(json::array({tri})))),
               rk = maxDiff(solve(base(json::array({ntri}), "rk")), solve(base(json::array({tri}), "rk")));
  std::printf("V25b triangle: RCWA %.2e, RK %.2e\n", rc, rk);
  EXPECT_LT(rc, 1e-12);
  EXPECT_LT(rk, 1e-10);
}

// Plan v2.54: a points profile with a horizontal piece inside (0, depth) (a step at 0.05 um): Q jumps there. Like the
// degree-1 NURBS (band boundary of M24), RK cuts at that height and takes the one-sided limits; equal 1e-10.
TEST(V25b, InteriorStepOfPointsEqualsDegree1InRK) {
  const Pts step{{0.0, 0.0}, {0.06, 0.0}, {0.09, 0.05}, {0.15, 0.05}, {0.2, kH}, {0.24, kH}, {kD, 0.0}};
  const json prof = {{"type", "profile1d"}, {"shape", "points"}, {"depth_um", kH},
                     {"x_um", {0.0, 0.06, 0.09, 0.15, 0.2, 0.24}}, {"y_um", {0.0, 0.0, 0.05, 0.05, kH, kH}},
                     {"below", "pillar"}, {"above", "air"}, {"factorization", "laurent"}};
  const json nl = nurbs(1, step);
  const double rc = maxDiff(solve(base(json::array({nl}))), solve(base(json::array({prof})))),
               rk = maxDiff(solve(base(json::array({nl}), "rk")), solve(base(json::array({prof}), "rk")));
  std::printf("V25b interior step: RCWA %.2e, RK %.2e\n", rc, rk);
  EXPECT_LT(rc, 1e-12);
  EXPECT_LT(rk, 1e-10);
}

TEST(V25c, RoundedLineAgainstTheAnalyticProfile) {
  // the rounded trapezoid of V24c exactly as a NURBS (rational quadratic arcs) against profile1d Function with the
  // analytic g(x) (RK, Laurent); band boundaries at the arc/wall joints
  const RoundedLine L;
  Pts pts;
  std::vector<double> w;
  L.nurbs(pts, w);
  const json cfg = base(json::array({nurbs(2, pts, w)}), "rk");
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto* l = std::get_if<fg::Sectioned1DLayer>(&c.stack.layers.front().spec);
  ASSERT_NE(l, nullptr);
  EXPECT_NEAR(l->provider->thickness_um(), L.h, 1e-15);
  ASSERT_EQ(l->yBreaks_um.size(), 2u);  // wall / foot arc, wall / head arc
  EXPECT_NEAR(l->yBreaks_um[0], L.r * (1 - L.s()), 1e-15);
  EXPECT_NEAR(l->yBreaks_um[1], L.h - L.r * (1 - L.s()), 1e-15);
  fg::Stack ref = c.stack;
  ref.layers.clear();
  fg::Profile1DLayer p;
  p.depth_um = L.h;
  p.shape = fg::Profile1DLayer::Shape::Function;
  p.g = [L](double x) { return L.g(x); };
  p.below = "pillar";
  p.above = "air";
  ref.layers.push_back({p, fg::Factorization::Laurent});
  const auto a = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
  const auto b = fg::makeSolver(c.solver)->solve(ref, c.scan.front());
  std::printf("V25c rounded line: |nurbs - profile1d Function| = %.3e\n", maxDiff(a, b));
  EXPECT_LT(maxDiff(a, b), 1e-10);
}

TEST(V25d, OverhangEnergyReciprocityAndConvergence) {
  const json layer = nurbs(2, chain(kMushroom), {}, "pillar", "auto");
  json rc = base(json::array({layer}), "rcwa", 8), rk = base(json::array({layer}), "rk", 8);
  rk["rk"] = {{"tol", 1e-11}};
  const fg::Result b = solve(rk);
  std::printf("V25d overhang: energy RK %.2e\n", energy(b) - 1);
  EXPECT_NEAR(energy(b), 1.0, 1e-8);
  double last = 1e300;
  for (int slices : {20, 80, 320}) {
    rc["rcwa"] = {{"slices", slices}};
    const fg::Result a = solve(rc);
    const double d = maxDiff(a, b);
    std::printf("V25d overhang: RCWA(%d) energy %.2e, |RCWA - RK| = %.3e\n", slices, energy(a) - 1, d);
    EXPECT_NEAR(energy(a), 1.0, 1e-10);
    EXPECT_LT(d, 0.5 * last);
    last = d;
  }
  // start of the curve at another piece (the same geometry): equal up to rounding
  const json shifted = nurbs(2, rotated(kMushroom, 2), {}, "pillar", "auto");
  for (const std::string solver : {"rcwa", "rk"}) {
    const double d = maxDiff(solve(base(json::array({layer}), solver, 8)), solve(base(json::array({shifted}), solver, 8)));
    std::printf("V25d shifted start %s: %.2e\n", solver.c_str(), d);
    EXPECT_LT(d, 1e-12) << solver;
  }
  // reciprocity e_n(alpha0) = e_n(-alpha_n) (as V24d): falls with the truncation
  const double lam = 0.4, th = 30.0, thRev = std::asin(-(std::sin(th * pi / 180) - lam / kD)) * 180 / pi;
  auto effR = [](const fg::Result& x, int n) {
    for (const auto& o : x.reflected)
      if (o.n == n) return o.efficiency;
    return -1.0;
  };
  double dev[2];
  const int Ns[2] = {8, 24};
  for (int k = 0; k < 2; ++k) {
    json c = base(json::array({layer}), "rcwa", Ns[k]);
    c["rcwa"] = {{"slices", 80}};
    c["incidence"] = {{"lambda_um", lam}, {"phi_deg", th}, {"psi_deg", 0}, {"delta_deg", 90}};
    json cr = c;
    cr["incidence"]["phi_deg"] = thRev;
    dev[k] = std::abs(effR(solve(c), -1) - effR(solve(cr), -1));
  }
  std::printf("V25d reciprocity: N = 8 %.3e, N = 24 %.3e\n", dev[0], dev[1]);
  EXPECT_LT(dev[1], 0.5 * dev[0]);
  // truncation: FFF_sym (default) is close to its N = 56 value already at N = 8 (not monotonic, measured
  // 2e-5); Laurent approaches it slowly and monotonically (6.4e-4, 2.3e-4 at N = 8, 24)
  auto r0 = [](const fg::Result& x) {
    for (const auto& o : x.reflected)
      if (o.n == 0) return o.efficiency;
    return -1.0;
  };
  json lau = layer;
  lau["factorization"] = "laurent";
  const double ref = r0(solve(base(json::array({layer}), "rcwa", 56)));
  const double f8 = std::abs(r0(solve(base(json::array({layer}), "rcwa", 8))) - ref),
               l8 = std::abs(r0(solve(base(json::array({lau}), "rcwa", 8))) - ref),
               l24 = std::abs(r0(solve(base(json::array({lau}), "rcwa", 24))) - ref);
  std::printf("V25d truncation (R0 against FFF_sym N = 56): FFF_sym N = 8 %.2e; Laurent N = 8 %.2e, N = 24 %.2e\n", f8,
              l8, l24);
  EXPECT_LT(f8, 0.1 * l8);
  EXPECT_LT(l24, 0.5 * l8);
}

TEST(V25e, SeamSmoothAndKinked) {
  // smooth at the seam (horizontal tangent at both ends) and with a kink at the seam (triangle wave from its
  // valley); each described from x = 0 and from an inner point: the same efficiencies
  const std::vector<Pts> smooth{{{0.0, 0.02}, {0.05, 0.02}, {0.1, 0.06}}, {{0.1, 0.06}, {0.15, 0.1}, {0.2, 0.06}},
                                {{0.2, 0.06}, {0.25, 0.02}, {0.3, 0.02}}};
  const std::vector<Pts> kinked{{{0.0, 0.0}, {0.075, 0.05}, {0.15, 0.1}}, {{0.15, 0.1}, {0.225, 0.05}, {0.3, 0.0}}};
  for (const auto* pieces : {&smooth, &kinked})
    for (const std::string solver : {"rcwa", "rk"}) {
      const double d = maxDiff(solve(base(json::array({nurbs(2, chain(*pieces))}), solver)),
                               solve(base(json::array({nurbs(2, rotated(*pieces, 1))}), solver)));
      std::printf("V25e %s %s: %.2e\n", pieces == &smooth ? "smooth" : "kinked", solver.c_str(), d);
      EXPECT_LT(d, 1e-12) << solver;
    }
}

TEST(V25f, ChebyshevAsrAndNearField) {
  json layer = nurbs(2, chain(kMushroom), {}, "glass", "auto");
  const fg::Result exact = solve(base(json::array({layer}), "rk", 8));
  layer["chebyshev_nodes"] = 16;
  const fg::Result cheb = solve(base(json::array({layer}), "rk", 8));
  std::printf("V25f Chebyshev (16 nodes per band) against exact sections: %.2e\n", maxDiff(cheb, exact));
  EXPECT_LT(maxDiff(cheb, exact), 1e-8);
  // ASR in x (breakpoints at the stem): lossless, and close to the result without it
  layer.erase("chebyshev_nodes");
  json asr = base(json::array({layer}), "rcwa", 16);
  asr["rcwa"] = {{"slices", 80}};
  json plain = asr;
  asr["coordinate_transform"] = {{"x", {{"type", "asr"}, {"eta", 0.5}, {"breakpoints_um", {0.08, 0.22}}}}};
  const fg::Result a = solve(asr), p = solve(plain);
  std::printf("V25f ASR: energy %.2e, |ASR - plain| = %.3e (N = 16)\n", energy(a) - 1, maxDiff(a, p));
  EXPECT_NEAR(energy(a), 1.0, 1e-6);  // ~N^-4 through the conversion of the orders (M17: 1e-7 at N = 40)
  EXPECT_LT(maxDiff(a, p), 2e-2);
  // near field: runs, finite, no warnings (a horizontal cut through the cap)
  json nf = plain;
  nf["nearfield"] = {{"enabled", true}, {"horizontal_y_um", 0.11}, {"nx", 31}, {"ny", 21}, {"y_range_um", {-0.05, 0.25}}};
  const auto c = fg::io::parseConfig(nf.dump());
  const fg::NearFieldResult r = fg::makeSolver(c.solver)->nearField(c.stack, c.scan.front(), c.field);
  ASSERT_FALSE(r.vertical.E[0].empty());
  for (int k = 0; k < 3; ++k)
    for (const auto& v : r.vertical.E[k]) ASSERT_TRUE(std::isfinite(v.real()) && std::isfinite(v.imag()));
  EXPECT_TRUE(r.warnings.empty());
}

TEST(V25g, ChecksGradientAndSweep) {
  const json good = nurbs(2, chain(kMushroom), {}, "glass", "auto");
  auto withLayer = [&](const json& l) {  // the nurbs layer as /layers/1 above a film
    return base(json::array({{{"type", "homogeneous"}, {"thickness_um", 0.05}, {"material", "glass"}}, l}));
  };
  auto parse = [](const json& c) { (void)fg::io::parseConfig(c.dump()); };
  json l = good;
  l["nurbs"]["knots"].erase(0);
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_INPUT]");
  l = good;
  l["nurbs"]["points_um"].back()[1] = 0.03;
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_NOT_PERIODIC]");
  l = good;
  for (auto& p : l["nurbs"]["points_um"]) p[0] = p[0].get<double>() * 1.1;
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_PERIOD]");
  l = good;
  l["depth_um"] = 0.15;  // the cap reaches 0.16
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_RANGE]");
  l = good;
  for (auto& p : l["nurbs"]["points_um"]) p[1] = p[1].get<double>() - 0.03;
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_RANGE]");
  l = nurbs(1, {{0.0, 0.1}, {0.2, 0.1}, {0.2, 0.3}, {0.15, 0.3}, {0.15, 0.0}, {0.3, 0.1}});
  expectError([&] { parse(withLayer(l)); }, "/layers/1: [NURBS_SELF_INTERSECTION]");
  // stacks
  json c2 = withLayer(good);
  c2["dimension"] = "2D";
  c2["period_um"]["z"] = 0.3;
  c2["truncation"]["Nz"] = 2;
  expectError([&] { parse(c2); }, "[NURBS_STACK]", "1D");
  json pml = withLayer(good);
  pml["aperiodic"] = {{"pml", {{"x", {{"width_um", 0.1}}}}}};
  expectError([&] { parse(pml); }, "[NURBS_STACK]", "PML");
  // FFF_rig needs a signed normal field
  l = good;
  l["factorization"] = "FFF_rig";
  expectError([&] { solve(withLayer(l)); }, "FFF_rig");
  // depth_um and the control points are no gradient parameters (shape gradients are not part of M24)
  fg::Objective obj;
  obj.weights.push_back({fg::OrderRef{'R', 0, 0}, 1.0});
  l = good;
  l["depth_um"] = 0.2;
  expectError([&] { (void)fg::gradient(withLayer(l).dump(), obj, {"/layers/1/depth_um"}); }, "depth_um");
  expectError([&] { (void)fg::gradient(withLayer(good).dump(), obj, {"/layers/1/nurbs/points_um/3/1"}); }, "points_um");
  // the material gradient works (the layer is a Sectioned1DLayer)
  const auto g = fg::gradient(withLayer(good).dump(), obj, {"/materials/glass/n"});
  EXPECT_TRUE(std::isfinite(g.grad[0]));
  // a sweep over a control point rebuilds the layer: equal to the new configuration
  fg::Simulation sim(withLayer(good).dump());
  sim.setParameter("/layers/1/nurbs/points_um/3/1", 0.11);
  json moved = withLayer(good);
  moved["layers"][1]["nurbs"]["points_um"][3][1] = 0.11;
  EXPECT_LT(maxDiff(sim.solve(), solve(moved)), 1e-14);
  EXPECT_GT(maxDiff(sim.solve(), solve(withLayer(good))), 1e-6);
}
