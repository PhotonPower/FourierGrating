// CAD slicer (Plan section "M14: Tests und Abnahme"): STEP files written at test time with
// Nurbforge's native writer (unit micrometre, named solids) -> cad layers -> V16b, V17, V18, two
// materials with a cavity, Chebyshev interpolation, faulty files with their codes.
#include <gtest/gtest.h>

#include "../support/std_bessel.hpp"

#include <cmath>
#include <nlohmann/json.hpp>
#include <numbers>
#include <nurbforge/geometry/primitives.hpp>
#include <nurbforge/geometry/transform.hpp>
#include <nurbforge/io/step.hpp>
#include <nurbforge/modeling/construct.hpp>
#include <nurbforge/topology/boolean.hpp>

#include "fg/adjoint/gradient.hpp"
#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
#include "fg/fourier/contour.hpp"
#include "fg/geometry/section.hpp"
#include "fg/io/config.hpp"

using nlohmann::json;
namespace {

constexpr double pi = std::numbers::pi;

// ------------------------------------------------------------------ bodies (CAD coordinates, um)
nf::TrimmedSurface disk(double cx, double cy, double z, double r) {
  return nf::TrimmedSurface(nf::make_rectangle({cx - r, cy - r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}),
                            nf::circle_loop(0.5, 0.5, 0.5));
}

nf::Shell prism(const std::vector<std::array<double, 2>>& xy, double z0, double z1) {
  std::vector<nf::TrimmedSurface> f;
  double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
  for (const auto& p : xy) {
    x0 = std::min(x0, p[0]);
    y0 = std::min(y0, p[1]);
    x1 = std::max(x1, p[0]);
    y1 = std::max(y1, p[1]);
  }
  for (std::size_t i = 0; i < xy.size(); ++i) {
    const auto& a = xy[i];
    const auto& b = xy[(i + 1) % xy.size()];
    f.emplace_back(nf::extrude(nf::make_line({a[0], a[1], z0}, {b[0], b[1], z0}), {0, 0, z1 - z0}));
  }
  std::vector<std::array<double, 2>> uv;
  for (const auto& p : xy) uv.push_back({(p[0] - x0) / (x1 - x0), (p[1] - y0) / (y1 - y0)});
  for (double z : {z0, z1})
    f.emplace_back(nf::make_rectangle({x0, y0, z}, {x1 - x0, 0, 0}, {0, y1 - y0, 0}),
                   nf::polygon_loop(uv));
  return nf::Shell::sew(std::move(f));
}

nf::Shell cylinder(double cx, double cy, double r, double z0, double z1) {
  return nf::Shell::sew({nf::TrimmedSurface(nf::extrude(nf::make_circle({cx, cy, z0}, {0, 0, 1}, r),
                                                        {0, 0, z1 - z0})),
                         disk(cx, cy, z0, r), disk(cx, cy, z1, r)});
}

nf::Shell frustum(double cx, double cy, double r0, double r1, double z0, double h) {
  const auto side = nf::revolve(nf::make_line({cx + r0, cy, z0}, {cx + r1, cy, z0 + h}),
                                {cx, cy, 0}, {0, 0, 1}, 2 * pi);
  return nf::Shell::sew({nf::TrimmedSurface(side), disk(cx, cy, z0, r0), disk(cx, cy, z0 + h, r1)});
}

nf::Shell cap(double cx, double cy, double rho, double hc, double z0) {
  const double cz = z0 + hc - rho, t0 = std::asin((z0 - cz) / rho);
  const auto side = nf::revolve(nf::make_arc({cx, cy, cz}, {1, 0, 0}, {0, 0, 1}, rho, t0, 0.5 * pi),
                                {cx, cy, 0}, {0, 0, 1}, 2 * pi);
  return nf::Shell::sew({nf::TrimmedSurface(side), disk(cx, cy, z0, rho * std::cos(t0))});
}

/// Line of a 1D grating as a CAD body: profile (X, Z) extruded along CAD Y over [0, L]; the lids
/// are trimmed out of rectangles enlarged by pad (carrier surfaces larger than the body).
nf::Shell line(const std::vector<std::array<double, 2>>& xz, double L, double pad = 0) {
  std::vector<nf::TrimmedSurface> f;
  double x0 = 1e300, z0 = 1e300, x1 = -1e300, z1 = -1e300;
  for (const auto& p : xz) {
    x0 = std::min(x0, p[0]);
    z0 = std::min(z0, p[1]);
    x1 = std::max(x1, p[0]);
    z1 = std::max(z1, p[1]);
  }
  for (std::size_t i = 0; i < xz.size(); ++i) {
    const auto& a = xz[i];
    const auto& b = xz[(i + 1) % xz.size()];
    f.emplace_back(nf::extrude(nf::make_line({a[0], 0, a[1]}, {b[0], 0, b[1]}), {0, L, 0}));
  }
  const double X0 = x0 - pad, Z0 = z0 - pad, W = x1 - x0 + 2 * pad, H = z1 - z0 + 2 * pad;
  std::vector<std::array<double, 2>> uv;
  for (const auto& p : xz) uv.push_back({(p[0] - X0) / W, (p[1] - Z0) / H});
  for (double y : {0.0, L})
    f.emplace_back(nf::make_rectangle({X0, y, Z0}, {W, 0, 0}, {0, 0, H}), nf::polygon_loop(uv));
  return nf::Shell::sew(std::move(f));
}

/// Trapezoidal line (CD 0.1 at half height, h = 0.1, 15 degree walls) in a 0.5 um cell.
const std::vector<std::array<double, 2>> kLine{
    {0.1866, 0.0}, {0.3134, 0.0}, {0.2866, 0.1}, {0.2134, 0.1}};

std::string step(const std::vector<nf::Shell>& s, const std::vector<std::string>& names,
                 nf::StepLengthUnit unit = nf::StepLengthUnit::Micrometre) {
  nf::StepOptions o;
  o.unit = unit;
  o.shell_names = names;
  o.timestamp = "2026-11-05T00:00:00";
  return nf::to_step(s, o);
}

// ------------------------------------------------------------------ configurations
/// 2D configuration (period 0.5 x 0.5) with the given layer list; FG (x, z) = (X, 0.5 - Y) for
/// cell_origin_um [0, -0.5].
json base(const json& layers, const std::string& solver = "rcwa", int N = 2) {
  return json{
      {"dimension", "2D"},
      {"solver", solver},
      {"truncation", {{"Nx", N}, {"Nz", N}}},
      {"period_um", {{"x", 0.5}, {"z", 0.5}}},
      {"incidence", {{"lambda_um", 0.633}, {"phi_deg", 15}, {"psi_deg", 30}, {"delta_deg", 35}}},
      {"rk", {{"tol", 1e-9}}},
      {"rcwa", {{"slices", 40}}},
      {"materials",
       {{"air", {{"n", 1.0}}},
        {"glass", {{"n", 1.5}}},
        {"pillar", {{"n", 2.0}}},
        {"lid", {{"n", 1.46}}}}},
      {"superstrate", "air"},
      {"substrate", "glass"},
      {"layers", layers}};
}

json cadLayer(const std::string& stepText, const json& bodies, const std::string& fac = "laurent") {
  return json{{"type", "cad"},
              {"step", stepText},
              {"cell_origin_um", {0.0, -0.5}},
              {"background", "superstrate"},
              {"bodies", bodies},
              {"factorization", fac},
              {"normal_field", {{"grid", 64}}}};
}

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

json report(const json& cfg) {
  return json::parse(fg::io::cadReport(cfg.dump()));
}

std::vector<std::string> codes(const json& r, const char* kind = "errors") {
  std::vector<std::string> c;
  for (const auto& e : r[kind]) c.push_back(e["code"]);
  return c;
}

const std::vector<std::array<double, 2>> kL{{0.08, 0.06}, {0.38, 0.06}, {0.38, 0.18},
                                            {0.20, 0.18}, {0.20, 0.42}, {0.08, 0.42}};

}  // namespace

// ------------------------------------------------------------------ V16b
TEST(CAD_V16b, LPrismIsPolygon2D) {
  const std::string s = step({prism(kL, 0.0, 0.2)}, {"L"});
  const json cad = cadLayer(s, json::array({{{"name", "L"}, {"material", "pillar"}}}), "FFF_sym");
  const auto c = fg::io::parseConfig(base(json::array({cad})).dump());
  ASSERT_EQ(c.stack.layers.size(), 1u);
  const auto* cl = std::get_if<fg::Contour2DLayer>(&c.stack.layers[0].spec);
  ASSERT_NE(cl, nullptr) << "an extruded prism is one y-invariant band";
  EXPECT_NEAR(cl->thickness_um, 0.2, 1e-14);
  ASSERT_EQ(cl->regions.size(), 1u);
  // coefficients = Polygon2D with the hand-mapped polygon (x, z) = (X, 0.5 - Y)
  std::vector<std::array<double, 2>> poly;
  for (const auto& p : kL) poly.push_back({p[0], 0.5 - p[1]});
  auto reg = cl->regions[0];
  reg.normalize(1e-9);
  const auto g = fg::fourier::regionCoefficients(reg, 0.5, 0.5, 8, 8);
  double d = 0;
  for (int p = -8; p <= 8; ++p)
    for (int q = -8; q <= 8; ++q)
      d = std::max(d, std::abs(g(p, q) - fg::fourier::polygonCoefficient(
                                             fg::geometry::Contour2D::polygon(poly).signedArea() > 0
                                                 ? poly
                                                 : std::vector(poly.rbegin(), poly.rend()),
                                             0.5, 0.5, p, q)));
  EXPECT_LT(d, 1e-13);
  // efficiencies: RK and RCWA, cad layer = polygon2d layer
  json pl{{"type", "polygon2d"},
          {"thickness_um", 0.2},
          {"background", "air"},
          {"inclusions", json::array({{{"material", "pillar"}, {"vertices_um", poly}}})},
          {"factorization", "FFF_sym"},
          {"normal_field", {{"grid", 64}}}};
  for (const char* solver : {"rcwa", "rk"}) {
    const auto a = solve(base(json::array({cad}), solver)),
               b = solve(base(json::array({pl}), solver));
    EXPECT_LT(maxDiff(a, b), 1e-12) << solver;
  }
  // mirrored (left-handed) axes are rejected
  json mirrored = cad;
  mirrored["axes"] = {{"x", "+X"}, {"y", "+Z"}, {"z", "+Y"}};
  EXPECT_THROW(fg::io::parseConfig(base(json::array({mirrored})).dump()), fg::Error);
  mirrored["axes"] = {
      {"x", "+X"}, {"y", "+Z"}, {"z", "-Y"}};  // the right-handed choice is accepted
  EXPECT_NO_THROW(fg::io::parseConfig(base(json::array({mirrored})).dump()));
}

// ------------------------------------------------------------------ V17
TEST(CAD_V17, CylinderIsExactCircle) {
  FG_SKIP_WITHOUT_STD_BESSEL();
  const std::string s = step({cylinder(0.25, 0.25, 0.15, 0.0, 0.2)}, {"Pillar"});
  const auto c = fg::io::parseConfig(
      base(json::array({cadLayer(s, json::array({{{"name", "Pillar"}, {"material", "pillar"}}}))}))
          .dump());
  ASSERT_EQ(c.stack.layers.size(), 1u);
  const auto* cl = std::get_if<fg::Contour2DLayer>(&c.stack.layers[0].spec);
  ASSERT_NE(cl, nullptr);
  auto reg = cl->regions.at(0);
  reg.normalize(1e-9);
  const auto g = fg::fourier::regionCoefficients(reg, 0.5, 0.5, 16, 16);
  double d = 0;
  for (int p = -16; p <= 16; ++p)
    for (int q = -16; q <= 16; ++q) {
      const double kx = 4 * pi * p, kz = 4 * pi * q, k = std::hypot(kx, kz), R = 0.15;
      const fg::cplx ref = (p == 0 && q == 0)
                               ? fg::cplx(pi * R * R / 0.25)
                               : 2 * pi * R * fgtest::cyl_bessel_j(1.0, k * R) / (k * 0.25) *
                                     std::exp(-fg::I * (kx * 0.25 + kz * 0.25));  // Plan eq:circle
      d = std::max(d, std::abs(g(p, q) - ref));
    }
  EXPECT_LT(d, 1e-12);
}

// ------------------------------------------------------------------ V18
namespace {
/// The M13 analytic frustum provider in FG coordinates (centre (0.25, 0.25)).
class FrustumProvider final : public fg::SectionProvider {
 public:
  FrustumProvider(double r0, double r1, double h) : r0_(r0), r1_(r1), h_(h) {}
  fg::CrossSection section(double y) const override {
    const double R = r0_ + (r1_ - r0_) * y / h_;
    fg::CrossSection cs;
    fg::geometry::Region r;
    r.material = "pillar";
    r.outer = fg::geometry::Contour2D::circle({0.25, 0.25}, R);
    cs.regions.push_back(r);
    fg::BoundaryChain ch;
    for (int s = 0; s < 720; ++s) {
      const double t = 2 * pi * s / 720;
      ch.samples.push_back({0.25 + R * std::cos(t), 0.25 + R * std::sin(t),
                            fg::Vec3D(std::cos(t), -(r1_ - r0_) / h_, std::sin(t)).normalized(),
                            false});
    }
    cs.chains.push_back(ch);
    return cs;
  }
  double thickness_um() const override { return h_; }

 private:
  double r0_, r1_, h_;
};
}  // namespace

TEST(CAD_V18, FrustumAgainstProfileAndAnalyticProvider) {
  const double r0 = 0.17, r1 = 0.07, h = 0.2;
  const std::string s = step({frustum(0.25, 0.25, r0, r1, 0.0, h)}, {"Frustum"});
  const json bodies = json::array({{{"name", "Frustum"}, {"material", "pillar"}}});
  const auto rep = report(base(json::array({cadLayer(s, bodies)})));
  ASSERT_TRUE(rep["ok"].get<bool>()) << rep.dump(2);
  ASSERT_EQ(rep["bands"].size(), 1u);
  EXPECT_EQ(rep["bands"][0]["kind"], "sectioned");
  const auto cadRK = solve(base(json::array({cadLayer(s, bodies)}), "rk"));
  EXPECT_NEAR(cadRK.diag.energyBalance, 1.0, 1e-10);
  // against the analytic provider (M13) in a hand-built stack
  auto c = fg::io::parseConfig(base(json::array(), "rk").dump());
  fg::SectionedLayer sl;
  sl.provider = std::make_shared<FrustumProvider>(r0, r1, h);
  sl.background = "air";
  sl.normalField.grid = 64;
  fg::LayerEntry le;
  le.spec = sl;
  le.factorization = fg::Factorization::Laurent;
  c.stack.layers.push_back(le);
  const auto ana = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
  EXPECT_LT(maxDiff(cadRK, ana), 1e-9);
  // against Profile2D with g(x, z) (material below g)
  auto cp = fg::io::parseConfig(base(json::array(), "rk").dump());
  fg::Profile2DLayer pl;
  pl.shape = fg::Profile2DLayer::Shape::Function;
  pl.depth_um = h;
  pl.below = "pillar";
  pl.above = "air";
  pl.g = [=](double x, double z) {
    return std::clamp(h * (r0 - std::hypot(x - 0.25, z - 0.25)) / (r0 - r1), 0.0, h);
  };
  fg::LayerEntry lp;
  lp.spec = pl;
  lp.factorization = fg::Factorization::Laurent;
  cp.stack.layers.push_back(lp);
  const auto pro = fg::makeSolver(cp.solver)->solve(cp.stack, cp.scan.front());
  EXPECT_LT(maxDiff(cadRK, pro), 1e-6);
  // RCWA staircase converges with the number of slices
  double prev = 1;
  for (int S : {10, 20, 40}) {
    json l = cadLayer(s, bodies);
    l["slicing"] = {{"rcwa_slices", S}};
    const double e = maxDiff(solve(base(json::array({l}), "rcwa")), cadRK);
    EXPECT_LT(e, prev) << S;
    prev = e;
  }
  EXPECT_LT(prev, 5e-3);
  // Chebyshev interpolation in y = exact sections
  json lc = cadLayer(s, bodies);
  lc["slicing"] = {{"y_interpolation", "chebyshev"}, {"chebyshev_nodes", 16}};
  EXPECT_LT(maxDiff(solve(base(json::array({lc}), "rk")), cadRK), 1e-9);
}

TEST(CAD_V18, SphericalCapAgainstProfile) {
  const double rho = 0.25, hc = 0.12;
  const std::string s = step({cap(0.25, 0.25, rho, hc, 0.0)}, {"Cap"});
  const auto cadRK = solve(base(
      json::array({cadLayer(s, json::array({{{"name", "Cap"}, {"material", "pillar"}}}))}), "rk"));
  EXPECT_NEAR(cadRK.diag.energyBalance, 1.0, 1e-10);
  auto cp = fg::io::parseConfig(base(json::array(), "rk").dump());
  fg::Profile2DLayer pl;
  pl.shape = fg::Profile2DLayer::Shape::Function;
  pl.depth_um = hc;
  pl.below = "pillar";
  pl.above = "air";
  const double yc = hc - rho;
  pl.g = [=](double x, double z) {
    const double r2 = (x - 0.25) * (x - 0.25) + (z - 0.25) * (z - 0.25);
    return r2 < rho * rho - yc * yc ? std::clamp(yc + std::sqrt(rho * rho - r2), 0.0, hc) : 0.0;
  };
  fg::LayerEntry lp;
  lp.spec = pl;
  lp.factorization = fg::Factorization::Laurent;
  cp.stack.layers.push_back(lp);
  EXPECT_LT(maxDiff(cadRK, fg::makeSolver(cp.solver)->solve(cp.stack, cp.scan.front())), 1e-6);
}

// ------------------------------------------------------------------ two materials and a cavity
TEST(CAD, TwoMaterialsAndCavity) {
  // pillar with a lid of another material, above it a block with a spherical cavity
  const nf::Shell block =
      nf::boolean(prism({{0.05, 0.05}, {0.45, 0.05}, {0.45, 0.45}, {0.05, 0.45}}, 0.3, 0.5),
                  nf::Shell::sew({nf::TrimmedSurface(nf::revolve(
                      nf::make_arc({0.25, 0.25, 0.4}, {0, 0, -1}, {1, 0, 0}, 0.06, 0.0, pi),
                      {0.25, 0.25, 0.4}, {0, 0, 1}, 2 * pi))}),
                  nf::BooleanOp::Difference);
  const std::string s =
      step({cylinder(0.25, 0.25, 0.12, 0.0, 0.15), cylinder(0.25, 0.25, 0.12, 0.15, 0.25), block},
           {"Pillar", "Lid", "Block"});
  const json bodies = json::array({{{"name", "Pillar"}, {"material", "pillar"}},
                                   {{"name", "Lid"}, {"material", "lid"}},
                                   {{"name", "Block"}, {"material", "pillar"}}});
  json l = cadLayer(s, bodies);
  l["slicing"] = {{"rcwa_slices", 30}};
  const auto rep = report(base(json::array({l})));
  ASSERT_TRUE(rep["ok"].get<bool>()) << rep.dump(2);
  std::vector<std::string> kinds;
  for (const auto& b : rep["bands"]) kinds.push_back(b["kind"]);
  // 0-0.15 pillar, 0.15-0.25 lid, 0.25-0.3 air, block: below / around / above the cavity
  ASSERT_EQ(kinds, (std::vector<std::string>{"extrusion", "extrusion", "homogeneous", "extrusion",
                                             "sectioned", "extrusion"}));
  EXPECT_EQ(rep["bands"][2]["material"], "air");
  const auto rk = solve(base(json::array({l}), "rk")), rc = solve(base(json::array({l}), "rcwa"));
  EXPECT_NEAR(rk.diag.energyBalance, 1.0, 1e-10);
  EXPECT_NEAR(rc.diag.energyBalance, 1.0, 1e-10);
  EXPECT_LT(maxDiff(rk, rc), 2e-3);  // RCWA staircase (30 slices) in the cavity band
}

// ------------------------------------------------------------------ faulty files
TEST(CAD, FaultyFilesGiveTheirCodes) {
  const json pillar = json::array({{{"name", "P"}, {"material", "pillar"}}});
  const auto codesOf = [&](const std::string& s, const json& bodies,
                           const json& extra = json::object()) {
    json l = cadLayer(s, bodies);
    for (const auto& [k, v] : extra.items()) l[k] = v;
    return codes(report(base(json::array({l}))));
  };
  using V = std::vector<std::string>;
  EXPECT_EQ(codesOf(step({cylinder(0.25, 0.25, 0.1, 0, 0.2)}, {"P"}), pillar), V{});
  // open shell: a box without its top face
  {
    std::vector<nf::TrimmedSurface> f;
    const double a = 0.1, b = 0.4;
    f.emplace_back(nf::make_rectangle({a, a, 0}, {b - a, 0, 0}, {0, b - a, 0}));
    f.emplace_back(nf::make_rectangle({a, a, 0}, {b - a, 0, 0}, {0, 0, 0.2}));
    f.emplace_back(nf::make_rectangle({a, b, 0}, {b - a, 0, 0}, {0, 0, 0.2}));
    f.emplace_back(nf::make_rectangle({a, a, 0}, {0, b - a, 0}, {0, 0, 0.2}));
    f.emplace_back(nf::make_rectangle({b, a, 0}, {0, b - a, 0}, {0, 0, 0.2}));
    EXPECT_EQ(codesOf(step({nf::Shell::sew(std::move(f))}, {"P"}), pillar), V{"SHELL_OPEN"});
  }
  // overlapping bodies
  EXPECT_EQ(codesOf(step({cylinder(0.2, 0.25, 0.1, 0, 0.2), cylinder(0.3, 0.25, 0.1, 0, 0.2)},
                         {"P", "Q"}),
                    json::array({{{"name", "P"}, {"material", "pillar"}},
                                 {{"name", "Q"}, {"material", "lid"}}})),
            V{"BODY_OVERLAP"});
  // a body overlapping its periodic image
  EXPECT_EQ(codesOf(step({cylinder(0.25, 0.25, 0.3, 0, 0.2)}, {"P"}), pillar),
            V{"PERIODIC_OVERLAP"});
  // micrometre coordinates declared as millimetres
  EXPECT_EQ(
      codesOf(step({cylinder(0.25, 0.25, 0.1, 0, 0.2)}, {"P"}, nf::StepLengthUnit::Millimetre),
              pillar),
      V{"UNIT_SUSPECT"});
  // ... which cad.unit = "um" corrects
  EXPECT_EQ(
      codesOf(step({cylinder(0.25, 0.25, 0.1, 0, 0.2)}, {"P"}, nf::StepLengthUnit::Millimetre),
              pillar, {{"unit", "um"}}),
      V{});
  // body without material
  EXPECT_EQ(codesOf(step({cylinder(0.15, 0.25, 0.05, 0, 0.2), cylinder(0.35, 0.25, 0.05, 0, 0.2)},
                         {"P", "Q"}),
                    pillar),
            V{"BODY_UNASSIGNED"});
  // structure outside the y range
  EXPECT_EQ(codesOf(step({cylinder(0.25, 0.25, 0.1, 0, 0.2)}, {"P"}), pillar,
                    {{"y_range_um", {0.0, 0.1}}}),
            V{"OUT_OF_RANGE"});
  // not a STEP file
  EXPECT_EQ(codesOf("ISO-10303-21;\nthis is not STEP", pillar), V{"STEP_UNSUPPORTED"});
}

TEST(CAD, ReportPreviewAndToolSpec) {
  const std::string s = step({cylinder(0.25, 0.25, 0.1, 0, 0.2)}, {"P"});
  const json cfg =
      base(json::array({cadLayer(s, json::array({{{"name", "P"}, {"material", "pillar"}}}))}));
  const auto r = report(cfg);
  ASSERT_EQ(r["bands"].size(), 1u);
  EXPECT_NEAR(r["bands"][0]["area_fractions"]["pillar"].get<double>(), pi * 0.01 / 0.25, 1e-12);
  const std::string svg = fg::io::cadPreview(cfg.dump(), "", "", 0.1);
  EXPECT_NE(svg.find("<svg"), std::string::npos);
  EXPECT_NE(svg.find("<title>pillar</title>"), std::string::npos);
  const json spec = json::parse(fg::cad::llmToolSpec());
  EXPECT_EQ(spec["name"], "fouriergrating_check_cad");
  EXPECT_TRUE(spec["parameters"]["required"].size() == 2);
}

TEST(CAD, FilesFromDisk) {
  // the STEP files of tests/data/cad (scripts/make_cad_examples.py), relative to the configuration
  json l{{"type", "cad"},
         {"file", "l_prism.step"},
         {"cell_origin_um", {0.0, -0.5}},
         {"bodies", json::array({{{"name", "L"}, {"material", "pillar"}}})},
         {"factorization", "FFF_sym"},
         {"normal_field", {{"grid", 64}}}};
  const std::string dir = std::string(FG_SOURCE_DIR) + "/tests/data/cad";
  const auto a = fg::io::parseConfig(base(json::array({l})).dump(), dir);
  ASSERT_EQ(a.stack.layers.size(), 1u);
  const json inline_ = cadLayer(step({prism(kL, 0.0, 0.2)}, {"L"}),
                                json::array({{{"name", "L"}, {"material", "pillar"}}}), "FFF_sym");
  EXPECT_LT(maxDiff(solve(base(json::array({inline_}))),
                    fg::makeSolver(a.solver)->solve(a.stack, a.scan.front())),
            1e-13);
  l["file"] = "cylinder.step";
  l["bodies"] = json::array({{{"name", "Pillar"}, {"material", "pillar"}}});
  const auto r = json::parse(fg::io::cadReport(base(json::array({l})).dump(), dir));
  EXPECT_TRUE(r["ok"].get<bool>()) << r.dump(2);
  EXPECT_EQ(r["bands"][0]["kind"], "extrusion");
}

// ------------------------------------------------------------------ V20g (M18)
// Contact hole as a CAD body (cylinder of air in the film, PML in x and z) = the same hole as Contour2D
// (Plan V20g, 1e-10), normal incidence (total field) and oblique incidence (scattered field).
TEST(CAD_V20g, ContactHoleEqualsContour2D) {
  const std::string s = step({cylinder(1.0, 1.0, 0.25, 0.0, 0.3)}, {"Hole"});
  for (double phi : {0.0, 20.0}) {
    json cfg{{"dimension", "2D"},
             {"truncation", {{"Nx", 4}, {"Nz", 4}}},
             {"period_um", {{"x", 2.0}, {"z", 2.0}}},
             {"incidence", {{"lambda_um", 0.6328}, {"phi_deg", phi}, {"delta_deg", 30}}},
             {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"film", {{"n", 2.1}}}}},
             {"superstrate", "air"},
             {"substrate", "glass"},
             {"layers", json::array({json{{"type", "cad"},
                                          {"step", s},
                                          {"cell_origin_um", {0.0, -2.0}},
                                          {"background", "film"},
                                          {"bodies", json::array({{{"name", "Hole"}, {"material", "air"}}})},
                                          {"factorization", "laurent"}}})},
             {"aperiodic", {{"pml", {{"x", {{"width_um", 0.3}}}, {"z", {{"width_um", 0.3}}}}}}}};
    const auto c = fg::io::parseConfig(cfg.dump());
    ASSERT_EQ(c.stack.layers.size(), 1u);
    const fg::Result a = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
    fg::Stack ref = c.stack;
    fg::Contour2DLayer l;
    l.thickness_um = 0.3;
    l.background = "film";
    fg::geometry::Region r;
    r.outer = fg::geometry::Contour2D::circle({1.0, 1.0}, 0.25);
    r.material = "air";
    r.normalize(1e-9);
    l.regions = {r};
    ref.layers = {fg::LayerEntry{l, fg::Factorization::Laurent}};
    const fg::Result b = fg::makeSolver(c.solver)->solve(ref, c.scan.front());
    ASSERT_TRUE(a.scattering.present && b.scattering.present);
    EXPECT_NEAR(a.scattering.up, b.scattering.up, 1e-10) << phi;
    EXPECT_NEAR(a.scattering.down, b.scattering.down, 1e-10) << phi;
    EXPECT_NEAR(a.scattering.lateral, b.scattering.lateral, 1e-10) << phi;
    EXPECT_NEAR(a.scattering.absorbed, b.scattering.absorbed, 1e-10) << phi;
    EXPECT_GT(a.scattering.total(), 0.1) << phi;
  }
}

// ------------------------------------------------------------------ M19: CAD gradient (V18)
// Lateral move of a CAD wall (frustum top radius, Plan eq:shape3d): the structure at theta +- delta comes
// from the re-sliced STEP body; adjoint = difference quotient of the forward computation (Laurent).
TEST(CAD_M19, FrustumTopRadiusGradient) {
  const double r0 = 0.15, r1 = 0.1, h = 0.2;
  auto cfgFor = [&](double top) {
    const std::string s = step({frustum(0.25, 0.25, r0, top, 0.0, h)}, {"Frustum"});
    return base(json::array({cadLayer(s, json::array({{{"name", "Frustum"}, {"material", "pillar"}}}))}), "rcwa", 3);
  };
  const auto c = fg::io::parseConfig(cfgFor(r1).dump());
  fg::Objective obj;
  obj.weights = {{fg::OrderRef{'R', 0, 0}, 1.0}, {fg::OrderRef{'T', -1, 0}, 0.7}};
  const auto g = fg::gradient(c.solver, c.stack, c.scan.front(), obj,
                              {fg::StackParameter{"r1", [&](double d) { return fg::io::parseConfig(cfgFor(r1 + d).dump()).stack; }}});
  auto L = [&](double top) {
    const auto r = solve(cfgFor(top));
    double v = 0;
    for (const auto& e : r.reflected) v += (e.n == 0 && e.m == 0) ? e.efficiency : 0.0;
    for (const auto& e : r.transmitted) v += (e.n == -1 && e.m == 0) ? 0.7 * e.efficiency : 0.0;
    return v;
  };
  const double dh = 1e-6, num = (L(r1 + dh) - L(r1 - dh)) / (2 * dh);
  std::printf("CAD frustum r1: adjoint %.10e fd %.10e\n", g.grad[0], num);
  EXPECT_NEAR(g.value, L(r1), 1e-12);
  EXPECT_NEAR(g.grad[0], num, 1e-5 * std::max(1.0, std::abs(num)));
}

// ------------------------------------------------------------------ trimmed bounds, 1D, wrapping
TEST(CAD, BoundsOfTrimmedLids) {
  // lids trimmed out of enlarged rectangles: y range and cell origin from the trimmed body
  // (before: from the carrier control points, adding 0.05 um background above and below)
  const json bodies = json::array({{{"name", "Line"}, {"material", "pillar"}}});
  json tight = cadLayer(step({line(kLine, 0.3)}, {"Line"}), bodies),
       padded = cadLayer(step({line(kLine, 0.3, 0.05)}, {"Line"}), bodies);
  tight.erase("cell_origin_um");
  padded.erase("cell_origin_um");
  const json a = report(base(json::array({tight}))), b = report(base(json::array({padded})));
  EXPECT_TRUE(b["ok"].get<bool>());
  ASSERT_EQ(b["bands"].size(), 1u);
  EXPECT_NEAR(b["bands"][0]["y0_um"].get<double>(), 0.0, 1e-12);
  EXPECT_NEAR(b["bands"][0]["y1_um"].get<double>(), 0.1, 1e-12);
  EXPECT_NEAR(b["cell_origin_um"][0].get<double>(), 0.1866, 1e-12);
  EXPECT_NEAR(b["cell_origin_um"][1].get<double>(), -0.3, 1e-12);
  EXPECT_EQ(a["cell_origin_um"], b["cell_origin_um"]);
  ASSERT_EQ(a["bands"].size(), 1u);
  EXPECT_EQ(a["bands"][0]["kind"], b["bands"][0]["kind"]);
  EXPECT_EQ(a["bands"][0]["y1_um"], b["bands"][0]["y1_um"]);
  EXPECT_NEAR(b["bands"][0]["area_fractions"]["pillar"].get<double>(), 0.12, 1e-12);
}

TEST(CAD, ZInvariantLine) {
  const json bodies = json::array({{{"name", "Line"}, {"material", "pillar"}}});
  json cfg = base(json::array({cadLayer(step({line(kLine, 0.5)}, {"Line"}), bodies)}));
  cfg["truncation"] = {{"Nx", 4}, {"Nz", 0}};
  // extruded over the whole z period: touches its image, which check reported as ok and solve
  // rejected ("regions overlap its periodic image")
  EXPECT_EQ(codes(report(cfg)), std::vector<std::string>{"PERIODIC_OVERLAP"});
  cfg["layers"][0]["periodic_union"] = true;  // cannot resolve it either
  EXPECT_EQ(codes(report(cfg)), std::vector<std::string>{"PERIODIC_OVERLAP"});
  EXPECT_THROW((void)solve(cfg), fg::Error);
  // a shorter line (a 2D grating of dashes) is fine
  cfg["layers"][0] = cadLayer(step({line(kLine, 0.3)}, {"Line"}), bodies);
  EXPECT_EQ(codes(report(cfg)), std::vector<std::string>{});
  EXPECT_NO_THROW((void)solve(cfg));
  // the same body in a 1D stack is a line grating (M23): no PERIODIC_OVERLAP, no NaN fractions
  cfg["dimension"] = "1D";
  cfg["period_um"] = {{"x", 0.5}};
  cfg["layers"][0] = cadLayer(step({line(kLine, 0.5)}, {"Line"}), bodies);
  cfg["layers"][0]["cell_origin_um"] = 0.0;
  cfg["layers"][0]["periodic_union"] = true;
  const json r1 = report(cfg);
  EXPECT_EQ(codes(r1), std::vector<std::string>{});
  EXPECT_EQ(r1["bands"][0]["kind"], "sectioned1d");
  EXPECT_NEAR(r1["bands"][0]["area_fractions"]["pillar"].get<double>(), 0.2, 1e-12);
}

// ------------------------------------------------------------------ M23: CAD import for 1D gratings (V24)
namespace {

constexpr double kD = 0.3, kH = 0.1;  // period and line height of V24

/// Closed chain of curves in the CAD X-Z plane (FG x, y): faces extruded along CAD Y over [0, L], lids trimmed
/// by the same curves in (u, v) (first loop outer, the others holes).
nf::Shell profileBody(const std::vector<std::vector<nf::NurbsCurve>>& loops, double L) {
  double X0 = 1e300, Z0 = 1e300, X1 = -1e300, Z1 = -1e300;
  for (const auto& l : loops)
    for (const auto& c : l)
      for (const auto& p : c.control_points()) {
        X0 = std::min(X0, p.x);
        Z0 = std::min(Z0, p.z);
        X1 = std::max(X1, p.x);
        Z1 = std::max(Z1, p.z);
      }
  const double W = X1 - X0, H = Z1 - Z0;
  const nf::Transform uv({1 / W, 0, 0, -X0 / W, 0, 0, 1 / H, -Z0 / H, 0, 0, 0, 0});
  std::vector<nf::TrimmedSurface> f;
  std::vector<nf::TrimLoop> trims;
  for (const auto& l : loops) {
    nf::TrimLoop t;
    for (const auto& c : l) {
      f.emplace_back(nf::extrude(c, {0, L, 0}));
      t.curves.push_back(c.transformed(uv));
    }
    trims.push_back(t);
  }
  const std::vector<nf::TrimLoop> holes(trims.begin() + 1, trims.end());
  for (double y : {0.0, L})
    f.emplace_back(nf::make_rectangle({X0, y, Z0}, {W, 0, 0}, {0, 0, H}), trims.front(), holes);
  return nf::Shell::sew(std::move(f));
}

std::vector<nf::NurbsCurve> polyCurves(const std::vector<std::array<double, 2>>& xz) {
  std::vector<nf::NurbsCurve> c;
  for (std::size_t i = 0; i < xz.size(); ++i) {
    const auto& a = xz[i];
    const auto& b = xz[(i + 1) % xz.size()];
    c.push_back(nf::make_line({a[0], 0, a[1]}, {b[0], 0, b[1]}));
  }
  return c;
}

/// The line of the report (v2.21): trapezoid, CD at half height, 15 degree walls, corners rounded with exact
/// arcs of radius r (head convex, foot concave and tangent to the substrate line).
struct RoundedLine {
  double xc = 0.15, cd = 0.1, h = kH, r = 0.01, a = 15 * pi / 180;
  double s() const { return std::sin(a); }
  double c() const { return std::cos(a); }
  double xf() const { return xc - 0.5 * cd - 0.5 * h * std::tan(a); }        // left foot of the trapezoid
  double cfx() const { return xf() - r * (1 - s()) / c(); }                   // foot arc centre (cfx, r)
  double ctx() const { return xf() + (s() * (h - r) + r) / c(); }             // head arc centre (ctx, h - r)
  double yWF() const { return r * (1 - s()); }                                // wall / foot arc
  double yWT() const { return h - r * (1 - s()); }                            // wall / head arc
  double left(double y) const {                                               // left boundary x(y)
    if (y < yWF()) return cfx() + std::sqrt(r * r - (r - y) * (r - y));
    if (y < yWT()) return xf() + y * std::tan(a);
    return ctx() - std::sqrt(r * r - (y - (h - r)) * (y - (h - r)));
  }
  double g(double x) const {  // profile y = g(x) (profile1d Function)
    const double u = xc - std::abs(x - xc);  // mirrored onto the left half
    if (u <= cfx()) return 0.0;
    if (u <= cfx() + r * c()) return r - std::sqrt(std::max(0.0, r * r - (u - cfx()) * (u - cfx())));
    if (u <= ctx() - r * c()) return (u - xf()) / std::tan(a);
    if (u <= ctx()) return h - r + std::sqrt(std::max(0.0, r * r - (u - ctx()) * (u - ctx())));
    return h;
  }
  std::vector<nf::NurbsCurve> curves() const {
    auto arc = [&](double cx, double cz, double a0, double a1) {
      return nf::make_arc({cx, 0, cz}, {1, 0, 0}, {0, 0, 1}, r, a0, a1);
    };
    auto start = [](const nf::NurbsCurve& k) { return k.evaluate(k.domain_start()); };
    auto end = [](const nf::NurbsCurve& k) { return k.evaluate(k.domain_end()); };
    const double m = 2 * xc, d2r = pi / 180, w = a / d2r;
    const auto footR = arc(m - cfx(), r, (180 + w) * d2r, 270 * d2r).reversed();
    const auto headR = arc(m - ctx(), h - r, w * d2r, 90 * d2r);
    const auto headL = arc(ctx(), h - r, 90 * d2r, (180 - w) * d2r);
    const auto footL = arc(cfx(), r, -90 * d2r, -w * d2r).reversed();
    return {nf::make_line(end(footL), start(footR)), footR, nf::make_line(end(footR), start(headR)), headR,
            nf::make_line(end(headR), start(headL)), headL, nf::make_line(end(headL), start(footL)), footL};
  }
};

json base1D(const json& layers, const std::string& solver = "rcwa", int N = 6) {
  json c = base(layers, solver, N);
  c["dimension"] = "1D";
  c["period_um"] = {{"x", kD}};
  c["truncation"] = {{"Nx", N}};
  c["incidence"] = {{"lambda_um", 0.633}, {"phi_deg", 20}, {"psi_deg", 30}, {"delta_deg", 35}};
  c["rk"] = {{"tol", 1e-10}};
  return c;
}

json cad1D(const std::string& stepText, const json& bodies, const std::string& fac = "laurent") {
  return json{{"type", "cad"},           {"step", stepText},   {"cell_origin_um", 0.0},
              {"background", "superstrate"}, {"bodies", bodies}, {"factorization", fac}};
}

const json kLineBody = json::array({{{"name", "Line"}, {"material", "pillar"}}});

const std::vector<std::array<double, 2>> kTrap{
    {0.0866, 0.0}, {0.2134, 0.0}, {0.2134 - kH * std::tan(15 * pi / 180), kH}, {0.0866 + kH * std::tan(15 * pi / 180), kH}};

fg::Stack stackOf(const json& cfg) { return fg::io::parseConfig(cfg.dump()).stack; }

/// Layers of the band decomposition (bottom to top) with their bottom heights.
/// y0: bottom of the y range (the band boundaries are critical heights, e.g. -1.5e-12 for a foot arc tangent to
/// y = 0, found by Newton in Nurbforge).
std::vector<std::pair<double, const fg::Sectioned1DLayer*>> sectioned1D(const fg::Stack& st, double y0) {
  std::vector<std::pair<double, const fg::Sectioned1DLayer*>> out;
  double y = y0;
  for (const auto& le : st.layers) {
    const double t = std::visit(
        [](const auto& l) -> double {
          using T = std::decay_t<decltype(l)>;
          if constexpr (std::is_same_v<T, fg::Sectioned1DLayer>) return l.provider->thickness_um();
          else if constexpr (requires { l.thickness_um; }) return l.thickness_um;
          else return 0.0;
        },
        le.spec);
    if (const auto* s = std::get_if<fg::Sectioned1DLayer>(&le.spec)) out.push_back({y, s});
    y += t;
  }
  return out;
}

}  // namespace

TEST(CAD_V24, RectangleIsLamellar) {
  const std::string s = step({line({{0.1, 0.0}, {0.2, 0.0}, {0.2, kH}, {0.1, kH}}, 0.3)}, {"Line"});
  for (const std::string solver : {"rcwa", "rk"}) {
    const json cad = base1D(json::array({cad1D(s, kLineBody, "auto")}), solver);
    const json r = report(cad);
    EXPECT_EQ(codes(r), std::vector<std::string>{});
    EXPECT_EQ(r["dimension"], "1D");
    ASSERT_EQ(r["bands"].size(), 1u);
    EXPECT_EQ(r["bands"][0]["kind"], "lamellar1d");
    const json lam = base1D(json::array({{{"type", "lamellar1d"}, {"thickness_um", kH},
                                          {"breakpoints_um", {0.0, 0.1, 0.2, kD}},
                                          {"materials", {"air", "pillar", "air"}}}}),
                            solver);
    EXPECT_LT(maxDiff(solve(cad), solve(lam)), 1e-12) << solver;
  }
}

TEST(CAD_V24, TrapezoidIsProfile1DWithLaurent) {
  const std::string s = step({line(kTrap, 0.3)}, {"Line"});
  const json prof = {{"type", "profile1d"}, {"shape", "points"}, {"depth_um", kH},
                     {"x_um", {0.0, kTrap[0][0], kTrap[3][0], kTrap[2][0], kTrap[1][0]}},
                     {"y_um", {0.0, 0.0, kH, kH, 0.0}}, {"below", "pillar"}, {"above", "air"},
                     {"factorization", "laurent"}};
  const json r = report(base1D(json::array({cad1D(s, kLineBody)})));
  ASSERT_EQ(r["bands"].size(), 1u);
  EXPECT_EQ(r["bands"][0]["kind"], "sectioned1d");
  const double mid = kTrap[1][0] - kTrap[0][0] - kH * std::tan(15 * pi / 180);  // width at half height
  EXPECT_NEAR(r["bands"][0]["area_fractions"]["pillar"].get<double>(), mid / kD, 1e-12);
  EXPECT_LT(maxDiff(solve(base1D(json::array({cad1D(s, kLineBody)}))), solve(base1D(json::array({prof})))), 1e-12);
  EXPECT_LT(maxDiff(solve(base1D(json::array({cad1D(s, kLineBody)}), "rk")), solve(base1D(json::array({prof}), "rk"))),
            1e-10);
}

TEST(CAD_V24, RoundedLineOfTheReport) {
  const RoundedLine L;
  const std::string s = step({profileBody({L.curves()}, 0.3)}, {"Line"});
  const json cfg = base1D(json::array({cad1D(s, kLineBody)}), "rk");
  const json r = report(cfg);
  EXPECT_EQ(codes(r), std::vector<std::string>{});
  ASSERT_EQ(r["bands"].size(), 3u);  // foot arcs, walls, head arcs
  EXPECT_NEAR(r["bands"][0]["y1_um"].get<double>(), L.yWF(), 1e-12);
  EXPECT_NEAR(r["bands"][1]["y1_um"].get<double>(), L.yWT(), 1e-12);
  for (const auto& b : r["bands"]) EXPECT_EQ(b["kind"], "sectioned1d");
  // the line sections: exact boundaries and radial normals on the arcs
  const fg::Stack st = stackOf(cfg);
  for (const auto& [y0, l] : sectioned1D(st, r["bands"][0]["y0_um"].get<double>())) {
    const double t = l->provider->thickness_um();
    for (double f : {0.01, 0.3, 0.5, 0.7, 0.99}) {
      const double y = y0 + f * t;
      const fg::LineSection sec = l->provider->section(f * t);
      ASSERT_EQ(sec.materials, (std::vector<std::string>{"air", "pillar", "air"})) << y;
      EXPECT_NEAR(sec.breaks[1], L.left(y), 1e-13) << y;
      EXPECT_NEAR(sec.breaks[2], 2 * L.xc - L.left(y), 1e-13) << y;
      if (y < L.yWF() || y > L.yWT()) {
        const double cz = y < L.yWF() ? L.r : L.h - L.r, cx = y < L.yWF() ? L.cfx() : L.ctx();
        const fg::Vec3D rad((sec.breaks[1] - cx) / L.r, (y - cz) / L.r, 0);
        EXPECT_LT((sec.normals[1] * sec.normals[1].transpose() - rad * rad.transpose()).cwiseAbs().maxCoeff(), 1e-12);
      }
    }
  }
  // against profile1d with the analytic g(x) (Laurent, RK)
  fg::Stack ref = st;
  ref.layers.clear();
  fg::Profile1DLayer p;
  p.depth_um = L.h;
  p.shape = fg::Profile1DLayer::Shape::Function;
  p.g = [L](double x) { return L.g(x); };
  p.below = "pillar";
  p.above = "air";
  ref.layers.push_back({p, fg::Factorization::Laurent});
  const auto c = fg::io::parseConfig(cfg.dump());
  const auto a = fg::makeSolver(c.solver)->solve(st, c.scan.front());
  const auto b = fg::makeSolver(c.solver)->solve(ref, c.scan.front());
  std::printf("V24c rounded line: |cad - profile1d| = %.3e\n", maxDiff(a, b));
  EXPECT_LT(maxDiff(a, b), 1e-8);
}

TEST(CAD_V24, OverhangCapAndCavity) {
  const auto lineLoops = std::vector<std::vector<nf::NurbsCurve>>{
      polyCurves({{0.11, 0.0}, {0.19, 0.0}, {0.21, kH}, {0.09, kH}}),
      polyCurves({{0.14, 0.03}, {0.16, 0.03}, {0.16, 0.06}, {0.14, 0.06}})};
  const std::string s = step({profileBody(lineLoops, 0.3), line({{0.09, kH}, {0.21, kH}, {0.21, 0.13}, {0.09, 0.13}}, 0.3)},
                             {"Line", "Cap"});
  const json bodies = json::array({{{"name", "Line"}, {"material", "pillar"}}, {{"name", "Cap"}, {"material", "lid"}}});
  json rc = base1D(json::array({cad1D(s, bodies, "auto")}), "rcwa", 8), rk = rc;
  rc["rcwa"] = {{"slices", 100}};
  rk["solver"] = "rk";
  rk["rk"] = {{"tol", 1e-8}};
  const json r = report(rc);
  EXPECT_EQ(codes(r), std::vector<std::string>{});
  std::vector<std::string> kinds;
  for (const auto& b : r["bands"]) kinds.push_back(b["kind"]);
  EXPECT_EQ(kinds, (std::vector<std::string>{"sectioned1d", "sectioned1d", "sectioned1d", "lamellar1d"}));
  const auto a = solve(rc), b = solve(rk);
  auto sum = [](const fg::Result& x) {
    double v = 0;
    for (const auto& o : x.reflected) v += o.efficiency;
    for (const auto& o : x.transmitted) v += o.efficiency;
    return v;
  };
  std::printf("V24d: |RCWA(100) - RK| = %.3e, energy RCWA %.3e RK %.3e\n", maxDiff(a, b), sum(a) - 1, sum(b) - 1);
  EXPECT_NEAR(sum(a), 1.0, 1e-10);
  EXPECT_NEAR(sum(b), 1.0, 1e-7);
  EXPECT_LT(maxDiff(a, b), 2e-3);
  // reciprocity e_n(alpha0) = e_n(-alpha_n) (as V3): exact up to truncation, the deviation falls with N
  const double lam = 0.4, th = 30.0, thRev = std::asin(-(std::sin(th * pi / 180) - lam / kD)) * 180 / pi;
  auto effR = [&](const fg::Result& x, int n) {
    for (const auto& o : x.reflected)
      if (o.n == n) return o.efficiency;
    return -1.0;
  };
  for (double delta : {0.0, 90.0}) {
    double dev[2];
    const int Ns[2] = {8, 24};
    for (int k = 0; k < 2; ++k) {
      json c = rc;
      c["truncation"] = {{"Nx", Ns[k]}};
      c["incidence"] = {{"lambda_um", lam}, {"phi_deg", th}, {"psi_deg", 0}, {"delta_deg", delta}};
      json cr = c;
      cr["incidence"]["phi_deg"] = thRev;
      dev[k] = std::abs(effR(solve(c), -1) - effR(solve(cr), -1));
    }
    std::printf("V24d reciprocity delta %.0f: N = 8 %.3e, N = 24 %.3e\n", delta, dev[0], dev[1]);
    EXPECT_LT(dev[1], 0.25 * dev[0]) << delta;
    EXPECT_LT(dev[1], 1e-3) << delta;
  }
}

TEST(CAD_V24, IndependentOfCutPlaneLengthAndOrigin) {
  const fg::Result ref = solve(base1D(json::array({cad1D(step({line(kTrap, 0.3)}, {"Line"}), kLineBody)})));
  for (double zc : {-0.01, -0.15, -0.29}) {
    json l = cad1D(step({line(kTrap, 0.3)}, {"Line"}), kLineBody);
    l["z_cut_um"] = zc;
    const json cfg = base1D(json::array({l}));
    EXPECT_NEAR(report(cfg)["z_cut_um"].get<double>(), zc, 1e-15);
    EXPECT_LT(maxDiff(solve(cfg), ref), 1e-13) << zc;
  }
  for (double len : {0.1, 1000.0}) {
    const json cfg = base1D(json::array({cad1D(step({line(kTrap, len)}, {"Line"}), kLineBody)}));
    EXPECT_EQ(codes(report(cfg)), std::vector<std::string>{}) << len;  // no UNIT_SUSPECT for long lines
    EXPECT_LT(maxDiff(solve(cfg), ref), 1e-13) << len;
  }
  // shifted over the cell boundary: only phases change
  std::vector<std::array<double, 2>> shifted = kTrap;
  for (auto& p : shifted) p[0] += 0.2;
  EXPECT_LT(maxDiff(solve(base1D(json::array({cad1D(step({line(shifted, 0.3)}, {"Line"}), kLineBody)}))), ref), 1e-12);
}

TEST(CAD_V24, Periodicity) {
  // residual layer and line as one body (spans the period: seam, no interface) = two touching bodies
  const double t = 0.02, w = kH * std::tan(15 * pi / 180);
  const std::string one = step({line({{0.0, 0.0}, {kD, 0.0}, {kD, t}, {0.2134, t}, {0.2134 - w, t + kH},
                                      {0.0866 + w, t + kH}, {0.0866, t}, {0.0, t}},
                                     0.3)},
                               {"Line"});
  std::vector<std::array<double, 2>> up = kTrap;
  for (auto& p : up) p[1] += t;
  const std::string two = step({line({{0.0, 0.0}, {kD, 0.0}, {kD, t}, {0.0, t}}, 0.3), line(up, 0.3)}, {"Slab", "Line"});
  const json both = json::array({{{"name", "Slab"}, {"material", "pillar"}}, {{"name", "Line"}, {"material", "pillar"}}});
  const json r1 = report(base1D(json::array({cad1D(one, kLineBody)})));
  EXPECT_EQ(codes(r1), std::vector<std::string>{});
  ASSERT_EQ(r1["bands"].size(), 2u);
  EXPECT_EQ(r1["bands"][0]["kind"], "homogeneous");
  EXPECT_EQ(r1["bands"][0]["material"], "pillar");
  EXPECT_LT(maxDiff(solve(base1D(json::array({cad1D(one, kLineBody)}))), solve(base1D(json::array({cad1D(two, both)})))),
            1e-12);
  // wider than the period: PERIODIC_OVERLAP; periodic_union gives the homogeneous layer
  json wide = cad1D(step({line({{0.0, 0.0}, {0.35, 0.0}, {0.35, 0.05}, {0.0, 0.05}}, 0.3)}, {"Line"}), kLineBody);
  EXPECT_EQ(codes(report(base1D(json::array({wide})))), std::vector<std::string>{"PERIODIC_OVERLAP"});
  wide["periodic_union"] = true;
  const json hom = {{"type", "homogeneous"}, {"thickness_um", 0.05}, {"material", "pillar"}};
  EXPECT_LT(maxDiff(solve(base1D(json::array({wide}))), solve(base1D(json::array({hom})))), 1e-13);
}

TEST(CAD_V24, Findings) {
  using V = std::vector<std::string>;
  // not z-invariant: a pillar
  EXPECT_EQ(codes(report(base1D(json::array({cad1D(step({cylinder(0.15, 0.25, 0.05, 0, 0.1)}, {"Line"}), kLineBody)})))),
            V{"NOT_Z_INVARIANT"});
  // cutting plane outside the line
  json l = cad1D(step({line(kTrap, 0.3)}, {"Line"}), kLineBody);
  l["z_cut_um"] = 0.1;
  EXPECT_EQ(codes(report(base1D(json::array({l})))), V{"OUT_OF_RANGE"});
  // overlapping bodies, resolved by priority
  const std::string two = step({line({{0.05, 0.0}, {0.15, 0.0}, {0.15, kH}, {0.05, kH}}, 0.3),
                                line({{0.10, 0.0}, {0.20, 0.0}, {0.20, kH}, {0.10, kH}}, 0.3)},
                               {"P", "Q"});
  const json pq = json::array({{{"name", "P"}, {"material", "pillar"}}, {{"name", "Q"}, {"material", "lid"}}});
  EXPECT_EQ(codes(report(base1D(json::array({cad1D(two, pq)})))), V{"BODY_OVERLAP"});
  json pr = cad1D(two, pq);
  pr["overlap"] = "priority";
  const json rp = report(base1D(json::array({pr})));
  EXPECT_EQ(codes(rp), V{});
  EXPECT_EQ(codes(rp, "warnings"), V{"BODY_OVERLAP"});
  const json lam = {{"type", "lamellar1d"}, {"thickness_um", kH}, {"breakpoints_um", {0.0, 0.05, 0.15, 0.2, kD}},
                    {"materials", {"air", "pillar", "lid", "air"}}, {"factorization", "laurent"}};
  EXPECT_LT(maxDiff(solve(base1D(json::array({pr}))), solve(base1D(json::array({lam})))), 1e-12);
  // keys of the other dimension
  json z = cad1D(step({line(kTrap, 0.3)}, {"Line"}), kLineBody);
  z["cell_origin_um"] = {0.0, 0.0};
  EXPECT_THROW((void)report(base1D(json::array({z}))), fg::Error);
}

TEST(CAD_V24, Gradients) {
  // material (adjointSlice) and shape (wall angle via StackParameter) against central differences
  auto trap = [](double top) {  // top width
    const double c = 0.15;
    return std::vector<std::array<double, 2>>{{0.0866, 0.0}, {0.2134, 0.0}, {c + 0.5 * top, kH}, {c - 0.5 * top, kH}};
  };
  fg::Objective obj;
  obj.weights = {{fg::OrderRef{'R', 0, 0}, 1.0}, {fg::OrderRef{'T', -1, 0}, 0.7}};
  auto L = [&](const json& cfg) {
    const auto r = solve(cfg);
    double v = 0;
    for (const auto& e : r.reflected) v += e.n == 0 ? e.efficiency : 0.0;
    for (const auto& e : r.transmitted) v += e.n == -1 ? 0.7 * e.efficiency : 0.0;
    return v;
  };
  for (const std::string solver : {"rcwa", "rk"}) {
    auto cfgFor = [&](double top, double n) {
      json c = base1D(json::array({cad1D(step({line(trap(top), 0.3)}, {"Line"}), kLineBody, "auto")}), solver, 5);
      c["materials"]["pillar"] = {{"n", n}};
      return c;
    };
    const double top = 0.0732, n = 2.0;
    const auto c = fg::io::parseConfig(cfgFor(top, n).dump());
    const auto g = fg::gradient(c.solver, c.stack, c.scan.front(), obj,
                                {fg::StackParameter{"top", [&](double d) { return stackOf(cfgFor(top + d, n)); }}});
    const double h = 1e-5, fd = (L(cfgFor(top + h, n)) - L(cfgFor(top - h, n))) / (2 * h);
    std::printf("V24h %s top width: adjoint %.10e fd %.10e\n", solver.c_str(), g.grad[0], fd);
    EXPECT_NEAR(g.grad[0], fd, 1e-5 * std::max(1.0, std::abs(fd)));
    const auto gm = fg::gradient(cfgFor(top, n).dump(), obj, {"/materials/pillar/n"});
    const double hn = 1e-6, fdn = (L(cfgFor(top, n + hn)) - L(cfgFor(top, n - hn))) / (2 * hn);
    std::printf("V24h %s material: adjoint %.10e fd %.10e\n", solver.c_str(), gm.grad[0], fdn);
    EXPECT_NEAR(gm.grad[0], fdn, 1e-6 * std::max(1.0, std::abs(fdn)));
  }
}

// ------------------------------------------------------------------ periodic_union in 2D, priority rounds (v2.24)
TEST(CAD, PeriodicUnionCannotResolveSelfOverlapIn2D) {
  // B overlaps B + t, so B u (B + t) overlaps its image again: the union wraps around the cell and is no closed
  // contour. Before v2.24 the check passed after four unions and only the solver rejected the layer.
  const json pillar = json::array({{{"name", "P"}, {"material", "pillar"}}});
  json l = cadLayer(step({cylinder(0.25, 0.25, 0.3, 0, 0.2)}, {"P"}), pillar);
  l["periodic_union"] = true;
  const json cfg = base(json::array({l}));
  EXPECT_EQ(codes(report(cfg)), std::vector<std::string>{"PERIODIC_OVERLAP"});
  EXPECT_THROW((void)solve(cfg), fg::Error);
}

TEST(CAD, PriorityResolvesEveryOverlap) {
  // five pairs of overlapping boxes: five Boolean differences, more than the four rounds of v2.21. The second box
  // of a pair lies strictly inside the first in y and z (a planar cut without coplanar faces; coplanar faces:
  // PriorityWithCoplanarLids)
  std::vector<nf::Shell> boxes;
  std::vector<std::string> names;
  json bodies = json::array();
  for (int k = 0; k < 5; ++k) {
    const double x = 0.02 + 0.16 * (k % 3), y = 0.05 + 0.22 * (k / 3);
    boxes.push_back(prism({{x, y}, {x + 0.08, y}, {x + 0.08, y + 0.15}, {x, y + 0.15}}, 0.0, 0.1));
    boxes.push_back(prism({{x + 0.05, y + 0.01}, {x + 0.13, y + 0.01}, {x + 0.13, y + 0.14}, {x + 0.05, y + 0.14}},
                          0.005, 0.095));
    for (int q = 0; q < 2; ++q) {
      names.push_back("B" + std::to_string(2 * k + q));
      bodies.push_back({{"name", names.back()}, {"material", q ? "lid" : "pillar"}});
    }
  }
  json l = cadLayer(step(boxes, names), bodies);
  l["overlap"] = "priority";
  const json r = report(base(json::array({l})));
  EXPECT_EQ(codes(r), std::vector<std::string>{});
  EXPECT_EQ(codes(r, "warnings").size(), 5u);
  EXPECT_NO_THROW((void)solve(base(json::array({l}))));
}

TEST(CAD, PriorityWithCoplanarLids) {
  // two boxes of the same y range overlap (coplanar lids and coplanar side faces): Nurbforge v0.6.0 returned an
  // open shell for this difference (BODY_OVERLAP since v2.24), v0.7.0 fixed it (Nurbforge #45). The earlier body
  // wins: the result equals the model in which the overlap is already cut out of Q
  const auto P = prism({{0.05, 0.2}, {0.15, 0.2}, {0.15, 0.3}, {0.05, 0.3}}, 0.0, 0.1),
             Q = prism({{0.12, 0.2}, {0.22, 0.2}, {0.22, 0.3}, {0.12, 0.3}}, 0.0, 0.1),
             Qcut = prism({{0.15, 0.2}, {0.22, 0.2}, {0.22, 0.3}, {0.15, 0.3}}, 0.0, 0.1);
  const json pq = json::array({{{"name", "P"}, {"material", "pillar"}}, {{"name", "Q"}, {"material", "lid"}}});
  json l = cadLayer(step({P, Q}, {"P", "Q"}), pq);
  l["overlap"] = "priority";
  const json r = report(base(json::array({l})));
  EXPECT_EQ(codes(r), std::vector<std::string>{});
  EXPECT_EQ(codes(r, "warnings"), std::vector<std::string>{"BODY_OVERLAP"});
  const json ref = cadLayer(step({P, Qcut}, {"P", "Q"}), pq);
  EXPECT_LT(maxDiff(solve(base(json::array({l}))), solve(base(json::array({ref})))), 1e-10);
}

TEST(CAD, OverlapWithCollinearEdgesIsFound) {
  // two boxes of the same z range: the boundaries only run along each other (no proper crossing); before v2.24
  // neither the slicer nor the solver noticed the overlap
  const std::string s = step({prism({{0.05, 0.2}, {0.15, 0.2}, {0.15, 0.3}, {0.05, 0.3}}, 0.0, 0.1),
                              prism({{0.12, 0.2}, {0.22, 0.2}, {0.22, 0.3}, {0.12, 0.3}}, 0.0, 0.1)},
                             {"P", "Q"});
  const json pq = json::array({{{"name", "P"}, {"material", "pillar"}}, {{"name", "Q"}, {"material", "lid"}}});
  EXPECT_EQ(codes(report(base(json::array({cadLayer(s, pq)})))), std::vector<std::string>{"BODY_OVERLAP"});
  // with priority the coplanar difference works since Nurbforge v0.7.0: CAD.PriorityWithCoplanarLids
  // touching boxes (common edge) are fine
  const std::string t = step({prism({{0.05, 0.2}, {0.15, 0.2}, {0.15, 0.3}, {0.05, 0.3}}, 0.0, 0.1),
                              prism({{0.15, 0.2}, {0.25, 0.2}, {0.25, 0.3}, {0.15, 0.3}}, 0.0, 0.1)},
                             {"P", "Q"});
  EXPECT_EQ(codes(report(base(json::array({cadLayer(t, pq)})))), std::vector<std::string>{});
  // the solver's check of contour2d regions (same configuration as polygons)
  auto rect = [](double x0, double x1, const std::string& m) {
    return json{{"material", m}, {"polygon_um", {{x0, 0.2}, {x1, 0.2}, {x1, 0.3}, {x0, 0.3}}}};
  };
  json c = base(json::array({{{"type", "contour2d"}, {"thickness_um", 0.1}, {"background", "air"},
                               {"regions", {rect(0.05, 0.15, "pillar"), rect(0.12, 0.22, "lid")}}}}));
  EXPECT_THROW((void)solve(c), fg::Error);
  c["layers"][0]["regions"] = {rect(0.05, 0.15, "pillar"), rect(0.15, 0.25, "lid")};
  EXPECT_NO_THROW((void)solve(c));
}

// ------------------------------------------------------------------ M20: solids of revolution
// cad with "axisymmetric": true in a cylindrical stack (Plan "M20: Schichttypen"): cylinder, frustum and
// spherical cap on one axis -> bands rings / rings_profile; the result equals the hand-written stack with
// rings and rings_profile (cone, sphere); sections that are no circles about the axis are NOT_AXISYMMETRIC.
namespace {

json cylBase(const json& layers) {
  return json{{"cylindrical", {{"spectrum", {{"kappa_max", 2}}}, {"factorization", "li_r"}}},
              {"materials", {{"air", {{"n", 1.0}}}, {"glass", {{"n", 1.5}}}, {"pillar", {{"n", 1.6}}}, {"lid", {{"n", 1.4}}}}},
              {"superstrate", "air"},
              {"substrate", "glass"},
              {"incidence",
               {{"type", "gaussian_beam"}, {"lambda_um", 1.0}, {"waist_um", 0.5}, {"focus_y_um", 0.9}, {"polarization", "x"}}},
              {"layers", layers}};
}

json axisymLayer(const std::string& stepText, const json& bodies, const json& extra = json::object()) {
  json l{{"type", "cad"}, {"step", stepText}, {"axisymmetric", true}, {"background", "superstrate"},
         {"bodies", bodies}, {"slicing", {{"rcwa_slices", 6}}}};
  for (const auto& [k, v] : extra.items()) l[k] = v;
  return l;
}

}  // namespace

TEST(CAD_M20, SolidOfRevolutionEqualsRingsAndRingsProfile) {
  // axis at CAD (X, Y) = (2, 3): FG axis (x, z) = (2, -3) (FG z = -CAD Y), the default origin
  const double cx = 2, cy = 3;
  const std::string s = step({cylinder(cx, cy, 0.4, 0.0, 0.2), frustum(cx, cy, 0.3, 0.15, 0.2, 0.25),
                              cap(cx, cy, 0.3, 0.2, 0.45)},
                             {"base", "cone", "dome"});
  const json bodies = json::array({{{"name", "base"}, {"material", "pillar"}}, {{"name", "cone"}, {"material", "lid"}},
                                   {{"name", "dome"}, {"material", "pillar"}}});
  const json cfg = cylBase(json::array({axisymLayer(s, bodies)}));
  const json rep = report(cfg);
  EXPECT_TRUE(rep["ok"].get<bool>()) << rep.dump(2);
  EXPECT_EQ(rep["dimension"], "axisymmetric");
  EXPECT_NEAR(rep["axis_um"][0].get<double>(), 2.0, 1e-12);
  EXPECT_NEAR(rep["axis_um"][1].get<double>(), -3.0, 1e-12);
  ASSERT_EQ(rep["bands"].size(), 3u);
  EXPECT_EQ(rep["bands"][0]["kind"], "rings");
  EXPECT_EQ(rep["bands"][1]["kind"], "rings_profile");
  EXPECT_EQ(rep["bands"][2]["kind"], "rings_profile");
  const auto c = fg::io::parseConfig(cfg.dump());
  ASSERT_EQ(c.stack.layers.size(), 3u);
  const auto* r0 = std::get_if<fg::RingsLayer>(&c.stack.layers[0].spec);
  ASSERT_NE(r0, nullptr);
  ASSERT_EQ(r0->radii_um.size(), 1u);
  EXPECT_NEAR(r0->radii_um[0], 0.4, 1e-12);
  // the cap's centre lies 0.1 below its band (z0 + hc - rho = 0.35)
  const json hand = cylBase(json::array(
      {{{"type", "rings"}, {"thickness_um", 0.2}, {"radii_um", {0.4}}, {"materials", {"pillar", "air"}}},
       {{"type", "rings_profile"}, {"thickness_um", 0.25}, {"slices", 6}, {"radii", {{{"cone", {0.3, 0.15}}}}},
        {"materials", {"lid", "air"}}},
       {{"type", "rings_profile"}, {"thickness_um", 0.2}, {"slices", 6},
        {"radii", {{{"sphere", {{"center_y_um", -0.1}, {"radius_um", 0.3}}}}}}, {"materials", {"pillar", "air"}}}}));
  const fg::Result a = solve(cfg), b = solve(hand);
  std::printf("CAD M20 solid of revolution: R = %.12f (hand %.1e), T = %.12f (hand %.1e)\n", a.cylindrical.reflected,
              a.cylindrical.reflected - b.cylindrical.reflected, a.cylindrical.transmitted,
              a.cylindrical.transmitted - b.cylindrical.transmitted);
  EXPECT_NEAR(a.cylindrical.reflected, b.cylindrical.reflected, 1e-11);
  EXPECT_NEAR(a.cylindrical.transmitted, b.cylindrical.transmitted, 1e-11);
}

TEST(CAD_M20, FindingsAndConfiguration) {
  const std::string off = step({cylinder(2, 3, 0.4, 0.0, 0.2)}, {"base"});
  const json bodies = json::array({{{"name", "base"}, {"material", "pillar"}}});
  // the cylinder is not about the given axis
  json rep = report(cylBase(json::array({axisymLayer(off, bodies, {{"cell_origin_um", {2.1, -3.0}}})})));
  EXPECT_FALSE(rep["ok"].get<bool>());
  EXPECT_EQ(codes(rep), std::vector<std::string>{"NOT_AXISYMMETRIC"});
  EXPECT_NEAR(rep["errors"][0]["where_um"][1].get<double>(), 0.1, 1e-9);  // band mid height
  // a square prism is no solid of revolution
  rep = report(cylBase(json::array({axisymLayer(step({prism(kL, 0.0, 0.2)}, {"base"}), bodies)})));
  EXPECT_EQ(codes(rep), std::vector<std::string>{"NOT_AXISYMMETRIC"});
  // an annulus (tube) is a ring
  const json tube = cylBase(json::array({axisymLayer(
      step({nf::boolean(cylinder(0, 0, 0.4, 0.0, 0.2), cylinder(0, 0, 0.2, -0.1, 0.3), nf::BooleanOp::Difference)},
           {"base"}),
      bodies)}));
  const auto t = fg::io::parseConfig(tube.dump());
  const auto* tr = std::get_if<fg::RingsLayer>(&t.stack.layers.at(0).spec);
  ASSERT_NE(tr, nullptr);
  ASSERT_EQ(tr->radii_um.size(), 2u);
  EXPECT_NEAR(tr->radii_um[0], 0.2, 1e-12);
  EXPECT_NEAR(tr->radii_um[1], 0.4, 1e-12);
  EXPECT_EQ(tr->materials, (std::vector<std::string>{"air", "pillar", "air"}));
  // axisymmetric only in a cylindrical stack, and there only axisymmetric
  json noFlag = cylBase(json::array({axisymLayer(off, bodies)}));
  noFlag["layers"][0].erase("axisymmetric");
  EXPECT_THROW(fg::io::parseConfig(noFlag.dump()), fg::Error);
  EXPECT_THROW(fg::io::parseConfig(base(json::array({axisymLayer(off, bodies)})).dump()), fg::Error);
}

// V25c (Plan section "M24: Tests und Abnahme"): the rounded line of the example as profile1d "nurbs" (exact
// rational arcs, one layer with band boundaries) equals the CAD line (STEP, three band layers); RK for both (the
// RCWA staircases differ: 30 slices per band against 30 over the layer). The step sequences differ (one layer
// against three), so the agreement is that of the RK tolerance 1e-10: measured 3.2e-10 (CI)
TEST(CAD_V25, NurbsProfileEqualsTheCadLine) {
  const std::string ex = std::string(FG_SOURCE_DIR) + "/examples/";
  fg::io::RunConfig cad = fg::io::loadConfig(ex + "cad_line_grating_1d.json"),
                    nrb = fg::io::loadConfig(ex + "profile1d_nurbs.json");
  for (auto* c : {&cad, &nrb}) {
    c->solver.kind = fg::SolverKind::RK;
    c->solver.rkTol = 1e-10;
  }
  ASSERT_EQ(nrb.stack.layers.size(), 1u);
  const auto a = fg::makeSolver(cad.solver)->solve(cad.stack, cad.scan.front());
  const auto b = fg::makeSolver(nrb.solver)->solve(nrb.stack, nrb.scan.front());
  std::printf("V25c: |profile1d nurbs - cad| = %.3e (RK)\n", maxDiff(a, b));
  EXPECT_LT(maxDiff(a, b), 1e-9);
}
