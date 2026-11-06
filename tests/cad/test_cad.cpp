// CAD slicer (Plan section "M14: Tests und Abnahme"): STEP files written at test time with
// Nurbforge's native writer (unit micrometre, named solids) -> cad layers -> V16b, V17, V18, two
// materials with a cavity, Chebyshev interpolation, faulty files with their codes.
#include <gtest/gtest.h>

#include <cmath>
#include <nlohmann/json.hpp>
#include <numbers>
#include <nurbforge/geometry/primitives.hpp>
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
                               : 2 * pi * R * std::cyl_bessel_j(1.0, k * R) / (k * 0.25) *
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
  const fg::Objective obj{{{fg::OrderRef{'R', 0, 0}, 1.0}, {fg::OrderRef{'T', -1, 0}, 0.7}}};
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
