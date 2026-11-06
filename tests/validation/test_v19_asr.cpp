// Validation M17 (Plan section "M17: Tests und Abnahme", docs/validation.md): coordinate
// transformations and the adaptive spatial resolution (ASR) of Granet. eta = 0 bit-identical; the
// transformed path with A = diag(alpha) equals the diagonal path; V1-V3, V9 with eta = 0.5; V19
// (gold lamellae, TM) convergence with and without ASR.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {

constexpr double DEG = fg::PI / 180.0;

fg::LayerEntry entry(fg::LayerSpec s, fg::Factorization f = fg::Factorization::Auto) {
  fg::LayerEntry e;
  e.spec = std::move(s);
  e.factorization = f;
  return e;
}

fg::Result run(const fg::Stack& st, const fg::Wavelength& wl, int Nx, int Nz = 0,
               fg::SolverKind k = fg::SolverKind::RCWA, int slices = 16) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = Nx;
  o.Nz = Nz;
  o.rcwaSlices = slices;
  o.rkTol = 1e-10;
  return fg::makeSolver(o)->solve(st, wl);
}

/// Largest efficiency difference over all orders of both results (missing orders count fully).
double diff(const fg::Result& a, const fg::Result& b) {
  double d = std::abs(a.absorption - b.absorption);
  auto side = [&](const std::vector<fg::OrderEfficiency>& x, const std::vector<fg::OrderEfficiency>& y) {
    for (const auto& o : x) {
      double e = 0;
      for (const auto& p : y)
        if (p.n == o.n && p.m == o.m) e = p.efficiency;
      d = std::max(d, std::abs(o.efficiency - e));
    }
  };
  side(a.reflected, b.reflected);
  side(b.reflected, a.reflected);
  side(a.transmitted, b.transmitted);
  side(b.transmitted, a.transmitted);
  return d;
}

/// Bit-identical efficiencies (and absorption).
bool identical(const fg::Result& a, const fg::Result& b) {
  if (a.reflected.size() != b.reflected.size() || a.transmitted.size() != b.transmitted.size()) return false;
  for (std::size_t k = 0; k < a.reflected.size(); ++k)
    if (a.reflected[k].efficiency != b.reflected[k].efficiency) return false;
  for (std::size_t k = 0; k < a.transmitted.size(); ++k)
    if (a.transmitted[k].efficiency != b.transmitted[k].efficiency) return false;
  return a.absorption == b.absorption;
}

double eff(const std::vector<fg::OrderEfficiency>& v, int n) {
  for (const auto& o : v)
    if (o.n == n) return o.efficiency;
  return 0.0;
}

fg::AxisTransformSpec asr(double eta, std::vector<double> br = {}) {
  fg::AxisTransformSpec a;
  a.type = fg::AxisTransformSpec::Type::ASR;
  a.eta = eta;
  a.breakpoints_um = std::move(br);
  return a;
}

fg::Stack base(double dx, double dz = 0.0) {
  fg::Stack st;
  st.dx_um = dx;
  st.dz_um = dz;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["layer"] = fg::Material::fromIndex("layer", cplx(2.1, 0.0));
  st.materials["metal"] = fg::Material::fromIndex("metal", cplx(0.6, 2.5));
  st.materials["al"] = fg::Material::fromIndex("al", cplx(1.38, 7.62));
  st.materials["au"] = fg::Material::fromIndex("au", cplx(0.2, 3.4));
  return st;
}

fg::Stack lamellar(cplx nBar, double dz = 0.0) {
  fg::Stack st = base(1.0, dz);
  st.materials["bar"] = fg::Material::fromIndex("bar", nBar);
  st.layers = {entry(fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"bar", "superstrate"}})};
  return st;
}

/// The lamella of lamellar() as a z-invariant polygon strip in a 2D stack.
fg::Stack strip(cplx nBar, double dz) {
  fg::Stack st = base(1.0, dz);
  st.materials["bar"] = fg::Material::fromIndex("bar", nBar);
  fg::Polygon2DLayer p;
  p.thickness_um = 0.4;
  p.background = "superstrate";
  p.inclusions.push_back({{{0.0, 0.0}, {0.45, 0.0}, {0.45, dz}, {0.0, dz}}, "bar"});
  st.layers = {entry(p)};
  return st;
}

fg::Profile1DLayer sinus1D() {
  fg::Profile1DLayer p;
  p.depth_um = 0.2;
  p.x0_um = 0.1;
  p.below = "al";
  p.above = "superstrate";
  return p;
}

/// V19: gold lamellae (nu = 0.2 + 3.4i, V13) on gold, period 1 um, width 0.5 um, height 0.3 um.
fg::Stack v19() {
  fg::Stack st = base(1.0);
  st.substrate = "au";
  st.layers = {entry(fg::Lamellar1DLayer{0.3, {0.0, 0.5, 1.0}, {"au", "superstrate"}})};
  return st;
}
const fg::Wavelength v19wl{0.6328, 20 * DEG, 0, 90 * DEG};  // TM

const fg::Wavelength conical{0.6328, 20 * DEG, 15 * DEG, 30 * DEG};

}  // namespace

// ---------------------------------------------------------------- eta = 0
// eta = 0 is no transformation: bit-identical results (Plan "eta = 0: bitgleich").
TEST(M17_Transform, EtaZeroBitIdentical) {
  for (cplx nb : {cplx(2.0, 0), cplx(0.5, 3.0)}) {
    auto st = lamellar(nb);
    const auto a = run(st, conical, 10);
    st.transform.x = asr(0.0);
    EXPECT_TRUE(identical(a, run(st, conical, 10)));
  }
  auto s2 = strip(2.0, 0.7);
  const auto a = run(s2, conical, 5, 2);
  s2.transform.x = asr(0.0);
  s2.transform.z = asr(0.0, {0.2});
  EXPECT_TRUE(identical(a, run(s2, conical, 5, 2)));
}

// The transformed path (dense A, G, numerical half-space modes, conversions eq:ct_inc /
// eq:ct_rayleigh) with f' = g' = 1 equals the diagonal path (Plan "Ersetzungsregel"). The breakpoint
// only fixes the quadrature panels; the exactly degenerate TE/TM pairs of the homogeneous media test
// the structured modes (homogeneousModes in solver.cpp).
TEST(M17_Transform, MatrixPathWithIdentityEqualsDiagonalPath) {
  auto check = [](fg::Stack st, int Nx, int Nz, fg::SolverKind k, double tol, const char* what) {
    const auto a = run(st, conical, Nx, Nz, k);
    st.transform.forceMatrixPath = true;
    st.transform.x = asr(0.0, {0.2});
    if (st.is2D()) st.transform.z = asr(0.0, {0.3});
    const auto b = run(st, conical, Nx, Nz, k);
    EXPECT_LT(diff(a, b), tol) << what;
  };
  auto flat = base(0.7);
  flat.layers = {entry(fg::HomogeneousLayer{0.37, "metal"})};
  check(flat, 3, 0, fg::SolverKind::RCWA, 1e-13, "V2 homogeneous");
  check(lamellar(2.0), 10, 0, fg::SolverKind::RCWA, 1e-12, "V3 dielectric");
  check(lamellar(cplx(0.5, 3.0)), 10, 0, fg::SolverKind::RCWA, 1e-12, "V3 metal");
  check(lamellar(2.0), 8, 0, fg::SolverKind::RK, 1e-10, "V3 RK");
  check(strip(2.0, 0.7), 5, 2, fg::SolverKind::RCWA, 1e-12, "V9 polygon strip");
  auto prof = base(0.8);
  prof.layers = {entry(sinus1D())};
  check(prof, 6, 0, fg::SolverKind::RCWA, 1e-12, "sinus profile");
}

// ---------------------------------------------------------------- V1, V2 with eta = 0.5
// Flat stacks: the transformed orders converge to the analytic result. One breakpoint gives an
// analytic f (exponential convergence), several give a jump of f''' (algebraic, ~N^-5).
TEST(M17_Transform, V1V2FlatStacks) {
  auto flat = base(0.7);
  auto layer = base(0.7);
  layer.layers = {entry(fg::HomogeneousLayer{0.37, "metal"})};
  for (fg::Stack* st : {&flat, &layer})
    for (double delta : {0.0, 90.0, 30.0}) {
      const fg::Wavelength wl{0.6328, 35 * DEG, 20 * DEG, delta * DEG};
      const auto ref = run(*st, wl, 3);
      fg::Stack t = *st;
      t.transform.x = asr(0.5, {0.2});
      EXPECT_LT(diff(ref, run(t, wl, 6)), 1e-12) << "one breakpoint, delta " << delta;
      t.transform.x = asr(0.5, {0.2, 0.45});
      EXPECT_LT(diff(ref, run(t, wl, 16)), 3e-8) << "two breakpoints, delta " << delta;
    }
}

// ---------------------------------------------------------------- V3 with eta = 0.5
TEST(M17_Transform, V3Lamellar) {
  for (cplx nb : {cplx(2.0, 0), cplx(0.5, 3.0)}) {
    auto st = lamellar(nb);
    const auto ref = run(st, conical, 120);
    st.transform.x = asr(0.5);
    const auto a = run(st, conical, 40);
    // same limit; the truncation error is smaller with ASR (measured: 1.1e-6 / 1.9e-4 vs 4.2e-6 / 4.8e-4)
    EXPECT_LT(diff(a, ref), nb.imag() == 0 ? 3e-6 : 4e-4) << nb;
    if (nb.imag() == 0) {
      EXPECT_NEAR(a.diag.energyBalance, 1.0, 1e-7);  // ~N^-4 through the conversion of the orders
      // RK = RCWA in transformed coordinates
      EXPECT_LT(diff(run(st, conical, 8, 0, fg::SolverKind::RK), run(st, conical, 8)), 1e-8);
    }
  }
}

// ---------------------------------------------------------------- V9 with eta = 0.5
// z-invariant 2D structures with the same x transformation (and a z transformation of the
// homogeneous z direction) reproduce the transformed 1D conical result.
TEST(M17_Transform, V9ZInvariant) {
  auto s1 = lamellar(2.0);
  s1.transform.x = asr(0.5);
  const auto r1 = run(s1, conical, 8);
  auto s2 = strip(2.0, 0.7);  // mapped polygon contour, resampled normal field (FFF_sym)
  s2.transform.x = asr(0.5);
  EXPECT_LT(diff(r1, run(s2, conical, 8, 0)), 1e-10);
  EXPECT_LT(diff(r1, run(s2, conical, 8, 2)), 1e-10);
  s2.transform.z = asr(0.5, {0.2});
  EXPECT_LT(diff(r1, run(s2, conical, 8, 2)), 1e-9);
  // sinus profile (1D) = profile in the 2D stack = sinus_x
  auto p1 = base(0.8);
  p1.layers = {entry(sinus1D())};
  p1.transform.x = asr(0.5, {0.1, 0.5});
  const auto q1 = run(p1, conical, 6);
  auto p2 = p1;
  p2.dz_um = 0.6;
  fg::Profile2DLayer sx;
  sx.shape = fg::Profile2DLayer::Shape::SinusX;
  sx.depth_um = 0.2;
  sx.x0_um = 0.1;
  sx.below = "al";
  sx.above = "superstrate";
  p2.layers = {entry(sx)};
  EXPECT_LT(diff(q1, run(p2, conical, 6, 0)), 1e-12);
  EXPECT_LT(diff(q1, run(p2, conical, 6, 1)), 1e-10);
  // and the transformed profile converges to the untransformed one (dielectric: the aluminium profile
  // converges too slowly in conical incidence for a meaningful comparison; measured 3.5e-4 at N = 20)
  auto d1 = p1;
  auto sp = sinus1D();
  sp.below = "layer";
  d1.layers = {entry(sp)};
  auto d0 = d1;
  d0.transform.x = {};
  EXPECT_LT(diff(run(d0, conical, 20), run(d1, conical, 20)), 1e-3);
}

// ---------------------------------------------------------------- near field
// Airy layer in transformed coordinates: the near field at physical points equals the untransformed one.
TEST(M17_Transform, NearFieldFlatLayer) {
  auto st = base(0.7);
  st.layers = {entry(fg::HomogeneousLayer{0.37, "layer"})};
  fg::NearFieldRequest req;
  req.nx = 7;
  req.ny = 5;
  req.y_min_um = -0.1;
  req.y_max_um = 0.5;
  req.horizontal_y_um = 0.2;
  fg::SolverOptions o;
  o.Nx = 3;
  const auto a = fg::makeSolver(o)->nearField(st, conical, req);
  st.transform.x = asr(0.5, {0.2});
  o.Nx = 8;
  const auto b = fg::makeSolver(o)->nearField(st, conical, req);
  double d = 0;
  for (int c = 0; c < 3; ++c)
    for (std::size_t k = 0; k < a.vertical.E[c].size(); ++k) {
      d = std::max(d, std::abs(a.vertical.E[c][k] - b.vertical.E[c][k]));
      d = std::max(d, std::abs(a.vertical.H[c][k] - b.vertical.H[c][k]));
    }
  EXPECT_LT(d, 1e-10);
}

// ---------------------------------------------------------------- V19
// Gold lamellae, TM: reference without ASR from N = 600, 800 (split polarizations, Richardson in 1/N):
// e_-1 = 0.013591, e_0 = 0.79435 (N = 800 itself: 0.0135760, 0.7942045). ASR (eta = 0.95) at N = 160:
// 0.0135853, 0.7942996. Convergence table in docs/validation.md.
TEST(M17_Transform, V19GoldLamellaeTM) {
  const double ref1 = 0.013591, ref0 = 0.79435;
  auto st = v19();
  st.transform.x = asr(0.8);
  const auto a = run(st, v19wl, 80);
  EXPECT_NEAR(eff(a.reflected, -1), ref1, 3e-5);
  EXPECT_NEAR(eff(a.reflected, 0), ref0, 3e-4);
  st.transform.x = asr(0.95);
  const auto b = run(st, v19wl, 40);
  EXPECT_NEAR(eff(b.reflected, -1), ref1, 1e-5);
  EXPECT_NEAR(eff(b.reflected, 0), ref0, 1e-4);
  // 1e-3 already at N = 10 with eta = 0.95 (without ASR: N > 50)
  const auto c = run(st, v19wl, 10);
  EXPECT_NEAR(eff(c.reflected, 0), ref0, 1e-3);
}

// ---------------------------------------------------------------- configuration
TEST(M17_Transform, Configuration) {
  const std::string j = R"({"dimension": "1D", "period_um": {"x": 1.0}, "truncation": {"Nx": 10},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "bar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"lambda_um": 0.6328, "phi_deg": 20},
    "layers": [{"type": "lamellar1d", "thickness_um": 0.4, "breakpoints_um": [0.0, 0.45, 1.0], "materials": ["bar", "air"]}],
    "coordinate_transform": {"x": {"type": "asr", "eta": 0.5, "breakpoints_um": "auto"}}})";
  const auto cfg = fg::io::parseConfig(j);
  EXPECT_EQ(cfg.stack.transform.x.type, fg::AxisTransformSpec::Type::ASR);
  EXPECT_DOUBLE_EQ(cfg.stack.transform.x.eta, 0.5);
  const auto tr = fg::resolveTransform(cfg.stack);
  ASSERT_TRUE(tr.x);
  ASSERT_EQ(tr.x->breakpoints().size(), 2u);  // material jumps at 0 and 0.45
  EXPECT_DOUBLE_EQ(tr.x->breakpoints()[1], 0.45);
  EXPECT_FALSE(tr.z);
  auto bad = j;
  bad.replace(bad.find("0.5, \"breakpoints_um\""), 3, "1.0");
  EXPECT_THROW(fg::io::parseConfig(bad), fg::Error);
}
