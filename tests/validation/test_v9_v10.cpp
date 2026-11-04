// Validation V9 (2D code with z-invariant structure = 1D conical) and V10 (Ag hole array,
// book IX.1.7); Plan section "Validierung und Tests", docs/validation.md.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <sstream>

#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {

constexpr double DEG = fg::PI / 180.0;

const fg::Wavelength conical{0.6328, 30 * DEG, 40 * DEG, 25 * DEG};

fg::Stack baseStack(double dz) {
  fg::Stack st;
  st.dx_um = 0.8;
  st.dz_um = dz;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["al"] = fg::Material::fromIndex("al", cplx(1.38, 7.62));
  st.superstrate = "air";
  st.substrate = "glass";
  return st;
}

fg::LayerEntry entry(fg::LayerSpec s, fg::Factorization f = fg::Factorization::Auto) {
  fg::LayerEntry e;
  e.spec = std::move(s);
  e.factorization = f;
  return e;
}

fg::Lamellar1DLayer lamella(const std::string& m) {
  return {0.35, {0.0, 0.3, 0.8}, {m, "air"}};
}

/// The same lamella as a z-invariant polygon (strip over the full z period).
fg::Polygon2DLayer strip(const std::string& m, double dz) {
  fg::Polygon2DLayer p;
  p.thickness_um = 0.35;
  p.background = "air";
  p.inclusions.push_back({{{0.0, 0.0}, {0.3, 0.0}, {0.3, dz}, {0.0, dz}}, m});
  return p;
}

fg::Profile1DLayer sinus1D() {
  fg::Profile1DLayer p;
  p.depth_um = 0.2;
  p.x0_um = 0.1;
  p.below = "al";
  p.above = "air";
  return p;
}

fg::Profile2DLayer sinusX() {
  fg::Profile2DLayer p;
  p.shape = fg::Profile2DLayer::Shape::SinusX;
  p.depth_um = 0.2;
  p.x0_um = 0.1;
  p.below = "al";
  p.above = "air";
  return p;
}

fg::Result run(const fg::Stack& st, fg::SolverKind k, int Nx, int Nz, const fg::Wavelength& wl = conical,
               int slices = 16, double tol = 1e-9) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = Nx;
  o.Nz = Nz;
  o.rcwaSlices = slices;
  o.rkTol = tol;
  return fg::makeSolver(o)->solve(st, wl);
}

/// Every 1D order n must appear as (n, 0) in the 2D result with the same efficiency; all other
/// 2D orders carry no power. Returns the largest deviation.
double compare(const fg::Result& r1, const fg::Result& r2) {
  double dev = 0;
  auto side = [&](const std::vector<fg::OrderEfficiency>& a, const std::vector<fg::OrderEfficiency>& b) {
    for (const auto& o1 : a) {
      bool found = false;
      for (const auto& o2 : b)
        if (o2.n == o1.n && o2.m == 0) {
          dev = std::max(dev, std::abs(o1.efficiency - o2.efficiency));
          found = true;
        }
      if (!found) dev = std::max(dev, 1.0);
    }
    for (const auto& o2 : b)
      if (o2.m != 0) dev = std::max(dev, std::abs(o2.efficiency));
  };
  side(r1.reflected, r2.reflected);
  side(r1.transmitted, r2.transmitted);
  return std::max(dev, std::abs(r1.absorption - r2.absorption));
}

}  // namespace

// ---------------------------------------------------------------- V9
// Lamellar1D layers inside 2D stacks (Nz = 0 and Nz = 2) reproduce the 1D conical result.
TEST(V9_ZInvariant, LamellarIn2DStack) {
  for (const char* m : {"glass", "al"}) {
    auto s1 = baseStack(0.0), s2 = baseStack(0.7);
    for (auto* s : {&s1, &s2}) s->layers = {entry(lamella(m))};
    for (auto k : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
      const auto r1 = run(s1, k, 8, 0);
      EXPECT_LT(compare(r1, run(s2, k, 8, 0)), 1e-12) << m;
      if (k == fg::SolverKind::RCWA) EXPECT_LT(compare(r1, run(s2, k, 8, 2)), 1e-10) << m;
    }
  }
}

// Polygon2D strip (exact polygon coefficients, Algorithm NF, FFF_sym / FFF_iso) = Li lamellar 1D.
TEST(V9_ZInvariant, PolygonStripEqualsLamellar) {
  for (const char* m : {"glass", "al"}) {
    auto s1 = baseStack(0.0);
    s1.layers = {entry(lamella(m))};
    for (auto k : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
      const auto r1 = run(s1, k, 8, 0);
      for (auto f : {fg::Factorization::Auto, fg::Factorization::FFF_iso}) {
        auto s2 = baseStack(0.7);
        s2.layers = {entry(strip(m, 0.7), f)};
        EXPECT_LT(compare(r1, run(s2, k, 8, 0)), 1e-10) << m;
        if (k == fg::SolverKind::RCWA) EXPECT_LT(compare(r1, run(s2, k, 8, 2)), 1e-10) << m;
      }
    }
  }
}

// Profile2D sinus_x (and a 1D profile in a 2D stack) = 1D conical sinus profile, RK and RCWA.
TEST(V9_ZInvariant, ProfileSinusX) {
  auto s1 = baseStack(0.0);
  s1.layers = {entry(sinus1D())};
  auto s2 = baseStack(0.6), s3 = baseStack(0.6);
  s2.layers = {entry(sinusX())};
  s3.layers = {entry(sinus1D())};
  for (auto k : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
    const auto r1 = run(s1, k, 5, 0, conical, 16, 1e-6);
    EXPECT_LT(compare(r1, run(s2, k, 5, 0, conical, 16, 1e-6)), 1e-12);
    EXPECT_LT(compare(r1, run(s3, k, 5, 0, conical, 16, 1e-6)), 1e-12);
    if (k == fg::SolverKind::RCWA) {
      EXPECT_LT(compare(r1, run(s2, k, 5, 1, conical, 16)), 1e-10);
    } else {
      // Nz > 0: the step control sees the decoupled m != 0 blocks -> agreement on the RK tolerance
      EXPECT_LT(compare(r1, run(s2, k, 5, 1, conical, 16, 1e-6)), 1e-4);
    }
  }
}

// Mixed stack: homogeneous layer, polygon strip and z-invariant profile in one 2D stack.
TEST(V9_ZInvariant, MixedStack) {
  auto s1 = baseStack(0.0), s2 = baseStack(0.5);
  s1.layers = {entry(sinus1D()), entry(fg::HomogeneousLayer{0.1, "glass"}), entry(lamella("glass"))};
  s2.layers = {entry(sinusX()), entry(fg::HomogeneousLayer{0.1, "glass"}), entry(strip("glass", 0.5))};
  EXPECT_LT(compare(run(s1, fg::SolverKind::RCWA, 6, 0), run(s2, fg::SolverKind::RCWA, 6, 0)), 1e-10);
  EXPECT_LT(compare(run(s1, fg::SolverKind::RCWA, 6, 0), run(s2, fg::SolverKind::RCWA, 6, 2)), 1e-10);
}

// ---------------------------------------------------------------- 2D sanity
// Lossless crossed grating (dielectric egg crate): energy within the FFF tolerance, RK = RCWA.
TEST(Profile2D, EggCrateEnergyAndRKvsRCWA) {
  fg::Stack st = baseStack(0.7);
  fg::Profile2DLayer p;
  p.shape = fg::Profile2DLayer::Shape::EggCrate;
  p.depth_um = 0.25;
  p.below = "glass";
  p.above = "air";
  st.layers = {entry(p)};
  const auto rk = run(st, fg::SolverKind::RK, 3, 3, conical, 0, 1e-6);
  EXPECT_NEAR(rk.diag.energyBalance, 1.0, 5e-4);
  const auto rc = run(st, fg::SolverKind::RCWA, 3, 3, conical, 120);
  EXPECT_NEAR(rc.diag.energyBalance, 1.0, 5e-4);
  double dev = 0;
  for (const auto& a : rk.reflected)
    for (const auto& b : rc.reflected)
      if (a.n == b.n && a.m == b.m) dev = std::max(dev, std::abs(a.efficiency - b.efficiency));
  EXPECT_LT(dev, 2e-3);  // staircase error ~ 1/K
}

// Mirror x <-> z (non-square period, asymmetric inclusion, conical incidence): orders (n, m) and
// (m, n) carry the same power; Psi -> 90 deg - Psi, s and p stay s and p. Checks the 2D index
// handling of polygon coefficients, Algorithm NF and the Toeplitz blocks.
TEST(Polygon2D, MirrorXZ) {
  auto make = [](bool swap) {
    fg::Stack st;
    st.dx_um = swap ? 0.6 : 0.8;
    st.dz_um = swap ? 0.8 : 0.6;
    st.materials["air"] = fg::Material::fromIndex("air", 1.0);
    st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
    st.materials["si"] = fg::Material::fromIndex("si", cplx(3.8, 0.02));
    st.superstrate = "air";
    st.substrate = "glass";
    fg::Polygon2DLayer g;
    g.thickness_um = 0.15;
    g.background = "air";
    std::vector<std::array<double, 2>> v{{0.1, 0.05}, {0.55, 0.15}, {0.2, 0.45}};
    if (swap)
      for (auto& p : v) std::swap(p[0], p[1]);
    g.inclusions.push_back({v, "si"});
    st.layers = {entry(g)};
    return st;
  };
  for (double delta : {0.0, 90.0}) {
    const auto a = run(make(false), fg::SolverKind::RCWA, 3, 2, {0.7, 25 * DEG, 30 * DEG, delta * DEG});
    const auto b = run(make(true), fg::SolverKind::RCWA, 2, 3, {0.7, 25 * DEG, 60 * DEG, delta * DEG});
    double dev = 0;
    int matched = 0;
    using Side = std::vector<fg::OrderEfficiency> fg::Result::*;
    for (Side side : {&fg::Result::reflected, &fg::Result::transmitted})
      for (const auto& oa : a.*side)
        for (const auto& ob : b.*side)
          if (oa.n == ob.m && oa.m == ob.n) {
            dev = std::max(dev, std::abs(oa.efficiency - ob.efficiency));
            ++matched;
          }
    EXPECT_GE(matched, 4);
    EXPECT_LT(dev, 1e-10) << "delta=" << delta;
    EXPECT_NEAR(a.absorption, b.absorption, 1e-10);
  }
}

// ---------------------------------------------------------------- V10
namespace {
fg::Stack holeArray() {
  fg::Stack st;
  st.dx_um = 0.9;
  st.dz_um = 0.9;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["ag"] = fg::Material::fromIndex("ag", cplx(0.1, 8.94));
  st.superstrate = "air";
  st.substrate = "glass";
  fg::Polygon2DLayer g;
  g.thickness_um = 0.2;
  g.background = "ag";
  g.inclusions.push_back({{{0.325, 0.325}, {0.575, 0.325}, {0.575, 0.575}, {0.325, 0.575}}, "air"});
  st.layers = {entry(g)};
  return st;
}
double totalT(const fg::Result& r) {
  double t = 0;
  for (const auto& o : r.transmitted) t += o.efficiency;
  return t;
}
}  // namespace

// SPP pole of the flat Ag/glass interface for the (+-1, 0) orders: lambda = d Re sqrt(e_m e_d/(e_m+e_d)).
TEST(V10_HoleArray, PlasmonPoleOfFlatInterface) {
  const cplx em = cplx(0.1, 8.94) * cplx(0.1, 8.94), ed = 2.25;
  EXPECT_NEAR(0.9 * std::sqrt(em * ed / (em + ed)).real(), 1.369406, 1e-6);
}

// Normal incidence on square holes: C4 symmetry -> T independent of the polarization.
TEST(V10_HoleArray, PolarizationIndependence) {
  const auto st = holeArray();
  const double t0 = totalT(run(st, fg::SolverKind::RCWA, 4, 4, {1.388, 0.0, fg::PI, 0.0}));
  const double t1 = totalT(run(st, fg::SolverKind::RCWA, 4, 4, {1.388, 0.0, fg::PI, 90 * DEG}));
  EXPECT_NEAR(t0, t1, 1e-9 * std::max(t0, 1e-3));
}

// Extraordinary transmission: resonance peak just above the SPP pole, far above the flat-film
// value; reference T ~ 7.5 % (book) is reached only for N >= 9 (docs/validation.md, M7).
TEST(V10_HoleArray, TransmissionPeakAboveSPPPole) {
  const auto st = holeArray();
  double best = 0, lbest = 0;
  for (double l = 1.376; l <= 1.396; l += 0.002) {
    const auto r = run(st, fg::SolverKind::RCWA, 5, 5, {l, 0.0, fg::PI, 0.0});
    const double t = totalT(r);
    EXPECT_GT(r.absorption, 0.0) << l;
    if (t > best) { best = t; lbest = l; }
  }
  const double tOff = totalT(run(st, fg::SolverKind::RCWA, 5, 5, {1.50, 0.0, fg::PI, 0.0}));
  EXPECT_GT(lbest, 1.369406);
  EXPECT_LT(lbest, 1.395);
  EXPECT_GT(best, 0.02);
  EXPECT_GT(best, 10 * tOff);
}

// The example configuration parses (polygon2d) and reproduces the direct setup.
TEST(V10_HoleArray, ExampleConfiguration) {
  std::ifstream in(std::string(FG_SOURCE_DIR) + "/examples/v10_hole_array.json");
  std::stringstream ss;
  ss << in.rdbuf();
  auto cfg = fg::io::parseConfig(ss.str());
  cfg.solver.Nx = cfg.solver.Nz = 3;
  ASSERT_EQ(cfg.scan.size(), 1u);
  const double t1 = totalT(fg::makeSolver(cfg.solver)->solve(cfg.stack, cfg.scan[0]));
  const double t2 = totalT(run(holeArray(), fg::SolverKind::RCWA, 3, 3, {1.388, 0.0, fg::PI, 0.0}));
  EXPECT_NEAR(t1, t2, 1e-12);
}
