// Validation M8: anisotropic 2D gratings - V8 (2D) and consistency of the anisotropic 2D paths
// (Plan section "Validierung", docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/physics/anisotropic.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::Mat3C;

namespace {

constexpr double DEG = fg::PI / 180.0;
const fg::Wavelength conical{0.6328, 30 * DEG, 40 * DEG, 25 * DEG};

Mat3C cobalt() {
  const cplx e(-8.19, 16.38), xz(-0.495, -0.106);
  Mat3C E;
  E << e, 0, xz, 0, e, 0, -xz, 0, e;
  return E;
}

fg::LayerEntry entry(fg::LayerSpec s, fg::Factorization f = fg::Factorization::Auto) {
  fg::LayerEntry e;
  e.spec = std::move(s);
  e.factorization = f;
  return e;
}

fg::Stack baseStack(double dz) {
  fg::Stack st;
  st.dx_um = 0.8;
  st.dz_um = dz;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.superstrate = "air";
  st.substrate = "glass";
  return st;
}

fg::Result run(const fg::Stack& st, fg::SolverKind k, int Nx, int Nz, const fg::Wavelength& wl = conical,
               int slices = 12, double tol = 1e-6) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = Nx;
  o.Nz = Nz;
  o.rcwaSlices = slices;
  o.rkTol = tol;
  return fg::makeSolver(o)->solve(st, wl);
}

/// Largest efficiency difference of matching orders; 1D orders n are matched to (n, 0), all other
/// 2D orders must be power-free (compareOrders = false: same order sets expected).
double compare(const fg::Result& a, const fg::Result& b) {
  double dev = 0;
  auto side = [&](const std::vector<fg::OrderEfficiency>& x, const std::vector<fg::OrderEfficiency>& y) {
    for (const auto& o1 : x) {
      bool found = false;
      for (const auto& o2 : y)
        if (o2.n == o1.n && o2.m == o1.m) {
          dev = std::max(dev, std::abs(o1.efficiency - o2.efficiency));
          found = true;
        }
      if (!found) dev = std::max(dev, std::abs(o1.efficiency));
    }
    for (const auto& o2 : y) {
      bool found = false;
      for (const auto& o1 : x) found = found || (o1.n == o2.n && o1.m == o2.m);
      if (!found) dev = std::max(dev, std::abs(o2.efficiency));
    }
  };
  side(a.reflected, b.reflected);
  side(a.transmitted, b.transmitted);
  return std::max(dev, std::abs(a.absorption - b.absorption));
}

fg::Polygon2DLayer square(const std::string& inc, const std::string& bg, double dz) {
  fg::Polygon2DLayer g;
  g.thickness_um = 0.12;
  g.background = bg;
  g.inclusions.push_back({{{0.15, 0.1}, {0.55, 0.15}, {0.5, 0.45 * dz / 0.6}, {0.2, 0.5 * dz / 0.6}}, inc});
  return g;
}

fg::Polygon2DLayer strip(const std::string& m, double dz) {
  fg::Polygon2DLayer p;
  p.thickness_um = 0.35;
  p.background = "air";
  p.inclusions.push_back({{{0.0, 0.0}, {0.3, 0.0}, {0.3, dz}, {0.0, dz}}, m});
  return p;
}

}  // namespace

// ---------------------------------------------------------------- V8 (2D)
// Laurent with an isotropic tensor = isotropic Laurent (polygon cross-section and egg crate).
TEST(V8_2D, LaurentTensorEqualsIsotropic) {
  const cplx n(0.5, 3.0);
  auto make = [&](bool tensor, bool profile) {
    fg::Stack st = baseStack(0.6);
    st.materials["m"] = tensor ? fg::Material::fromTensor("m", n * n * Mat3C::Identity()) : fg::Material::fromIndex("m", n);
    if (profile) {
      fg::Profile2DLayer p;
      p.shape = fg::Profile2DLayer::Shape::EggCrate;
      p.depth_um = 0.2;
      p.below = "m";
      p.above = "air";
      st.layers = {entry(p, fg::Factorization::Laurent)};
    } else {
      st.layers = {entry(square("m", "air", 0.6), fg::Factorization::Laurent)};
    }
    return st;
  };
  for (bool profile : {false, true}) {
    const auto a = run(make(false, profile), fg::SolverKind::RCWA, 3, 2);
    const auto b = run(make(true, profile), fg::SolverKind::RCWA, 3, 2);
    EXPECT_LT(compare(a, b), 1e-10) << "profile=" << profile;
  }
}

// FFF_B on a z-invariant polygon strip: rank-1 field = x exactly, so variant B with an isotropic
// tensor is the Li lamellar formula (Plan eq:QaniB with constant P), as in 1D.
TEST(V8_2D, ProjectorVariantOnStripEqualsLi) {
  const cplx n(1.38, 7.62);
  auto s1 = baseStack(0.0);
  s1.materials["m"] = fg::Material::fromIndex("m", n);
  s1.layers = {entry(fg::Lamellar1DLayer{0.35, {0.0, 0.3, 0.8}, {"m", "air"}})};
  const auto r1 = run(s1, fg::SolverKind::RCWA, 8, 0);
  for (bool tensor : {false, true})
    for (int Nz : {0, 2}) {
      auto s2 = baseStack(0.7);
      s2.materials["m"] = tensor ? fg::Material::fromTensor("m", n * n * Mat3C::Identity()) : fg::Material::fromIndex("m", n);
      s2.layers = {entry(strip("m", 0.7), fg::Factorization::FFF_B)};
      const auto r2 = run(s2, fg::SolverKind::RCWA, 8, Nz);
      double dev = 0;
      for (const auto* v : {&r1.reflected, &r1.transmitted})
        for (const auto& o1 : *v)
          for (const auto* w : {&r2.reflected, &r2.transmitted})
            for (const auto& o2 : *w)
              if (o2.n == o1.n && o2.m == 0 && (v == &r1.reflected) == (w == &r2.reflected))
                dev = std::max(dev, std::abs(o1.efficiency - o2.efficiency));
      EXPECT_LT(dev, 1e-10) << "tensor=" << tensor << " Nz=" << Nz;
      EXPECT_NEAR(r1.absorption, r2.absorption, 1e-10);
    }
}

// Anisotropic z-invariant structures in 2D stacks = 1D anisotropic code: cobalt strip (FFF_B, Laurent)
// and anisotropic sinus_x profile (FFF_A, FFF_B; RK and RCWA).
TEST(V8_2D, ZInvariantAnisotropicEquals1D) {
  auto mk = [](double dz) {
    fg::Stack st = baseStack(dz);
    st.materials["co"] = fg::Material::fromTensor("co", cobalt());
    return st;
  };
  for (auto f : {fg::Factorization::FFF_B, fg::Factorization::Laurent}) {
    auto s1 = mk(0.0), s2 = mk(0.7);
    s1.layers = {entry(fg::Lamellar1DLayer{0.35, {0.0, 0.3, 0.8}, {"co", "air"}}, f)};
    s2.layers = {entry(strip("co", 0.7), f)};
    EXPECT_LT(compare(run(s1, fg::SolverKind::RCWA, 6, 0), run(s2, fg::SolverKind::RCWA, 6, 0)), 1e-10);
  }
  for (auto f : {fg::Factorization::FFF_A, fg::Factorization::FFF_B}) {
    auto s1 = mk(0.0), s2 = mk(0.6);
    fg::Profile1DLayer p1;
    p1.depth_um = 0.2;
    p1.below = "co";
    p1.above = "air";
    fg::Profile2DLayer p2;
    p2.shape = fg::Profile2DLayer::Shape::SinusX;
    p2.depth_um = 0.2;
    p2.below = "co";
    p2.above = "air";
    s1.layers = {entry(p1, f)};
    s2.layers = {entry(p2, f)};
    for (auto k : {fg::SolverKind::RCWA, fg::SolverKind::RK})
      // 1D code against the z-invariant 2D code: different code paths, agreement at rounding level
      // (1.6e-12 without FMA contraction in normal_field.cpp, below 1e-12 with it; both are rounding)
      EXPECT_LT(compare(run(s1, k, 5, 0), run(s2, k, 5, 0)), 5e-12);
    // 1D anisotropic profile inside a 2D stack (z-invariant)
    auto s3 = mk(0.6);
    s3.layers = {entry(p1, f)};
    EXPECT_LT(compare(run(s1, fg::SolverKind::RCWA, 5, 0), run(s3, fg::SolverKind::RCWA, 5, 0)), 1e-12);
  }
}

// x <-> z mirror of an anisotropic cross-section (cobalt quadrilateral, FFF_B): the tensor is mirrored
// as well (eps' = S eps S, S swaps x and z). Checks the 2D index handling of the projector path.
TEST(V8_2D, MirrorXZAnisotropic) {
  auto make = [](bool swap) {
    fg::Stack st;
    st.dx_um = swap ? 0.6 : 0.8;
    st.dz_um = swap ? 0.8 : 0.6;
    st.materials["air"] = fg::Material::fromIndex("air", 1.0);
    st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
    Mat3C S = Mat3C::Zero();
    S(0, 2) = S(2, 0) = S(1, 1) = 1.0;
    st.materials["co"] = fg::Material::fromTensor("co", swap ? Mat3C(S * cobalt() * S) : cobalt());
    st.superstrate = "air";
    st.substrate = "glass";
    std::vector<std::array<double, 2>> v{{0.1, 0.05}, {0.55, 0.15}, {0.2, 0.45}};
    if (swap)
      for (auto& p : v) std::swap(p[0], p[1]);
    fg::Polygon2DLayer g;
    g.thickness_um = 0.1;
    g.background = "air";
    g.inclusions.push_back({v, "co"});
    g.normalField.grid = 256;
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
    EXPECT_LT(dev, 1e-9) << "delta=" << delta;
  }
}

// ---------------------------------------------------------------- convergence (M8)
// Cobalt pillars (square, 0.3 um, period 0.6 um, height 0.1 um) on glass, conical incidence:
// FFF_B and Laurent converge to the same efficiencies; FFF_B converges faster (docs/validation.md).
TEST(M8_Convergence, CobaltPillarsLaurentAndProjector) {
  auto make = [](fg::Factorization f) {
    fg::Stack st;
    st.dx_um = st.dz_um = 0.6;
    st.materials["air"] = fg::Material::fromIndex("air", 1.0);
    st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
    st.materials["co"] = fg::Material::fromTensor("co", cobalt());
    st.superstrate = "air";
    st.substrate = "glass";
    fg::Polygon2DLayer g;
    g.thickness_um = 0.1;
    g.background = "air";
    g.inclusions.push_back({{{0.15, 0.15}, {0.45, 0.15}, {0.45, 0.45}, {0.15, 0.45}}, "co"});
    st.layers = {entry(g, f)};
    return st;
  };
  const fg::Wavelength wl{0.6328, 30 * DEG, 20 * DEG, 60 * DEG};
  auto T = [&](fg::Factorization f, int N) {
    double t = 0;
    for (const auto& o : run(make(f), fg::SolverKind::RCWA, N, N, wl).transmitted) t += o.efficiency;
    return t;
  };
  // N = 2..7: FFF_B 0.548 0.543 0.549 0.549 0.549 0.551; Laurent 0.417 0.460 0.474 0.477 0.496 0.491
  const double b3 = T(fg::Factorization::FFF_B, 3), b5 = T(fg::Factorization::FFF_B, 5), b6 = T(fg::Factorization::FFF_B, 6);
  const double l3 = T(fg::Factorization::Laurent, 3), l6 = T(fg::Factorization::Laurent, 6);
  EXPECT_LT(std::abs(b6 - b5), 2e-3);
  EXPECT_LT(std::abs(l6 - b6), 0.7 * std::abs(l3 - b3));  // Laurent approaches FFF_B slowly
  EXPECT_NEAR(b6, 0.549, 3e-3);
}

// Anisotropic egg crate (rotated uniaxial tensor, lossless): FFF_A, FFF_B and Laurent agree, RK = RCWA,
// energy within the FFF tolerance.
TEST(M8_Convergence, AnisotropicEggCrate) {
  auto make = [](fg::Factorization f) {
    fg::Stack st;
    st.dx_um = st.dz_um = 0.6;
    st.materials["air"] = fg::Material::fromIndex("air", 1.0);
    st.materials["cal"] = fg::Material::fromTensor("cal", fg::rotatedTensor(2.75, 2.75, 2.21, 30 * DEG, 40 * DEG, 0));
    st.superstrate = "air";
    st.substrate = "cal";
    fg::Profile2DLayer p;
    p.shape = fg::Profile2DLayer::Shape::EggCrate;
    p.depth_um = 0.25;
    p.below = "cal";
    p.above = "air";
    st.layers = {entry(p, f)};
    return st;
  };
  const fg::Wavelength wl{0.6328, 30 * DEG, 20 * DEG, 60 * DEG};
  // r00 for N = 2..6 (20 slices): FFF_A 0.00610 0.00603 0.00602 0.00601 0.00600,
  // FFF_B 0.00609 0.00609 0.00602 0.00601, Laurent 0.00561 0.00569 0.00574 0.00578 (slower)
  const auto a = run(make(fg::Factorization::FFF_A), fg::SolverKind::RCWA, 3, 3, wl, 20);
  const auto b = run(make(fg::Factorization::FFF_B), fg::SolverKind::RCWA, 3, 3, wl, 20);
  const auto l = run(make(fg::Factorization::Laurent), fg::SolverKind::RCWA, 3, 3, wl, 20);
  EXPECT_LT(compare(a, b), 1e-3);
  EXPECT_LT(compare(a, l), 6e-3);
  EXPECT_NEAR(a.diag.energyBalance, 1.0, 1e-3);
  const auto rk = run(make(fg::Factorization::FFF_A), fg::SolverKind::RK, 2, 2, wl, 0, 1e-6);
  const auto rc = run(make(fg::Factorization::FFF_A), fg::SolverKind::RCWA, 2, 2, wl, 60);
  EXPECT_LT(compare(rk, rc), 2e-3);
}
