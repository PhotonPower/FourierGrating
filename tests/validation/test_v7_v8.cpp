// Validation M5: anisotropic 1D gratings (Plan section "Validierung", docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/physics/anisotropic.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::Mat3C;

namespace {
constexpr double DEG = fg::PI / 180.0;

const fg::OrderEfficiency* order(const std::vector<fg::OrderEfficiency>& v, int n) {
  for (const auto& o : v)
    if (o.n == n) return &o;
  return nullptr;
}
double sum(const std::vector<fg::OrderEfficiency>& v) {
  double s = 0;
  for (const auto& o : v) s += o.efficiency;
  return s;
}

Mat3C cobalt() {
  const cplx e(-8.19, 16.38), xz(-0.495, -0.106);
  Mat3C E;
  E << e, 0, xz, 0, e, 0, -xz, 0, e;
  return E;
}

fg::Stack cobaltStack(fg::Factorization f) {
  fg::Stack st;
  st.dx_um = 0.6;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromTensor("co", cobalt());
  fg::Profile1DLayer p;
  p.depth_um = 0.5;
  p.below = "substrate";
  p.above = "superstrate";
  fg::LayerEntry le;
  le.spec = p;
  le.factorization = f;
  st.layers.push_back(le);
  return st;
}

fg::Result run(const fg::Stack& st, fg::SolverKind k, int N, const fg::Wavelength& wl, int slices) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rcwaSlices = slices;
  o.rkTol = 1e-6;
  return fg::makeSolver(o)->solve(st, wl);
}

const fg::Wavelength TM30{0.6328, 30 * DEG, 0, 90 * DEG};
}  // namespace

// V7: cobalt sinus grating (book chapter VIII), TM incidence; internal reference N = 20, 100 slices
TEST(V7_Cobalt, VariantAMatchesReference) {
  const auto r = run(cobaltStack(fg::Factorization::FFF_A), fg::SolverKind::RCWA, 20, TM30, 100);
  const auto* o0 = order(r.reflected, 0);
  const auto* o1 = order(r.reflected, -1);
  ASSERT_TRUE(o0 && o1);
  EXPECT_NEAR(o0->eff_p, 0.0567, 0.03 * 0.0567);
  EXPECT_NEAR(o1->eff_p, 0.1435, 0.03 * 0.1435);
  EXPECT_NEAR(o0->eff_s, 1.3e-7, 0.1 * 1.3e-7);
  EXPECT_NEAR(o1->eff_s, 2.37e-6, 0.05 * 2.37e-6);
}

TEST(V7_Cobalt, VariantBAgreesWithVariantA) {
  const auto a = run(cobaltStack(fg::Factorization::FFF_A), fg::SolverKind::RCWA, 12, TM30, 60);
  const auto b = run(cobaltStack(fg::Factorization::FFF_B), fg::SolverKind::RCWA, 12, TM30, 60);
  for (int n : {0, -1}) {
    EXPECT_NEAR(order(b.reflected, n)->eff_p, order(a.reflected, n)->eff_p, 0.03 * order(a.reflected, n)->eff_p);
    EXPECT_NEAR(order(b.reflected, n)->eff_s, order(a.reflected, n)->eff_s, 0.15 * order(a.reflected, n)->eff_s);
  }
}

TEST(V7_Cobalt, RKAgreesWithRCWA) {
  const auto st = cobaltStack(fg::Factorization::FFF_A);
  const auto a = run(st, fg::SolverKind::RK, 8, TM30, 0);
  const auto b = run(st, fg::SolverKind::RCWA, 8, TM30, 400);
  for (int n : {0, -1}) {
    EXPECT_NEAR(order(a.reflected, n)->eff_p, order(b.reflected, n)->eff_p, 3e-4);
    EXPECT_NEAR(order(a.reflected, n)->eff_s, order(b.reflected, n)->eff_s, 0.02 * order(b.reflected, n)->eff_s);
  }
}

// V8: anisotropic code with an isotropic tensor
TEST(V8_IsotropicTensor, LamellarEqualsIsotropicCode) {
  auto make = [](bool asTensor) {
    fg::Stack st;
    st.dx_um = 1.0;
    st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
    st.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
    const cplx e = cplx(0.5, 3.0) * cplx(0.5, 3.0);
    st.materials["bar"] = asTensor ? fg::Material::fromTensor("bar", e * Mat3C::Identity())
                                   : fg::Material::fromIndex("bar", cplx(0.5, 3.0));
    fg::LayerEntry le;
    le.spec = fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"bar", "superstrate"}};
    st.layers.push_back(le);
    return st;
  };
  for (const fg::Wavelength wl : {fg::Wavelength{0.6328, 20 * DEG, 0, 90 * DEG},
                                  fg::Wavelength{0.6328, 25 * DEG, 30 * DEG, 45 * DEG}}) {
    const auto a = run(make(false), fg::SolverKind::RCWA, 12, wl, 1);
    const auto b = run(make(true), fg::SolverKind::RCWA, 12, wl, 1);
    ASSERT_EQ(a.reflected.size(), b.reflected.size());
    for (std::size_t k = 0; k < a.reflected.size(); ++k)
      EXPECT_NEAR(a.reflected[k].efficiency, b.reflected[k].efficiency, 1e-10);
  }
}

TEST(V8_IsotropicTensor, ProfileAgreesWithIsotropicFFF) {
  // Al sinus (V4) as a tensor: FFF_A vs isotropic FFF_rig, same staircase (100 slices), N = 20
  // (both converge to the same value; at N = 15 they still differ by 4e-3)
  auto make = [](bool asTensor) {
    fg::Stack st;
    st.dx_um = 0.5;
    st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
    const cplx n(1.38, 7.62);
    st.materials["substrate"] = asTensor ? fg::Material::fromTensor("al", n * n * Mat3C::Identity())
                                         : fg::Material::fromIndex("al", n);
    fg::Profile1DLayer p;
    p.depth_um = 0.2;
    p.below = "substrate";
    p.above = "superstrate";
    fg::LayerEntry le;
    le.spec = p;
    st.layers.push_back(le);
    return st;
  };
  const auto a = run(make(false), fg::SolverKind::RCWA, 20, TM30, 100);
  const auto b = run(make(true), fg::SolverKind::RCWA, 20, TM30, 100);
  EXPECT_NEAR(order(a.reflected, -1)->efficiency, order(b.reflected, -1)->efficiency, 1e-3);
}

// Anisotropic homogeneous media: analytic uniaxial substrate, lossless biaxial slab
TEST(M5_AnisotropicMedia, UniaxialSubstrateAnalytic) {
  const double eo = 2.4, ee = 3.1;  // optic axis along y (normal to the interface)
  Mat3C e = Mat3C::Zero();
  e(0, 0) = eo;
  e(1, 1) = ee;
  e(2, 2) = eo;
  fg::Stack st;
  st.dx_um = 0.7;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromTensor("uni", e);
  for (double phi : {0.0, 30.0, 70.0}) {
    const double a = std::sin(phi * DEG), b0 = std::cos(phi * DEG);
    const double kyo = std::sqrt(eo - a * a), kye = std::sqrt(eo * (1 - a * a / ee));
    const double rs = (b0 - kyo) / (b0 + kyo), rp = (b0 - kye / eo) / (b0 + kye / eo);
    const auto s = run(st, fg::SolverKind::RCWA, 2, {0.6328, phi * DEG, 0, 0}, 1);
    const auto p = run(st, fg::SolverKind::RCWA, 2, {0.6328, phi * DEG, 0, 90 * DEG}, 1);
    EXPECT_NEAR(sum(s.reflected), rs * rs, 1e-12) << "phi=" << phi;
    EXPECT_NEAR(sum(p.reflected), rp * rp, 1e-12) << "phi=" << phi;
    EXPECT_NEAR(s.diag.energyBalance, 1.0, 1e-12);
    EXPECT_NEAR(p.diag.energyBalance, 1.0, 1e-12);
  }
}

TEST(M5_AnisotropicMedia, LosslessBiaxialSlabAllPaths) {
  // rotated biaxial lossless slab on glass, conical incidence: energy conservation, and the
  // homogeneous-layer path equals the 'lamellar of one material' path (RCWA and RK, FFF_A)
  const Mat3C e = fg::rotatedTensor(2.1, 2.6, 3.3, 0.4, 0.9, -0.3);
  auto make = [&](bool lamellar) {
    fg::Stack st;
    st.dx_um = 0.8;
    st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
    st.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
    st.materials["bi"] = fg::Material::fromTensor("bi", e);
    fg::LayerEntry le;
    if (lamellar)
      le.spec = fg::Lamellar1DLayer{0.35, {0.0, 0.3, 0.8}, {"bi", "bi"}};
    else
      le.spec = fg::HomogeneousLayer{0.35, "bi"};
    st.layers.push_back(le);
    return st;
  };
  const fg::Wavelength wl{0.6328, 35 * DEG, 25 * DEG, 30 * DEG};
  const auto h = run(make(false), fg::SolverKind::RCWA, 3, wl, 1);
  const auto l = run(make(true), fg::SolverKind::RCWA, 3, wl, 1);
  const auto k = run(make(true), fg::SolverKind::RK, 3, wl, 1);
  EXPECT_NEAR(h.diag.energyBalance, 1.0, 1e-12);
  EXPECT_NEAR(sum(l.reflected), sum(h.reflected), 1e-10);
  EXPECT_NEAR(sum(k.reflected), sum(h.reflected), 1e-8);
  EXPECT_NEAR(sum(k.transmitted), sum(h.transmitted), 1e-8);
}
