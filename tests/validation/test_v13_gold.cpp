// V13/V14: weakly absorbing metal (gold) - factorization variants (Plan section
// "Varianten der isotropen FFF"). Reference V13: C-method (coordinate transformation, no
// factorization problem), converged to 5 digits for N >= 5 (internal computation).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {
constexpr double DEG = fg::PI / 180.0;
constexpr double REF_EM1 = 0.26614, REF_E0 = 0.02362;

double eff(const std::vector<fg::OrderEfficiency>& v, int n) {
  for (const auto& o : v)
    if (o.n == n) return o.efficiency;
  return 0.0;
}
double sum(const std::vector<fg::OrderEfficiency>& v) {
  double s = 0;
  for (const auto& o : v) s += o.efficiency;
  return s;
}

fg::Stack gold(fg::Profile1DLayer p, fg::Factorization f) {
  fg::Stack st;
  st.dx_um = 1.0;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("au", cplx(0.2, 3.4));
  p.depth_um = 0.25;
  p.below = "substrate";
  p.above = "superstrate";
  fg::LayerEntry le;
  le.spec = p;
  le.factorization = f;
  st.layers.push_back(le);
  return st;
}

fg::Result solve(const fg::Stack& st, fg::SolverKind k, int N, int slices = 80) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rkTol = 1e-5;
  o.rcwaSlices = slices;
  return fg::makeSolver(o)->solve(st, {0.6328, 20 * DEG, 0, 90 * DEG});
}
}  // namespace

TEST(V13_GoldSinus, RigorousFFFMatchesCMethod) {
  const auto r = solve(gold({}, fg::Factorization::FFF_rig), fg::SolverKind::RK, 16);
  EXPECT_NEAR(eff(r.reflected, -1), REF_EM1, 1e-3);
  EXPECT_NEAR(eff(r.reflected, 0), REF_E0, 3e-4);
}

TEST(V13_GoldSinus, SymmetricFFFMatchesCMethod) {
  const auto r = solve(gold({}, fg::Factorization::FFF_sym), fg::SolverKind::RK, 16);
  EXPECT_NEAR(eff(r.reflected, -1), REF_EM1, 3e-3);
  EXPECT_NEAR(eff(r.reflected, 0), REF_E0, 1.5e-3);
}

TEST(V13_GoldSinus, BookFFFFailsAndWarns) {
  // Documents the failure of book eq. (IV.32) for weakly absorbing metals (Plan, Abschnitt
  // "Varianten der isotropen FFF"): broken +-lambda pairing creates unphysical modes.
  const auto r = solve(gold({}, fg::Factorization::FFF_iso), fg::SolverKind::RK, 16);
  EXPECT_GT(std::abs(eff(r.reflected, -1) - REF_EM1), 2e-2);
  EXPECT_FALSE(r.diag.warnings.empty());
}

TEST(V14_GoldTriangle, SymmetricFFFIsPhysical) {
  // Corners: the signed normal field jumps at apex and valley -> default FFF_sym. There is no
  // reliable reference yet (Laurent oscillates: 0.607/0.622/0.664/0.635 for N = 8/16/32/64), so
  // only physical bounds are checked; the case is listed as open in docs/validation.md.
  fg::Profile1DLayer p;
  p.shape = fg::Profile1DLayer::Shape::PiecewiseLinear;
  p.xs_um = {0.0, 0.7};
  p.ys_um = {0.0, 0.25};
  const auto st = gold(p, fg::Factorization::Auto);
  for (int N : {8, 16, 32}) {
    const auto r = solve(st, fg::SolverKind::RCWA, N);
    EXPECT_LT(sum(r.reflected), 1.0) << "N=" << N;
    EXPECT_GT(eff(r.reflected, -1), 0.0) << "N=" << N;
  }
}
