// Validation cases V4-V6 and conical profile checks (Plan section "Validierung", docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {
constexpr double DEG = fg::PI / 180.0;

double eff(const std::vector<fg::OrderEfficiency>& v, int n) {
  for (const auto& o : v)
    if (o.n == n) return o.efficiency;
  return 0.0;
}

fg::Stack sinusStack(double d, double a, cplx nSub, double x0 = 0.0,
                     fg::Factorization f = fg::Factorization::Auto) {
  fg::Stack st;
  st.dx_um = d;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("sub", nSub);
  fg::Profile1DLayer p;
  p.depth_um = a;
  p.x0_um = x0;
  p.below = "substrate";
  p.above = "superstrate";
  fg::LayerEntry le;
  le.spec = p;
  le.factorization = f;
  st.layers.push_back(le);
  return st;
}

fg::Result run(const fg::Stack& st, fg::SolverKind k, int N, const fg::Wavelength& wl, int slices = 120,
               double tol = 1e-5) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rkTol = tol;
  o.rcwaSlices = slices;
  return fg::makeSolver(o)->solve(st, wl);
}
}  // namespace

// V4: Al sinus grating, reference from the C-method (internal): TE 0.3235, TM 0.8696
TEST(V4_AlSinus, DifferentialMethodRK) {
  const auto st = sinusStack(0.5, 0.2, cplx(1.38, 7.62));
  const auto te = run(st, fg::SolverKind::RK, 20, {0.6328, 30 * DEG, 0, 0});
  const auto tm = run(st, fg::SolverKind::RK, 20, {0.6328, 30 * DEG, 0, 90 * DEG});
  EXPECT_NEAR(eff(te.reflected, -1), 0.3235, 2e-3);
  EXPECT_NEAR(eff(tm.reflected, -1), 0.8696, 2e-3);
}

TEST(V4_AlSinus, RCWAStaircaseWithFFF) {
  // staircase with the smooth normal field converges like 1/K towards the RK value
  const auto st = sinusStack(0.5, 0.2, cplx(1.38, 7.62));
  const fg::Wavelength tm{0.6328, 30 * DEG, 0, 90 * DEG};
  const double e100 = eff(run(st, fg::SolverKind::RCWA, 20, tm, 100).reflected, -1);
  const double e300 = eff(run(st, fg::SolverKind::RCWA, 20, tm, 300).reflected, -1);
  EXPECT_NEAR(e300, 0.8696, 2e-3);
  EXPECT_LT(std::abs(e300 - 0.8696), std::abs(e100 - 0.8696));
  const double te = eff(run(st, fg::SolverKind::RCWA, 20, {0.6328, 30 * DEG, 0, 0}, 100).reflected, -1);
  EXPECT_NEAR(te, 0.3235, 2e-3);
}

// V5: deep glass grating (depth/period = 2) in transmission, order -1, TE: 0.820 (book Fig. III.4)
TEST(V5_GlassDeep, TransmissionTE) {
  const auto st = sinusStack(0.5, 1.0, 1.5, 0.25);
  const fg::Wavelength te{0.6328, 30 * DEG, 0, 0};
  const auto a = run(st, fg::SolverKind::RCWA, 10, te, 120);
  const auto b = run(st, fg::SolverKind::RK, 10, te, 0, 1e-6);
  EXPECT_NEAR(eff(a.transmitted, -1), 0.820, 5e-3);
  EXPECT_NEAR(eff(b.transmitted, -1), 0.820, 5e-3);
  EXPECT_NEAR(a.diag.energyBalance, 1.0, 1e-10);  // TE: Laurent-type factorization, exact
  EXPECT_NEAR(b.diag.energyBalance, 1.0, 1e-7);
}

TEST(V5_GlassDeep, EnergyBalanceTM) {
  // FFF is not Hermitian in the truncated space: energy only approximately conserved
  const auto st = sinusStack(0.5, 1.0, 1.5, 0.25);
  const auto r = run(st, fg::SolverKind::RCWA, 15, {0.6328, 30 * DEG, 0, 90 * DEG}, 120);
  EXPECT_NEAR(r.diag.energyBalance, 1.0, 5e-4);
}

// V6: dielectric sinus nu = 2.5, a/d = 1, TM: FFF converges much faster than Laurent
TEST(V6_DielectricSinus, FFFConvergesFasterThanLaurent) {
  const fg::Wavelength tm{0.6328, 30 * DEG, 0, 90 * DEG};
  auto t = [&](fg::Factorization f, int N) {
    return eff(run(sinusStack(1.0, 1.0, 2.5, 0.0, f), fg::SolverKind::RCWA, N, tm, 120).transmitted, -1);
  };
  const double dF = std::abs(t(fg::Factorization::Auto, 8) - t(fg::Factorization::Auto, 25));
  const double dL = std::abs(t(fg::Factorization::Laurent, 8) - t(fg::Factorization::Laurent, 25));
  EXPECT_LT(dF, 5e-4);
  EXPECT_LT(dF, 0.2 * dL);
}

// M4: conical incidence on a profile (full 4 n_f system, TE/TM coupled)
TEST(M4_ConicalProfile, RKAgreesWithRCWAAndEnergy) {
  const auto st = sinusStack(0.5, 0.25, 1.5, 0.0);
  const fg::Wavelength wl{0.6328, 30 * DEG, 40 * DEG, 30 * DEG};
  const auto a = run(st, fg::SolverKind::RK, 8, wl, 0, 1e-6);
  const auto b = run(st, fg::SolverKind::RCWA, 8, wl, 400);
  EXPECT_NEAR(a.diag.energyBalance, 1.0, 5e-4);
  for (int n : {-1, 0})
    EXPECT_NEAR(eff(a.reflected, n), eff(b.reflected, n), 1e-4) << "order " << n;
  for (int n : {-1, 0})
    EXPECT_NEAR(eff(a.transmitted, n), eff(b.transmitted, n), 1e-4) << "order " << n;
}
