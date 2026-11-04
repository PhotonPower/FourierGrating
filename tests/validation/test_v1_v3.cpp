// Validation cases V1-V3 (Plan section "Validierung und Tests", docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {

constexpr double DEG = fg::PI / 180.0;

cplx branch(cplx z) {
  cplx b = std::sqrt(z);
  if (b.imag() < 0 || (b.imag() == 0 && b.real() < 0)) b = -b;
  return b;
}

/// Reflectance of a flat stack: superstrate n0 | optional layer (n1, t) | substrate n2, s and p.
std::pair<double, double> airy(cplx n0, cplx n1, double t_um, cplx n2, double theta, double lambda) {
  const cplx s2 = n0 * n0 * std::sin(theta) * std::sin(theta);
  const cplx q0 = branch(n0 * n0 - s2), q1 = branch(n1 * n1 - s2), q2 = branch(n2 * n2 - s2);
  auto rs = [](cplx qa, cplx qb) { return (qa - qb) / (qa + qb); };
  auto rp = [](cplx na, cplx qa, cplx nb, cplx qb) {
    return (nb * nb * qa - na * na * qb) / (nb * nb * qa + na * na * qb);
  };
  const cplx ph = std::exp(fg::I * (2.0 * q1 * 2.0 * fg::PI / lambda * t_um));
  const cplx r01s = rs(q0, q1), r12s = rs(q1, q2);
  const cplx r01p = rp(n0, q0, n1, q1), r12p = rp(n1, q1, n2, q2);
  const cplx Rs = (r01s + r12s * ph) / (1.0 + r01s * r12s * ph);
  const cplx Rp = (r01p + r12p * ph) / (1.0 + r01p * r12p * ph);
  return {std::norm(Rs), std::norm(Rp)};
}

fg::Stack flatStack(cplx nSub) {
  fg::Stack st;
  st.dx_um = 0.7;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("sub", nSub);
  return st;
}

double sumEff(const std::vector<fg::OrderEfficiency>& v) {
  double s = 0;
  for (const auto& o : v) s += o.efficiency;
  return s;
}

double eff(const std::vector<fg::OrderEfficiency>& v, int n) {
  for (const auto& o : v)
    if (o.n == n) return o.efficiency;
  return 0.0;
}

fg::SolverOptions opts(fg::SolverKind k, int N) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rkTol = 1e-11;
  return o;
}

}  // namespace

// ---------------------------------------------------------------- V1: Fresnel
TEST(V1_Fresnel, AllAnglesAzimuthsPolarizations) {
  for (cplx nSub : {cplx(1.5, 0), cplx(1.38, 7.62), cplx(0.2, 3.4)}) {
    const auto st = flatStack(nSub);
    for (double phi : {0.0, 25.0, 60.0, 85.0})
      for (double psi : {0.0, 30.0, 180.0})
        for (double delta : {0.0, 40.0, 90.0}) {
          const auto [Rs, Rp] = airy(1.0, 1.0, 0.0, nSub, phi * DEG, 0.6328);
          const double cd = std::cos(delta * DEG), sd = std::sin(delta * DEG);
          const double Rref = cd * cd * Rs + sd * sd * Rp;
          const auto res = fg::makeSolver(opts(fg::SolverKind::RCWA, 2))
                               ->solve(st, {0.6328, phi * DEG, psi * DEG, delta * DEG});
          EXPECT_NEAR(sumEff(res.reflected), Rref, 1e-12)
              << "n=" << nSub << " phi=" << phi << " psi=" << psi << " delta=" << delta;
          if (nSub.imag() == 0) EXPECT_NEAR(res.diag.energyBalance, 1.0, 1e-12);
        }
  }
}

// ---------------------------------------------------------------- V2: Airy
TEST(V2_Airy, HomogeneousLayerBothSolvers) {
  for (cplx nLay : {cplx(2.1, 0), cplx(0.6, 2.5)}) {
    auto st = flatStack(1.5);
    st.materials["layer"] = fg::Material::fromIndex("layer", nLay);
    fg::LayerEntry le;
    le.spec = fg::HomogeneousLayer{0.37, "layer"};
    st.layers.push_back(le);
    // the same layer as a "lamellar" layer of one material exercises the RK/RCWA slice paths
    auto st2 = st;
    fg::Lamellar1DLayer lam{0.37, {0.0, 0.3, 0.7}, {"layer", "layer"}};
    st2.layers[0].spec = lam;
    for (double phi : {0.0, 35.0, 70.0})
      for (double delta : {0.0, 90.0}) {
        const auto [Rs, Rp] = airy(1.0, nLay, 0.37, 1.5, phi * DEG, 0.6328);
        const double Rref = delta == 0.0 ? Rs : Rp;
        const fg::Wavelength wl{0.6328, phi * DEG, 20 * DEG, delta * DEG};
        // conical azimuth (psi = 20 deg) keeps s/p decoupled for flat layers
        const double R1 = sumEff(fg::makeSolver(opts(fg::SolverKind::RCWA, 3))->solve(st, wl).reflected);
        const double R2 = sumEff(fg::makeSolver(opts(fg::SolverKind::RCWA, 3))->solve(st2, wl).reflected);
        const double R3 = sumEff(fg::makeSolver(opts(fg::SolverKind::RK, 3))->solve(st2, wl).reflected);
        EXPECT_NEAR(R1, Rref, 1e-12) << nLay << " phi=" << phi << " delta=" << delta;
        EXPECT_NEAR(R2, Rref, 1e-10);
        EXPECT_NEAR(R3, Rref, 1e-8);
      }
  }
}

// ---------------------------------------------------------------- V3: lamellar gratings
namespace {
fg::Stack lamellarStack(cplx nBar) {
  fg::Stack st;
  st.dx_um = 1.0;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["bar"] = fg::Material::fromIndex("bar", nBar);
  fg::LayerEntry le;
  le.spec = fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"bar", "superstrate"}};
  st.layers.push_back(le);
  return st;
}
}  // namespace

TEST(V3_Lamellar, EnergyConservationLossless) {
  const auto st = lamellarStack(2.0);
  for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK})
    for (const fg::Wavelength wl : {fg::Wavelength{0.6328, 20 * DEG, 0, 0},
                                    fg::Wavelength{0.6328, 20 * DEG, 0, 90 * DEG},
                                    fg::Wavelength{0.6328, 35 * DEG, 30 * DEG, 45 * DEG}}) {
      const auto res = fg::makeSolver(opts(kind, 15))->solve(st, wl);
      EXPECT_NEAR(res.diag.energyBalance, 1.0, kind == fg::SolverKind::RCWA ? 1e-10 : 1e-8);
      EXPECT_GE(res.reflected.size(), 2u);
    }
}

TEST(V3_Lamellar, RKEqualsRCWA) {
  for (cplx nBar : {cplx(2.0, 0), cplx(0.5, 3.0)}) {
    const auto st = lamellarStack(nBar);
    for (double delta : {0.0, 90.0}) {
      const fg::Wavelength wl{0.6328, 20 * DEG, 15 * DEG, delta * DEG};
      const auto a = fg::makeSolver(opts(fg::SolverKind::RCWA, 12))->solve(st, wl);
      const auto b = fg::makeSolver(opts(fg::SolverKind::RK, 12))->solve(st, wl);
      ASSERT_EQ(a.reflected.size(), b.reflected.size());
      for (std::size_t k = 0; k < a.reflected.size(); ++k)
        EXPECT_NEAR(a.reflected[k].efficiency, b.reflected[k].efficiency, 1e-8)
            << "n=" << nBar << " delta=" << delta << " order " << a.reflected[k].n;
    }
  }
}

TEST(V3_Lamellar, ReciprocityInReflection) {
  // e_n(alpha0) = e_n(-alpha_n) holds for the exact problem. With a symmetric truncation
  // |n| <= N around each incidence the two retained order sets differ by one order, so the
  // FMM satisfies reciprocity only up to truncation: the deviation must decrease with N.
  for (cplx nBar : {cplx(2.0, 0), cplx(0.5, 3.0)}) {
    const auto st = lamellarStack(nBar);
    const double lambda = 0.6328, theta = 20 * DEG;
    const double K = lambda / st.dx_um;  // normalized grating vector (nu_sup = 1)
    const int n = -1;
    const double thetaRev = std::asin(-(std::sin(theta) + n * K));
    for (double delta : {0.0, 90.0}) {
      double dev[2];
      const int Ns[2] = {10, 40};
      for (int k = 0; k < 2; ++k) {
        const auto a = fg::makeSolver(opts(fg::SolverKind::RCWA, Ns[k]))->solve(st, {lambda, theta, 0, delta * DEG});
        const auto b = fg::makeSolver(opts(fg::SolverKind::RCWA, Ns[k]))->solve(st, {lambda, thetaRev, 0, delta * DEG});
        dev[k] = std::abs(eff(a.reflected, n) - eff(b.reflected, n));
      }
      EXPECT_LT(dev[1], 0.25 * dev[0]) << "n=" << nBar << " delta=" << delta;
      EXPECT_LT(dev[1], 1e-4) << "n=" << nBar << " delta=" << delta;
    }
  }
}
