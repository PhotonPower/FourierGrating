// M20 building blocks: Bessel functions, spectral grid, Hankel-Galerkin matrices, factorization in
// r, propagation matrix, mode basis and flux (Plan sections sec:cylgrid, sec:cylgal, sec:cylfact).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/core/quadrature.hpp"
#include "fg/cyl/bessel.hpp"
#include "fg/cyl/excitation.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;
namespace cyl = fg::cyl;

namespace {

double relErr(const MatC& a, const MatC& b) {
  return (a - b).cwiseAbs().maxCoeff() / b.cwiseAbs().maxCoeff();
}

// int_a^b f(r) dr by composite Gauss-Legendre with many panels (reference quadrature)
template <class F>
cplx integrate(F f, double a, double b, int panels = 60) {
  const auto& gl = fg::quadrature::gaussLegendre(20);
  cplx s = 0.0;
  const double h = (b - a) / panels;
  for (int p = 0; p < panels; ++p)
    for (std::size_t q = 0; q < gl.x.size(); ++q)
      s += 0.5 * h * gl.w[q] * f(a + h * (p + 0.5 * (gl.x[q] + 1)));
  return s;
}

cyl::SpectralGrid realGrid(double kmax, double Rs, std::vector<double> branch = {}, int nodes = 8) {
  cyl::SpectrumOptions o;
  o.nodesPerPanel = nodes;
  o.kmax = kmax;
  o.Rs = Rs;
  o.branchPoints = std::move(branch);
  return cyl::spectralGrid(o);
}

cyl::SpectralGrid deformedGrid(int nodes = 8) {
  cyl::SpectrumOptions o;
  o.nodesPerPanel = nodes;
  o.kmax = 9.0;
  o.Rs = 2.0;
  o.branchPoints = {1.0, 1.5};
  o.deformA = 1.5;
  o.deformB = 2.4;
  o.delta = 0.05 * (o.deformB - o.deformA);
  return cyl::spectralGrid(o);
}

}  // namespace

TEST(CylBessel, RealArgumentMatchesStd) {
  double worst = 0;
  for (int n = 0; n <= 40; ++n)
    for (double x : {0.01, 0.3, 1.0, 4.7, 13.0, 39.5, 41.0, 120.3}) {
      const double ref = std::cyl_bessel_j(double(n), x);
      const cplx v = cyl::besselJ(n, cplx(x, 0.0));
      worst = std::max(worst, std::abs(v - ref) / std::max(std::abs(ref), 1e-300));
      EXPECT_LT(std::abs(v.imag()), 1e-300);
    }
  // relative accuracy, including the deeply decaying region n >> x (large x: scipy table below;
  // std::cyl_bessel_j of libstdc++ is only ~1e-11 there)
  EXPECT_LT(worst, 1e-12);
}

TEST(CylBessel, ComplexArgumentMatchesScipy) {
  struct Ref {
    int n;
    cplx z, J;
  };
  // scipy.special.jv 1.17.1
  const Ref refs[] = {
      {0, {3, 0.2}, {-0.2675350164376209, -0.06804814469901127}},
      {1, {3, 0.2}, {0.34260740511462073, -0.075047162575212836}},
      {5, {12.5, -0.7}, {0.046464737255462116, -0.15875284683099644}},
      {-3, {0.4, 0.1}, {-0.0010780635425370628, -0.00096463544432910436}},
      {2, {50.3, 1.5}, {-0.20221427761761918, -0.15363370624577349}},
      {17, {40, -2}, {-0.36753764630519453, -0.18583479788013191}},
      {40, {35, 0.5}, {0.014424762544017825, 0.0043963968786753498}},
      {0, {450.7, 0.8}, {-0.039487406202048693, 0.020672626034987693}},
      {30, {451.2, -1.1}, {-0.035071588741836621, 0.041467911795106019}},
      {3, {1e-3, 1e-3}, {-4.1666661458333166e-11, 4.1666671874999821e-11}},
      {60, {10, 0.3}, {-1.4379696252895342e-41, 6.9539971443014824e-41}},
      {19, {499.9, 0}, {-0.0013478578717085498342, 0}},  // mpmath (scipy is ~1e-15 off here)
      {1, {499.9, 0}, {0.013827097596907280418, 0}},     // mpmath
      {40, {499.9, 0}, {-0.012875507655243129343, 0}},   // mpmath
  };
  // relative to max(|J|, sqrt(2/(pi|z|))): Miller's recurrence loses ~|z| eps of the local
  // amplitude
  for (const auto& r : refs) {
    const double scale =
        std::max(std::abs(r.J), std::abs(r.z) > 1 ? std::sqrt(2 / (fg::PI * std::abs(r.z))) : 0.0);
    EXPECT_LT(std::abs(cyl::besselJ(r.n, r.z) - r.J) / scale, 1e-13) << r.n << " " << r.z;
  }
  // J_n' from the recursion
  const cplx z(7.3, 0.4);
  const double h = 1e-5;
  const cplx fd = (cyl::besselJ(2, z + h) - cyl::besselJ(2, z - h)) / (2 * h);
  EXPECT_LT(std::abs(cyl::besselJp(2, z) - fd), 1e-9);
}

TEST(CylSpectrum, BranchPointSubstitutionAndPath) {
  // int_0^kh sqrt(kh^2 - k^2) dk = pi kh^2 / 4: spectral thanks to the substitution (16 nodes per
  // panel; without it the square-root end point limits the error to ~1e-5)
  const double kh = 1.5;
  const auto g = realGrid(6.0, 2.0, {1.0, kh}, 16);
  cplx s = 0.0;
  for (Eigen::Index j = 0; j < g.k.size(); ++j)
    if (g.k(j).real() < kh) s += g.w(j) * std::sqrt(kh * kh - g.k(j) * g.k(j));
  EXPECT_NEAR(s.real(), fg::PI * kh * kh / 4, 1e-13);
  // panel edges: 0, branch points, k_max; widths <= pi / (2 R_s)
  EXPECT_EQ(g.edges, (std::vector<double>{0.0, 1.0, kh, 6.0}));
  const auto g8 = realGrid(6.0, 2.0, {1.0, kh});
  EXPECT_EQ(g8.k.size(),
            8 * (2 + 2 + 6));  // [0,1]: 2 panels, [1,1.5]: 2 (two branch ends), [1.5,6]: 6
  EXPECT_GT((g.k.real().tail(g.k.size() - 1) - g.k.real().head(g.k.size() - 1)).minCoeff(), 0.0);
  // analytic integrand: deformed path = real path (Cauchy)
  const auto d = deformedGrid(16);
  auto f = [](cplx k) { return k * std::exp(-k * k / 8.0) * std::cos(k); };
  cplx sr = 0.0, sd = 0.0;
  for (Eigen::Index j = 0; j < d.k.size(); ++j) sd += d.w(j) * f(d.k(j));
  const auto r = realGrid(9.0, 2.0, {1.0, 1.5}, 16);
  for (Eigen::Index j = 0; j < r.k.size(); ++j) sr += r.w(j) * f(r.k(j));
  EXPECT_LT(std::abs(sd - sr), 1e-13);
  EXPECT_LT(d.k.imag().minCoeff(), -0.04);  // the path dips into the lower half plane
}

TEST(CylGalerkin, HankelPairGaussian) {
  // f = r^nu exp(-r^2) <-> f^(k) = k^nu / 2^(nu+1) exp(-k^2/4) (Plan eq:hankel); discrete inverse
  // 1e-12
  const auto g = realGrid(14.0, 4.0);
  for (int nu : {0, 1, 3}) {
    double worst = 0;
    for (double r = 0.0; r <= 3.0; r += 0.137) {
      cplx s = 0.0;
      for (Eigen::Index j = 0; j < g.k.size(); ++j) {
        const cplx k = g.k(j);
        s += g.w(j) * k * std::pow(k, nu) / std::pow(2.0, nu + 1) * std::exp(-k * k / 4.0) *
             cyl::besselJ(nu, k * r);
      }
      worst = std::max(worst, std::abs(s - std::pow(r, nu) * std::exp(-r * r)));
    }
    EXPECT_LT(worst, 1e-12) << nu;
  }
}

TEST(CylGalerkin, LommelEqualsQuadrature) {
  // same-order matrices in closed form (Plan eq:lommel) against Gauss quadrature, real and complex
  // k
  const cyl::RadialSteps steps{{0.4, 1.1, 2.0}, {cplx(3.0, 0.1), 0.5, cplx(-2.0, 0.3)}};
  auto step = [&](double r) {
    for (std::size_t i = 0; i < steps.radii.size(); ++i)
      if (r < steps.radii[i]) return steps.value[i];
    return cplx(0.0);
  };
  for (const auto& grid : {realGrid(6.0, 2.0, {1.0, 1.5}), deformedGrid()})
    for (auto [nu, nup] : {std::pair{0, 0}, {1, 1}, {3, 3}, {-2, -2}, {1, -1}, {-3, 3}}) {
      const MatC L = cyl::galerkin(grid, nu, nup, steps);
      const MatC Q = cyl::galerkinQuadrature(grid, nu, nup, steps.radii, step);
      // closed form for |a - b| r >= 1e-2 cancels to ~eps / 1e-2 relative to the pair
      EXPECT_LT(relErr(L, Q), 5e-14) << nu << " " << nup;
    }
}

TEST(CylGalerkin, ActsAsMultiplication) {
  // [[g]]_{nu<-nu'} f^ = (g f)^ for band-limited f (background identity of Plan eq:cyl_gal):
  // f = r^nu' exp(-r^2) with its exact spectrum, g = 1 + steps; reference by quadrature in r.
  const cyl::RadialSteps steps{{0.7, 1.6}, {cplx(2.0, 0.3), -0.4}};
  const cplx gInf = 1.0;
  const auto grid = realGrid(14.0, 4.0);
  for (auto [nu, nup] : {std::pair{2, 2}, {2, 0}, {1, -1}, {0, 1}}) {
    const int a = std::abs(nup);
    const double sg =
        (nup < 0 && a % 2) ? -1.0 : 1.0;  // spectrum of f = r^a exp(-r^2) in order nup
    VecC fh(grid.k.size());
    for (Eigen::Index j = 0; j < grid.k.size(); ++j)
      fh(j) = sg * std::pow(grid.k(j), a) / std::pow(2.0, a + 1) *
              std::exp(-grid.k(j) * grid.k(j) / 4.0);
    // background only for equal orders (mixed orders: compact functions only, Plan sec:cylgal)
    VecC out = cyl::galerkin(grid, nu, nup, steps) * fh;
    if (nu == nup) out += gInf * fh;
    double worst = 0;
    for (Eigen::Index j = 0; j < grid.k.size(); ++j) {
      const cplx k = grid.k(j);
      auto integrand = [&](double r) {
        const cplx gg = r < 0.7 ? steps.value[0] : (r < 1.6 ? steps.value[1] : cplx(0.0));
        return gg * std::pow(r, a) * std::exp(-r * r) * cyl::besselJ(nu, k * r) * r;
      };
      cplx ref = integrate(integrand, 0, 0.7) + integrate(integrand, 0.7, 1.6);
      if (nu == nup) ref += gInf * fh(j);
      worst = std::max(worst, std::abs(out(j) - ref));
    }
    EXPECT_LT(worst, 1e-12) << nu << " " << nup;
  }
}

TEST(CylFactorization, LiRWithConstantEpsIsIdentity) {
  // eps == eps_inf: frak A = 1, frak B = eps_inf 1 (Plan eq:cyl_li), so Q = eps_inf 1, Q_+- = 0
  const auto grid = realGrid(6.0, 2.0, {1.0});
  const cplx e(2.25, 0.1);
  const cyl::Rings rings{{0.5, 2.0}, {e, e}, e};
  for (int m : {0, 1, 2}) {
    const auto Q = cyl::cylQ(grid, m, rings, cyl::CylFactorization::LiR);
    const MatC Id = MatC::Identity(grid.k.size(), grid.k.size());
    EXPECT_LT((Q.Qpp - e * Id).cwiseAbs().maxCoeff(), 1e-15);
    EXPECT_LT((Q.Qmm - e * Id).cwiseAbs().maxCoeff(), 1e-15);
    EXPECT_LT(Q.Qpm.cwiseAbs().maxCoeff(), 1e-15);
    EXPECT_LT((Q.Qyy - e * Id).cwiseAbs().maxCoeff(), 1e-15);
  }
}

TEST(CylPropagationMatrix, HomogeneousMediumEigenvectors) {
  // Plan, check point test_cyl_propagation_matrix: M^2 = (k^2 - eps) per node, Psi eigenvectors
  const cplx eps(2.25, 0.05);
  const auto grid = deformedGrid();
  const Eigen::Index n = grid.k.size();
  const MatC Id = MatC::Identity(n, n);
  const cyl::CylQ Q{eps * Id, MatC::Zero(n, n), MatC::Zero(n, n), eps * Id, eps * Id};
  const MatC M = cyl::propagationMatrix(grid.k, Q);
  VecC d(4 * n);
  for (int b = 0; b < 4; ++b) d.segment(b * n, n) = grid.k.cwiseProduct(grid.k).array() - eps;
  EXPECT_LT(((M * M) - MatC(d.asDiagonal())).cwiseAbs().maxCoeff(), 1e-13);
  const auto mb = cyl::homogeneousModes(eps, grid.k);
  EXPECT_LT((M * mb.Psi - mb.Psi * mb.lambda.asDiagonal()).cwiseAbs().maxCoeff(), 1e-12);
  EXPECT_LT((mb.Psi * mb.PsiInv - MatC::Identity(4 * n, 4 * n)).cwiseAbs().maxCoeff(), 1e-13);
  // TE (e_+ = e_-) has E_y = 0; k = 0 down-going: H_+ = i sqrt(eps) E_+ (Plan eq:cyl_Psi)
  const auto m0 = cyl::homogeneousModes(eps, VecC::Constant(1, 0.0));
  EXPECT_LT(std::abs(m0.Psi(2, 0) - fg::I * std::sqrt(eps)), 1e-15);
}

TEST(CylFlux, ClosedFormPartialWave) {
  // Plan eq:cyl_flux_closed against eq:cyl_flux for single partial waves in a lossless medium
  const double eps = 2.25;
  const auto grid = realGrid(4.0, 2.0, {1.5});
  const Eigen::Index n = grid.k.size();
  const auto mb = cyl::homogeneousModes(eps, grid.k);
  const cplx ap(0.7, -0.2), am(-0.3, 0.9);
  int checked = 0;
  for (Eigen::Index j = 0; j < n; ++j) {
    const double k = grid.k(j).real();
    if (k >= std::sqrt(eps)) continue;
    for (int dir : {0, 1}) {  // 0: down (a), 1: up (b)
      VecC V = VecC::Zero(4 * n);
      V(2 * n * dir + j) = ap;
      V(2 * n * dir + n + j) = am;
      const VecC F = mb.Psi * V;
      const double b = mb.beta(j).real(), wk = grid.w(j).real() * k;
      const double closed = (dir ? 1 : -1) * 0.5 * fg::PI * wk / b *
                            ((eps - 0.5 * k * k) * (std::norm(ap) + std::norm(am)) -
                             k * k * std::real(ap * std::conj(am)));
      EXPECT_NEAR(cyl::flux(grid, F), closed, 1e-15 * std::abs(closed) + 1e-18);
      ++checked;
    }
  }
  EXPECT_GT(checked, 10);
}

TEST(CylExcitation, BesselProjectionEqualsQuadrature) {
  // column of Plan eq:cyl_gal at an extra node k' without weight (Lommel; mixed orders by quadrature)
  // against quadrature; k'
  // = 0, generic, and on a node (near-coincident pairs by quadrature)
  const cyl::RadialSteps steps{{0.4, 1.1, 2.0}, {cplx(3.0, 0.1), 0.5, cplx(-2.0, 0.3)}};
  auto step = [&](double r) {
    for (std::size_t i = 0; i < steps.radii.size(); ++i)
      if (r < steps.radii[i]) return steps.value[i];
    return cplx(0.0);
  };
  for (const auto& grid : {realGrid(6.0, 2.0, {1.0, 1.5}), deformedGrid()})
    for (auto [nu, nup] : {std::pair{0, 0}, {1, 1}, {-2, -2}, {3, 3}, {1, -1}, {2, 0}, {-1, 1}, {0, -2}})
      for (const double kp : {0.0, 0.83, grid.k(5).real()}) {
        const VecC P = cyl::besselProjection(grid, nu, nup, steps, kp);
        double worst = 0, scale = 0;
        for (Eigen::Index j = 0; j < grid.k.size(); ++j) {
          cplx ref = 0.0;
          for (std::size_t i = 0; i < steps.radii.size(); ++i)
            ref += integrate([&](double r) { return step(r) * cyl::besselJ(nu, grid.k(j) * r) * cyl::besselJ(nup, kp * r) * r; },
                             i ? steps.radii[i - 1] : 0.0, steps.radii[i]);
          worst = std::max(worst, std::abs(P(j) - ref));
          scale = std::max(scale, std::abs(ref));
        }
        // k' = 0, nu' != 0: J_nu'(0) = 0, the column vanishes exactly
        EXPECT_LT(worst / std::max(scale, 1.0), 5e-14) << nu << " " << nup << " " << kp;
      }
}

TEST(CylExcitation, JacobiAngerSource) {
  // Plan eq:cyl_jacobi: the spectra of P_+- = (eps - eps_inf) E_+-,m of a plane wave, with the
  // harmonic E_+,m(r) = 1/(2 pi) int e^{-i phi} (E_x + i E_z) e^{-i m phi} dphi taken numerically
  const cyl::Rings rings{{0.6, 1.3}, {cplx(2.25, 0.1), 3.0}, 1.2};
  const auto grid = realGrid(5.0, 1.3, {1.0});
  const double kpar = 0.7, Psi = 0.4;
  for (const int m : {-2, 0, 1, 3}) {
    const cyl::CylQ Q = cyl::cylQ(grid, m, rings, cyl::CylFactorization::Laurent);
    const cyl::PlaneWaveSource B = cyl::planeWaveSource(grid, m, rings, kpar, Psi, cyl::CylFactorization::Laurent, Q);
    const Eigen::Index n = grid.k.size();
    auto harmonic = [&](int shift, double r) {  // E_x + i E_z = e^{i k r cos(phi - Psi)}; shift = +-1 (E_pm) or 0 (E_y)
      const int K = 64;
      cplx s = 0.0;
      for (int q = 0; q < K; ++q) {
        const double phi = 2 * fg::PI * q / K;
        s += std::exp(cplx(0, kpar * r * std::cos(phi - Psi) - (m + shift) * phi));
      }
      return s / double(K);
    };
    double worst = 0;
    for (Eigen::Index j = 0; j < grid.k.size(); j += 7)
      for (const int shift : {1, -1, 0}) {
        cplx ref = 0.0;
        for (std::size_t i = 0; i < rings.radii.size(); ++i)
          ref += (rings.eps[i] - rings.epsInf) *
                 integrate([&](double r) { return harmonic(shift, r) * cyl::besselJ(m + shift, grid.k(j) * r) * r; },
                           i ? rings.radii[i - 1] : 0.0, rings.radii[i]);
        const cplx got = shift == 1 ? B.Bpm(j, 0) : shift == -1 ? B.Bpm(n + j, 1) : B.By(j);
        worst = std::max(worst, std::abs(got - ref));
      }
    EXPECT_LT(worst, 1e-13) << m;
  }
  EXPECT_EQ(cyl::harmonicCutoff(0.0), 2);
  EXPECT_EQ(cyl::harmonicCutoff(8.0), 18);
}

TEST(CylExcitation, LiConsistentSourceOnNode) {
  // Li in r: P = b^ - Q a^ is (Q - eps_inf) E_bg (Plan eq:cyl_scat, consistency). With kpar on a
  // node l, E_bg = e J(k_l r) lies on the grid (spectrum e / (w_l k_l) at node l), so the identity
  // can be checked against the matrices.
  const cyl::Rings rings{{0.6, 1.3}, {cplx(2.25, 0.1), 3.0}, 1.2};
  const auto grid = realGrid(5.0, 1.3, {1.0});
  const Eigen::Index n = grid.k.size(), l = 9;
  const double kpar = grid.k(l).real(), Psi = 0.3;
  for (const int m : {0, 2, -1}) {
    const cyl::CylQ Q = cyl::cylQ(grid, m, rings, cyl::CylFactorization::LiR);
    const cyl::PlaneWaveSource B = cyl::planeWaveSource(grid, m, rings, kpar, Psi, cyl::CylFactorization::LiR, Q);
    MatC Qf(2 * n, 2 * n);
    Qf << Q.Qpp, Q.Qpm, Q.Qmp, Q.Qmm;
    const cplx cp(0.4, -1.1), cm(-0.7, 0.2);
    auto phase = [Psi](int nu) { return std::polar(1.0, nu * (0.5 * fg::PI - Psi)); };
    VecC E = VecC::Zero(2 * n);
    E(l) = cp * phase(m + 1) / (grid.w(l) * grid.k(l));
    E(n + l) = cm * phase(m - 1) / (grid.w(l) * grid.k(l));
    const VecC ref = Qf * E - rings.epsInf * E, got = B.Bpm * Eigen::Vector2cd(cp, cm);
    EXPECT_LT((got - ref).cwiseAbs().maxCoeff() / ref.cwiseAbs().maxCoeff(), 1e-12) << m;
  }
}
