// M20 building blocks: Bessel functions, spectral grid, Hankel-Galerkin matrices, factorization in
// r, propagation matrix, mode basis and flux (Plan sections sec:cylgrid, sec:cylgal, sec:cylfact).
#include <gtest/gtest.h>

#include "../support/std_bessel.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <tuple>
#include <cmath>

#include "fg/core/quadrature.hpp"
#include "fg/cyl/bessel.hpp"
#include "fg/cyl/coupled.hpp"
#include "fg/cyl/excitation.hpp"
#include "fg/cyl/field.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/smatrix/smatrix.hpp"
#include "fg/solvers/rcwa.hpp"

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
  FG_SKIP_WITHOUT_STD_BESSEL();
  double worst = 0;
  for (int n = 0; n <= 40; ++n)
    for (double x : {0.01, 0.3, 1.0, 4.7, 13.0, 39.5, 41.0, 120.3}) {
      const double ref = fgtest::cyl_bessel_j(double(n), x);
      const cplx v = cyl::besselJ(n, cplx(x, 0.0));
      worst = std::max(worst, std::abs(v - ref) / std::max(std::abs(ref), 1e-300));
      EXPECT_LT(std::abs(v.imag()), 1e-300);
    }
  // relative accuracy, including the deeply decaying region n >> x (large x: scipy table below;
  // fgtest::cyl_bessel_j of libstdc++ is only ~1e-11 there)
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

TEST(CylSpectrum, GradedPanelsAroundAResonance) {
  // a pole p just above the real section of the path (leaky plasmon between the branch points, V21h):
  // int_0^kmax dk / (k - p) = Log(kmax - p) - Log(-p); uniform panels of width pi / (2 R_s) miss the
  // Lorentzian of width Im p, graded panels Re p +- Im p 2^j resolve it
  const cplx p(1.2, 2e-3);
  const double kmax = 6.0;
  const cplx exact = std::log(kmax - p) - std::log(-p);
  auto integral = [&](const cyl::SpectralGrid& g) {
    cplx s = 0.0;
    for (Eigen::Index j = 0; j < g.k.size(); ++j) s += g.w(j) / (g.k(j) - p);
    return s;
  };
  cyl::SpectrumOptions o;
  o.nodesPerPanel = 16;
  o.kmax = kmax;
  o.Rs = 2.0;
  o.branchPoints = {1.0, 1.5};
  const auto plain = cyl::spectralGrid(o);
  o.resonances = {p};
  const auto graded = cyl::spectralGrid(o);
  EXPECT_GT(std::abs(integral(plain) - exact), 1e-3);
  EXPECT_LT(std::abs(integral(graded) - exact), 1e-10);
  // the grading stays between the branch points (edges kept) and adds few panels
  for (double b : {0.0, 1.0, 1.5, kmax})
    EXPECT_NE(std::find(graded.edges.begin(), graded.edges.end(), b), graded.edges.end());
  EXPECT_TRUE(std::is_sorted(graded.edges.begin(), graded.edges.end()));
  EXPECT_LT(graded.k.size(), plain.k.size() + 16 * 20);
  const auto inner = std::count_if(graded.edges.begin(), graded.edges.end(),
                                   [](double e) { return e > 1.0 && e < 1.5; });
  EXPECT_GE(inner, 8);
  // a resonance on the deformed section is left to the complex path; one with Im p = 0 gets a floor width
  o.deformA = 1.5;
  o.deformB = 2.4;
  o.delta = 0.05;
  o.resonances = {cplx(2.0, 1e-3)};
  auto d = o;
  d.resonances.clear();
  EXPECT_EQ(cyl::spectralGrid(o).edges, cyl::spectralGrid(d).edges);
  o.resonances = {cplx(1.2, 0.0)};
  EXPECT_LT(cyl::spectralGrid(o).k.size(), plain.k.size() + 16 * 40);
}

TEST(CylSpectrum, PlanarBoundModesOfAFilm) {
  // film n = 2.1, 0.3 um on glass under air at 633 nm (the V21d contact hole background): TE0, TM0, TE1 from
  // the dispersion relation k_f h = atan(p_c / k_f) + atan(p_s / k_f) + m pi (p = gamma, TM: gamma eps_f / eps)
  const double k0 = 2 * fg::PI / 0.6328, h = 0.3 * k0, ef = 4.41, es = 2.25, ec = 1.0;
  auto dispersion = [&](double b, bool tm, int m) {
    const double kf = std::sqrt(ef - b * b), gs = std::sqrt(b * b - es), gc = std::sqrt(b * b - ec);
    const double ps = tm ? gs * ef / es : gs, pc = tm ? gc * ef / ec : gc;
    return kf * h - std::atan(pc / kf) - std::atan(ps / kf) - m * fg::PI;
  };
  auto root = [&](bool tm, int m) {
    double a = 1.5 + 1e-12, b = 2.1 - 1e-12;
    for (int it = 0; it < 200; ++it) {
      const double c = 0.5 * (a + b);
      (dispersion(a, tm, m) * dispersion(c, tm, m) <= 0 ? b : a) = c;
    }
    return 0.5 * (a + b);
  };
  const double kb = 1.2 * 2.1;
  const auto modes = cyl::planarBoundModes(es, ec, {ef}, {h}, 1.5, kb);
  ASSERT_EQ(modes.size(), 3u);
  EXPECT_NEAR(modes[0], root(false, 1), 1e-12);  // TE1 near cut-off
  EXPECT_NEAR(modes[1], root(true, 0), 1e-12);
  EXPECT_NEAR(modes[2], root(false, 0), 1e-12);
  EXPECT_NEAR(modes[0], 1.5666, 1e-4);
  // the same film split into two layers; superstrate and substrate exchanged (same modes)
  EXPECT_EQ(cyl::planarBoundModes(es, ec, {ef, ef}, {0.4 * h, 0.6 * h}, 1.5, kb).size(), 3u);
  const auto flipped = cyl::planarBoundModes(ec, es, {ef}, {h}, 1.5, kb);
  ASSERT_EQ(flipped.size(), 3u);
  for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(flipped[i], modes[i], 1e-12);
  // no guiding layer, a lossy layer: none
  EXPECT_TRUE(cyl::planarBoundModes(es, ec, {1.45 * 1.45}, {h}, 1.5, kb).empty());
  EXPECT_TRUE(cyl::planarBoundModes(es, ec, {cplx(ef, 1e-3)}, {h}, 1.5, kb).empty());
}

TEST(CylSpectrum, GradedPanelsAroundBoundModes) {
  // a real pole p on the deformed section, close to k_a (guided mode near cut-off): the path passes below,
  // int_path dk / (k - p) = ln((k_max - p) / p) + i pi. The sine path is only delta sin(pi (p - k_a) / L)
  // away from p; graded panels with that width resolve it (V21d, Plan v2.38)
  const double p = 1.5666, kmax = 8.4;
  const cplx exact = std::log((kmax - p) / p) + cplx(0, fg::PI);
  auto integral = [&](const cyl::SpectralGrid& g) {
    cplx s = 0.0;
    for (Eigen::Index j = 0; j < g.k.size(); ++j) s += g.w(j) / (g.k(j) - p);
    return s;
  };
  cyl::SpectrumOptions o;
  o.kmax = kmax;
  o.Rs = 2.48;  // r = 0.25 um at 633 nm
  o.branchPoints = {1.5, 1.0};
  o.deformA = 1.5;
  o.deformB = 2.52;
  o.delta = 0.05 * (o.deformB - o.deformA);
  const auto plain = cyl::spectralGrid(o);
  o.boundModes = {p, 1.9026, 1.9655};
  const auto graded = cyl::spectralGrid(o);
  std::printf("pole near k_a: plain %.1e (%ld nodes), graded %.1e (%ld nodes)\n", std::abs(integral(plain) - exact),
              static_cast<long>(plain.k.size()), std::abs(integral(graded) - exact), static_cast<long>(graded.k.size()));
  EXPECT_GT(std::abs(integral(plain) - exact), 1.0);  // measured 1.8
  EXPECT_LT(std::abs(integral(graded) - exact), 1e-6);  // measured 2.4e-7 (8 nodes per panel)
  EXPECT_TRUE(std::is_sorted(graded.edges.begin(), graded.edges.end()));
  for (double b : {0.0, 1.0, 1.5, 2.52, kmax})
    EXPECT_NE(std::find(graded.edges.begin(), graded.edges.end(), b), graded.edges.end());
  // bound modes are ignored on the real path and outside the deformed section
  auto r = o;
  r.delta = 0;
  auto r0 = r;
  r0.boundModes.clear();
  EXPECT_EQ(cyl::spectralGrid(r).edges, cyl::spectralGrid(r0).edges);
  o.boundModes = {1.2, 3.0};
  EXPECT_EQ(cyl::spectralGrid(o).edges, plain.edges);
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
  // relative to max|k^2 - eps|: the rounding of M*M depends on FMA contraction (arm64 clang contracts by default)
  const double scale = std::max(1.0, d.cwiseAbs().maxCoeff());
  EXPECT_LT(((M * M) - MatC(d.asDiagonal())).cwiseAbs().maxCoeff(), 1e-13 * scale);
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

// Near-field building blocks (Plan "M20: S-Matrix, Loeser und Nahfeld"): synthesis, radial rule,
// frak A for G_pm and the background column a^.
TEST(CylField, SynthesisOfTheGaussianHankelPair) {
  // order nu: r^nu e^{-r^2} <-> k^nu / 2^{nu+1} e^{-k^2/4} (Plan eq:hankel); harmonic m = nu + 1 uses
  // T_m for order nu, m = nu - 1 uses T_p, m = nu uses T_0; negative orders through J_{-n} = (-1)^n J_n.
  // 12 nodes per panel: the Nystroem sum is then exact to rounding up to r = 1.2 R_s (8 nodes, the
  // default: 1e-12 for r < R_s / 2, 1e-10 at r = 1.2 R_s)
  const auto grid = realGrid(14.0, 3.0, {1.0, 1.5}, 12);
  const std::vector<double> r{0.0, 0.3, 1.1, 2.4, 3.7};
  for (const int nu : {0, 1, 2}) {
    VecC f(grid.k.size());
    for (Eigen::Index j = 0; j < f.size(); ++j)
      f(j) = std::pow(grid.k(j), nu) / std::pow(2.0, nu + 1) * std::exp(-grid.k(j) * grid.k(j) / 4.0);
    for (const int m : {nu + 1, nu - 1, nu, -nu}) {
      const cyl::Synthesis S = cyl::synthesis(grid, m, r);
      const MatC& T = m == nu + 1 ? S.Tm : m == nu - 1 ? S.Tp : S.T0;
      const double sign = (m == -nu && nu % 2 == 1) ? -1.0 : 1.0;  // J_{-nu} = (-1)^nu J_nu
      const VecC got = T * f;
      for (std::size_t s = 0; s < r.size(); ++s)
        EXPECT_NEAR(std::abs(got(static_cast<Eigen::Index>(s)) - sign * std::pow(r[s], nu) * std::exp(-r[s] * r[s])), 0.0,
                    1e-12)
            << nu << " " << m << " " << r[s];
    }
  }
}

TEST(CylField, RadialRuleIntegratesBesselProducts) {
  const std::vector<double> radii{0.4, 1.0, 2.5};
  const cyl::RadialRule q = cyl::radialRule(radii, 9.0);
  ASSERT_EQ(q.r.size(), q.w.size());
  double area = 0;
  for (std::size_t s = 0; s < q.r.size(); ++s) {
    area += q.w[s];
    EXPECT_TRUE(q.r[s] < radii[q.ring[s]] && q.r[s] > (q.ring[s] ? radii[q.ring[s] - 1] : 0.0));
  }
  EXPECT_NEAR(area, 0.5 * 2.5 * 2.5, 1e-14);
  // the Lommel integral of the galerkin matrices on the rule (J_3(7.9 r) J_3(8.6 r) r over [0, 2.5])
  cplx s = 0.0;
  for (std::size_t i = 0; i < q.r.size(); ++i) s += q.w[i] * cyl::besselJ(3, 7.9 * q.r[i]) * cyl::besselJ(3, 8.6 * q.r[i]);
  const cplx ref = integrate([](double r) { return cyl::besselJ(3, 7.9 * r) * cyl::besselJ(3, 8.6 * r) * r; }, 0.0, 2.5);
  EXPECT_LT(std::abs(s - ref), 1e-14);
}

TEST(CylField, FrakAAndBackgroundColumnOnNode) {
  // Li in r: A E_bg = E_bg + Apm (c_+; c_-) (Plan eq:cyl_li, the delta parts of a plane wave); with kpar
  // on a node E_bg lies on the grid, as in LiConsistentSourceOnNode
  const cyl::Rings rings{{0.6, 1.3}, {cplx(2.25, 0.1), 3.0}, 1.2};
  const auto grid = realGrid(5.0, 1.3, {1.0});
  const Eigen::Index n = grid.k.size(), l = 9;
  const double kpar = grid.k(l).real(), Psi = 0.3;
  for (const int m : {0, 2, -1}) {
    const cyl::CylQ Q = cyl::cylQ(grid, m, rings, cyl::CylFactorization::LiR);
    ASSERT_EQ(Q.A.rows(), 2 * n);
    const cyl::PlaneWaveSource B = cyl::planeWaveSource(grid, m, rings, kpar, Psi, cyl::CylFactorization::LiR, Q);
    ASSERT_EQ(B.Apm.rows(), 2 * n);
    const cplx cp(0.4, -1.1), cm(-0.7, 0.2);
    auto phase = [Psi](int nu) { return std::polar(1.0, nu * (0.5 * fg::PI - Psi)); };
    VecC E = VecC::Zero(2 * n);
    E(l) = cp * phase(m + 1) / (grid.w(l) * grid.k(l));
    E(n + l) = cm * phase(m - 1) / (grid.w(l) * grid.k(l));
    const VecC ref = Q.A * E - E, got = B.Apm * Eigen::Vector2cd(cp, cm);
    EXPECT_LT((got - ref).cwiseAbs().maxCoeff() / ref.cwiseAbs().maxCoeff(), 1e-12) << m;
    // Q = B A^{-1}: Q A = B, so (Q A - eps_inf) E_bg = b^ (the other half of the consistent source)
    MatC Qf(2 * n, 2 * n);
    Qf << Q.Qpp, Q.Qpm, Q.Qmp, Q.Qmm;
    const VecC bh = B.Bpm * Eigen::Vector2cd(cp, cm) + Qf * got;
    EXPECT_LT((Qf * (Q.A * E) - rings.epsInf * E - bh).cwiseAbs().maxCoeff() / bh.cwiseAbs().maxCoeff(), 1e-10) << m;
  }
  EXPECT_EQ(cyl::cylQ(grid, 1, rings, cyl::CylFactorization::Laurent).A.size(), 0);
}

TEST(CylField, FieldsFromGInvertTheAuxiliaryFields) {
  // G_pm = D_r / eps_inf +- i E_phi (Plan eq:cyl_G) from E_r, E_phi and back
  const cplx eps(-9.0, 1.3), ei = 1.7, Er(0.3, -0.8), Ephi(1.1, 0.4);
  const cplx gp = eps * Er / ei + fg::I * Ephi, gm = eps * Er / ei - fg::I * Ephi;
  cplx ep, em;
  cyl::fieldsFromG(eps, ei, gp, gm, ep, em);
  EXPECT_LT(std::abs(ep - (Er + fg::I * Ephi)), 1e-15);
  EXPECT_LT(std::abs(em - (Er - fg::I * Ephi)), 1e-15);
}

// ---- reverse-mode rules of the cylindrical building blocks (gradients V21i): central difference quotients

namespace {

double reDot(const MatC& Xbar, const MatC& dX) { return (Xbar.conjugate().array() * dX.array()).sum().real(); }

MatC randomMat(Eigen::Index r, Eigen::Index c, unsigned seed) {
  std::srand(seed);
  return MatC::Random(r, c);
}

}  // namespace

TEST(CylAdjoint, FluxBar) {
  const auto g = deformedGrid();
  const Eigen::Index n = g.k.size();
  const VecC F = randomMat(4 * n, 1, 3);
  const VecC Fb = cyl::fluxBar(g, F, 0.7);
  const double h = 1e-6;
  for (const Eigen::Index i : {Eigen::Index(1), n + 2, 2 * n + 1, 3 * n + 2, 3 * n + n - 1})
    for (const cplx dir : {cplx(1, 0), cplx(0, 1)}) {
      VecC a = F, b = F;
      a(i) += h * dir;
      b(i) -= h * dir;
      const double fd = 0.7 * (cyl::flux(g, a) - cyl::flux(g, b)) / (2 * h);
      EXPECT_NEAR(fd, std::real(std::conj(Fb(i)) * dir), 1e-8 * (1 + std::abs(fd))) << i;
    }
}

// Half-space materials (Plan v2.48): eps of homogeneousModes through Psi and PsiInv, lossy and lossless media, nodes
// on the deformed path
TEST(CylAdjoint, HomogeneousModesBar) {
  const auto g = deformedGrid();
  const Eigen::Index n = g.k.size();
  const MatC Pb = randomMat(4 * n, 4 * n, 5), Qb = randomMat(4 * n, 4 * n, 6);
  // lossless (superstrate: real direction only; Im would cross the branch cut of beta) off the branch points of the
  // grid, and lossy
  for (const cplx eps : {cplx(2.89, 0.0), cplx(-4.0, 0.6)}) {
    auto L = [&](cplx e) {
      const fg::ModeBasis mb = cyl::homogeneousModes(e, g.k);
      return (Pb.conjugate().cwiseProduct(mb.Psi).sum() + Qb.conjugate().cwiseProduct(mb.PsiInv).sum()).real();
    };
    const cplx bar = cyl::homogeneousModesBar(eps, g.k, Pb, Qb);
    const double h = 1e-6;
    for (const cplx dir : {cplx(1, 0), cplx(0, 1)}) {
      if (eps.imag() == 0.0 && dir.imag() != 0.0) continue;
      const double fd = (L(eps + h * dir) - L(eps - h * dir)) / (2 * h);
      EXPECT_NEAR(std::real(std::conj(bar) * dir), fd, 1e-7 * (1 + std::abs(fd))) << eps << " " << dir;
    }
  }
}

TEST(CylAdjoint, PropagationMatrixAndFactorization) {
  cyl::SpectrumOptions o;
  o.kmax = 4.0;
  o.Rs = 2.0;
  o.branchPoints = {1.0, 1.5};
  o.deformA = 1.5;
  o.deformB = 2.4;
  o.delta = 0.05;
  const auto g = cyl::spectralGrid(o);
  const Eigen::Index n = g.k.size();
  const cyl::Rings rings{{0.7, 1.6}, {cplx(4.0, 0.3), cplx(-6.0, 1.0)}, cplx(2.25, 0.0)};
  const double h = 1e-6;
  for (const auto f : {cyl::CylFactorization::Laurent, cyl::CylFactorization::LiR})
    for (const int m : {0, 1, -2, 3}) {
      const MatC Mb = randomMat(4 * n, 4 * n, 11 + m);
      auto L = [&](const cyl::Rings& r) {
        return reDot(Mb, cyl::propagationMatrix(g.k, cyl::cylQ(g, m, r, f)));
      };
      const cyl::CylQ Q = cyl::cylQ(g, m, rings, f);
      const cyl::CylQ Qb = cyl::propagationMatrixBar(g.k, Q, Mb);
      const cyl::CylQBar b = cyl::cylQBar(g, m, rings, f, Q, Qb);
      auto check = [&](const std::function<void(cyl::Rings&, double)>& move, double analytic, const char* what) {
        cyl::Rings a = rings, c = rings;
        move(a, h);
        move(c, -h);
        const double fd = (L(a) - L(c)) / (2 * h);
        EXPECT_NEAR(analytic, fd, 1e-6 * (1 + std::abs(fd))) << what << " m = " << m << " f = " << int(f);
      };
      for (std::size_t i = 0; i < 2; ++i) {
        check([i](cyl::Rings& r, double d) { r.eps[i] += d; }, b.eps[i].real(), "Re eps");
        check([i](cyl::Rings& r, double d) { r.eps[i] += cplx(0, d); }, b.eps[i].imag(), "Im eps");
        check([i](cyl::Rings& r, double d) { r.radii[i] += d; }, b.radii[i], "radius");
      }
      check([](cyl::Rings& r, double d) { r.epsInf += d; }, b.epsInf.real(), "Re eps_inf");
      check([](cyl::Rings& r, double d) { r.epsInf += cplx(0, d); }, b.epsInf.imag(), "Im eps_inf");
    }
}

TEST(CylAdjoint, PlaneWaveSource) {
  cyl::SpectrumOptions o;
  o.kmax = 4.0;
  o.Rs = 2.0;
  o.branchPoints = {1.0, 1.5};
  o.deformA = 1.5;
  o.deformB = 2.4;
  o.delta = 0.05;
  const auto g = cyl::spectralGrid(o);
  const Eigen::Index n = g.k.size();
  const cyl::Rings rings{{0.7, 1.6}, {cplx(4.0, 0.3), cplx(-6.0, 1.0)}, cplx(2.25, 0.0)};
  const double kpar = 0.8, Psi = 0.6, h = 1e-6;
  for (const auto f : {cyl::CylFactorization::Laurent, cyl::CylFactorization::LiR})
    for (const int m : {0, 2, -1}) {
      const MatC Bb = randomMat(2 * n, 2, 21 + m);
      const VecC yb = randomMat(n, 1, 31 + m);
      auto L = [&](const cyl::Rings& r) {
        const auto B = cyl::planeWaveSource(g, m, r, kpar, Psi, f, cyl::cylQ(g, m, r, f));
        return reDot(Bb, B.Bpm) + reDot(yb, B.By);
      };
      const cyl::CylQ Q = cyl::cylQ(g, m, rings, f);
      const auto B = cyl::planeWaveSource(g, m, rings, kpar, Psi, f, Q);
      cyl::CylQ Qb{MatC::Zero(n, n), MatC::Zero(n, n), MatC::Zero(n, n), MatC::Zero(n, n), MatC::Zero(n, n), MatC()};
      cyl::CylQBar b;
      cyl::planeWaveSourceBar(g, m, rings, kpar, Psi, f, Q, B, Bb, yb, Qb, b);
      const cyl::CylQBar c = cyl::cylQBar(g, m, rings, f, Q, Qb);  // the path through Q (Li in r)
      for (std::size_t i = 0; i < 2; ++i) {
        b.eps[i] += c.eps[i];
        b.radii[i] += c.radii[i];
      }
      b.epsInf += c.epsInf;
      auto check = [&](const std::function<void(cyl::Rings&, double)>& move, double analytic, const char* what) {
        cyl::Rings a = rings, d = rings;
        move(a, h);
        move(d, -h);
        const double fd = (L(a) - L(d)) / (2 * h);
        EXPECT_NEAR(analytic, fd, 1e-6 * (1 + std::abs(fd))) << what << " m = " << m << " f = " << int(f);
      };
      for (std::size_t i = 0; i < 2; ++i) {
        check([i](cyl::Rings& r, double d) { r.eps[i] += d; }, b.eps[i].real(), "Re eps");
        check([i](cyl::Rings& r, double d) { r.eps[i] += cplx(0, d); }, b.eps[i].imag(), "Im eps");
        check([i](cyl::Rings& r, double d) { r.radii[i] += d; }, b.radii[i], "radius");
      }
      check([](cyl::Rings& r, double d) { r.epsInf += d; }, b.epsInf.real(), "Re eps_inf");
    }
}

// Grouped S-matrices (M21 acceleration): homogeneous layers, the substrate transfer and their star products are
// block-diagonal over the G = n_f spectral nodes of a class (4 x 4 per node); per group they equal the dense
// stableLayerSMatrix, fromTransfer and star
TEST(CylGrouped, EqualTheDenseOperations) {
  VecC k(7);
  for (Eigen::Index j = 0; j < k.size(); ++j) k(j) = cplx(0.25 + 0.4 * j, j == 3 ? -0.05 : 0.0);  // one deformed node
  const VecC kt = cyl::tiledK(k, 3);  // a class of three harmonics
  const Eigen::Index G = kt.size();
  const fg::ModeBasis ref = cyl::homogeneousModes(1.0, kt), sub = cyl::homogeneousModes(2.25, kt);
  const fg::LayerModes film = fg::layerModesHomogeneous(cyl::homogeneousModes(cplx(2.1, 0.3), kt)),
                       glass = fg::layerModesHomogeneous(cyl::homogeneousModes(2.25, kt));
  auto diff = [](const fg::SMatrix& a, const fg::SMatrix& b) {
    return std::max({(a.S11 - b.S11).cwiseAbs().maxCoeff(), (a.S12 - b.S12).cwiseAbs().maxCoeff(),
                     (a.S21 - b.S21).cwiseAbs().maxCoeff(), (a.S22 - b.S22).cwiseAbs().maxCoeff()});
  };
  fg::RCWASliceRecord ra, rb;
  const fg::SMatrix Lf = fg::stableLayerSMatrix(film, 0.7, ref, &ra), Lg = fg::groupedLayerSMatrix(film, 0.7, ref, G, &rb);
  const fg::SMatrix Gf = fg::stableLayerSMatrix(glass, 1.3, ref), Gg = fg::groupedLayerSMatrix(glass, 1.3, ref, G);
  const fg::SMatrix Tf = fg::fromTransfer(ref.PsiInv * sub.Psi), Tg = fg::groupedFromTransfer(ref.PsiInv, sub.Psi, G);
  const double dl = std::max(diff(Lf, Lg), diff(Gf, Gg)), dt = diff(Tf, Tg);
  const double dG = (ra.G - rb.G).cwiseAbs().maxCoeff();
  const double ds = diff(fg::star(fg::star(Tf, Lf), Gf), fg::groupedStar(fg::groupedStar(Tg, Lg, G), Gg, G));
  std::printf("grouped: layer %.1e, record G %.1e, transfer %.1e, star %.1e\n", dl, dG, dt, ds);
  EXPECT_LT(dl, 1e-14);
  EXPECT_LT(dG, 1e-14);
  EXPECT_LT(dt, 1e-14);
  EXPECT_LT(ds, 1e-14);
  // the reference product of a dense slice: grouped Psi^-1 (W) equals the dense product
  // (modes of the substrate with a small dense perturbation: a well-conditioned slice that couples every node)
  MatC W(4 * G, 4 * G);
  for (Eigen::Index i = 0; i < W.rows(); ++i)
    for (Eigen::Index j = 0; j < W.cols(); ++j) W(i, j) = 0.05 * cplx(std::sin(0.3 * i + j), std::cos(i - 0.7 * j));
  W += sub.Psi;
  fg::LayerModes m = fg::layerModesHomogeneous(sub);
  m.Wd = W.leftCols(2 * G);
  m.Wu = W.rightCols(2 * G);
  const double dd = diff(fg::stableLayerSMatrix(m, 0.4, ref), fg::stableLayerSMatrix(m, 0.4, ref, nullptr, G));
  std::printf("grouped reference product of a dense slice: %.1e\n", dd);
  EXPECT_LT(dd, 1e-12);
  // star products with one grouped factor (flag SMatrix::groups) equal the dense ones, also with sources
  const fg::SMatrix D = fg::stableLayerSMatrix(m, 0.4, ref);
  fg::SMatrix Ld = Lg;
  Ld.groups = 0;
  const double dm = std::max(diff(fg::star(Lg, D), fg::star(Ld, D)), diff(fg::star(D, Lg), fg::star(D, Ld)));
  const VecC u = VecC::LinSpaced(2 * G, cplx(0.1, 0.2), cplx(-0.3, 0.4)), v = u.reverse();
  const fg::SourcedSMatrix a{Lg, u, v}, ad{Ld, u, v}, b{D, v, u}, z{Gg, VecC::Zero(2 * G), VecC::Zero(2 * G)};
  double dsrc = 0;
  for (const auto& [x, y, xd, yd] : {std::tuple{a, b, ad, b}, std::tuple{b, a, b, ad}, std::tuple{a, z, ad, z}}) {
    const fg::SourcedSMatrix c = fg::star(x, y), cd = fg::star(xd, yd);
    dsrc = std::max({dsrc, diff(c.S, cd.S), (c.u - cd.u).cwiseAbs().maxCoeff(), (c.v - cd.v).cwiseAbs().maxCoeff()});
  }
  std::printf("star with a grouped factor: %.1e, with sources %.1e\n", dm, dsrc);
  EXPECT_LT(dm, 1e-13);
  EXPECT_LT(dsrc, 1e-13);
}
