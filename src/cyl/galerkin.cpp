#include "fg/cyl/galerkin.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/cyl/bessel.hpp"

namespace fg::cyl {

namespace {

double sgnOrder(int nu) {
  return (nu < 0 && (-nu) % 2 == 1) ? -1.0 : 1.0;
}  // J_{-n} = (-1)^n J_n

double maxAbsK(const SpectralGrid& grid) {
  return grid.k.cwiseAbs().maxCoeff();
}

/// Composite Gauss-Legendre nodes r_s and weights W_s on [0, breaks.back()].
void radialRule(const std::vector<double>& breaks, double kabs, int nodes, std::vector<double>& r,
                std::vector<double>& W) {
  const auto& gl = quadrature::gaussLegendre(nodes);
  const double hmax = PI / (2 * std::max(kabs, 1e-300));
  double a = 0.0;
  for (double b : breaks) {
    if (b <= a) continue;
    const int np = std::max(1, static_cast<int>(std::ceil((b - a) / hmax - 1e-12)));
    const double h = (b - a) / np;
    for (int p = 0; p < np; ++p)
      for (std::size_t q = 0; q < gl.x.size(); ++q) {
        r.push_back(a + h * (p + 0.5 * (gl.x[q] + 1)));
        W.push_back(0.5 * h * gl.w[q]);
      }
    a = b;
  }
}

/// T(j, s) = J_nu(k_j r_s).
MatC besselTable(const VecC& k, const std::vector<double>& r, int nu) {
  const int a = std::abs(nu);
  const double sg = sgnOrder(nu);
  MatC T(k.size(), static_cast<Eigen::Index>(r.size()));
  for (Eigen::Index j = 0; j < k.size(); ++j)
    for (std::size_t s = 0; s < r.size(); ++s)
      T(j, static_cast<Eigen::Index>(s)) =
          sg * besselJAll(a, k(j) * r[s])[static_cast<std::size_t>(a)];
  return T;
}

}  // namespace

MatC galerkinQuadrature(const SpectralGrid& grid, int nu, int nup,
                        const std::vector<double>& breaks, const std::function<cplx(double)>& g) {
  FG_REQUIRE(!breaks.empty(), "galerkinQuadrature: breaks");
  std::vector<double> r, W;
  radialRule(breaks, maxAbsK(grid), 16, r, W);
  const MatC Tn = besselTable(grid.k, r, nu), Tp = besselTable(grid.k, r, nup);
  VecC d(static_cast<Eigen::Index>(r.size()));
  for (std::size_t s = 0; s < r.size(); ++s)
    d(static_cast<Eigen::Index>(s)) = g(r[s]) * r[s] * W[s];
  const VecC wk = grid.w.cwiseProduct(grid.k);
  return (Tn * d.asDiagonal() * Tp.transpose()) * wk.asDiagonal();
}

MatC galerkin(const SpectralGrid& grid, int nu, int nup, const RadialSteps& g) {
  FG_REQUIRE(g.radii.size() == g.value.size() && !g.radii.empty(), "galerkin: radii and values");
  for (std::size_t i = 0; i < g.radii.size(); ++i)
    FG_REQUIRE(g.radii[i] > (i ? g.radii[i - 1] : 0.0), "galerkin: radii increasing and > 0");
  auto step = [&](double r) {
    for (std::size_t i = 0; i < g.radii.size(); ++i)
      if (r < g.radii[i]) return g.value[i];
    return cplx(0.0);
  };
  if (std::abs(nu) != std::abs(nup)) return galerkinQuadrature(grid, nu, nup, g.radii, step);

  // Same order up to sign: Lommel integrals (Plan eq:lommel) of the antiderivative F_jl(r).
  const int n = std::abs(nu);
  const double sg = sgnOrder(nu) * sgnOrder(nup);
  const Eigen::Index nk = grid.k.size();
  const std::size_t nr = g.radii.size();
  MatC J(nk, static_cast<Eigen::Index>(nr)), Jp(nk, static_cast<Eigen::Index>(nr));
  for (Eigen::Index j = 0; j < nk; ++j)
    for (std::size_t i = 0; i < nr; ++i) {
      const auto v = besselJAll(n + 1, grid.k(j) * g.radii[i]);
      const cplx Jm1 = n == 0 ? -v[1] : v[static_cast<std::size_t>(n - 1)];
      J(j, static_cast<Eigen::Index>(i)) = v[static_cast<std::size_t>(n)];
      Jp(j, static_cast<Eigen::Index>(i)) = 0.5 * (Jm1 - v[static_cast<std::size_t>(n + 1)]);
    }
  // Antiderivative F_jl(r) with F(0) = 0 (every integer order); for near-coincident pairs
  // (|a - b| r < 1e-2) the closed form cancels, use Gauss quadrature of int_0^r instead.
  auto quadF = [&](cplx a, cplx b, double r1) {
    std::vector<double> rq, Wq;
    radialRule({r1}, std::max(std::abs(a), std::abs(b)), 16, rq, Wq);
    cplx q = 0.0;
    for (std::size_t s = 0; s < rq.size(); ++s)
      q += besselJ(n, a * rq[s]) * besselJ(n, b * rq[s]) * rq[s] * Wq[s];
    return q;
  };
  MatC G = MatC::Zero(nk, nk);
  for (Eigen::Index j = 0; j < nk; ++j)
    for (Eigen::Index l = 0; l < nk; ++l) {
      const cplx a = grid.k(j), b = grid.k(l);
      cplx sum = 0.0, F0 = 0.0;
      for (std::size_t i = 0; i < nr; ++i) {
        const double r1 = g.radii[i];
        const auto ii = static_cast<Eigen::Index>(i);
        cplx F1;
        if (j == l) {
          const cplx x = a * r1, Ja = J(j, ii), Jpa = Jp(j, ii);
          F1 = 0.5 * r1 * r1 * (Jpa * Jpa + (1.0 - double(n * n) / (x * x)) * Ja * Ja);
        } else if (std::abs(a - b) * r1 < 1e-2) {
          F1 = quadF(a, b, r1);
        } else {
          F1 = r1 * (b * J(j, ii) * Jp(l, ii) - a * Jp(j, ii) * J(l, ii)) / (a * a - b * b);
        }
        sum += g.value[i] * (F1 - F0);
        F0 = F1;
      }
      G(j, l) = sg * sum * grid.w(l) * grid.k(l);
    }
  return G;
}

}  // namespace fg::cyl

namespace fg::cyl {

CylQ cylQ(const SpectralGrid& grid, int m, const Rings& rings, CylFactorization f) {
  FG_REQUIRE(rings.radii.size() == rings.eps.size(), "cylQ: one permittivity per ring");
  const Eigen::Index n = grid.k.size();
  const MatC Id = MatC::Identity(n, n);
  const int p = m + 1, q = m - 1;
  const cplx ei = rings.epsInf;
  RadialSteps d{rings.radii, {}};  // eps - eps_inf
  for (cplx e : rings.eps) d.value.push_back(e - ei);
  CylQ Q;
  Q.Qyy = ei * Id + galerkin(grid, m, m, d);  // E_y tangential at walls r = const: Laurent
  if (f == CylFactorization::Laurent) {
    Q.Qpp = ei * Id + galerkin(grid, p, p, d);
    Q.Qmm = ei * Id + galerkin(grid, q, q, d);
    Q.Qpm = MatC::Zero(n, n);
    Q.Qmp = MatC::Zero(n, n);
    return Q;
  }
  // Li in r (Plan eq:cyl_li): a = (rho - 1)/2, b = (eps - eps_inf)/2, rho = eps_inf/eps
  RadialSteps a{rings.radii, {}}, b{rings.radii, {}};
  for (cplx e : rings.eps) {
    a.value.push_back(0.5 * (ei / e - 1.0));
    b.value.push_back(0.5 * (e - ei));
  }
  MatC A(2 * n, 2 * n), B(2 * n, 2 * n);
  A << Id + galerkin(grid, p, p, a), galerkin(grid, p, q, a), galerkin(grid, q, p, a),
      Id + galerkin(grid, q, q, a);
  B << ei * Id + galerkin(grid, p, p, b), -galerkin(grid, p, q, b), -galerkin(grid, q, p, b),
      ei * Id + galerkin(grid, q, q, b);
  // Q = B A^{-1} via LU of A^T (no explicit inverse)
  const MatC QT = A.transpose().partialPivLu().solve(B.transpose());
  const MatC Qf = QT.transpose();
  Q.Qpp = Qf.topLeftCorner(n, n);
  Q.Qpm = Qf.topRightCorner(n, n);
  Q.Qmp = Qf.bottomLeftCorner(n, n);
  Q.Qmm = Qf.bottomRightCorner(n, n);
  Q.A = std::move(A);
  return Q;
}

}  // namespace fg::cyl

namespace fg::cyl {

void galerkinBar(const SpectralGrid& grid, int nu, int nup, const RadialSteps& g, const MatC& Gbar,
                 std::vector<cplx>& valueBar, std::vector<double>& radiiBar) {
  const std::size_t nr = g.radii.size();
  FG_REQUIRE(valueBar.size() == nr && radiiBar.size() == nr, "galerkinBar: sizes");
  const Eigen::Index nk = grid.k.size();
  RadialSteps unit{g.radii, std::vector<cplx>(nr, 0.0)};
  for (std::size_t i = 0; i < nr; ++i) {
    unit.value.assign(nr, 0.0);
    unit.value[i] = 1.0;
    valueBar[i] += (Gbar.array() * galerkin(grid, nu, nup, unit).conjugate().array()).sum();
  }
  auto J = [](int order, cplx z) {  // J_order(z), any integer order
    const int a = std::abs(order);
    return sgnOrder(order) * besselJAll(a, z)[static_cast<std::size_t>(a)];
  };
  for (std::size_t i = 0; i < nr; ++i) {
    const double r = g.radii[i];
    const cplx jump = g.value[i] - (i + 1 < nr ? g.value[i + 1] : cplx(0.0));
    VecC u(nk), v(nk);
    for (Eigen::Index j = 0; j < nk; ++j) {
      u(j) = J(nu, grid.k(j) * r);
      v(j) = grid.w(j) * grid.k(j) * J(nup, grid.k(j) * r) * r;
    }
    radiiBar[i] += std::real(jump * (u.transpose() * Gbar.conjugate() * v)(0));
  }
}

CylQBar cylQBar(const SpectralGrid& grid, int m, const Rings& rings, CylFactorization f, const CylQ& Q,
                const CylQ& Qbar) {
  const std::size_t nr = rings.radii.size();
  const Eigen::Index n = grid.k.size();
  const int p = m + 1, q = m - 1;
  const cplx ei = rings.epsInf;
  CylQBar out;
  out.eps.assign(nr, 0.0);
  out.radii.assign(nr, 0.0);
  // Q_yy = eps_inf 1 + [[d]]_{m<-m}, d = eps - eps_inf (both factorizations)
  RadialSteps d{rings.radii, {}};
  for (cplx e : rings.eps) d.value.push_back(e - ei);
  std::vector<cplx> db(nr, 0.0);
  out.epsInf += Qbar.Qyy.trace();
  galerkinBar(grid, m, m, d, Qbar.Qyy, db, out.radii);
  if (f == CylFactorization::Laurent) {
    out.epsInf += Qbar.Qpp.trace() + Qbar.Qmm.trace();
    galerkinBar(grid, p, p, d, Qbar.Qpp, db, out.radii);
    galerkinBar(grid, q, q, d, Qbar.Qmm, db, out.radii);
  } else {
    // Q = B A^{-1}: Bbar = Qbar A^{-H}, Abar = -Q^H Bbar
    MatC Qf(2 * n, 2 * n), Qfb(2 * n, 2 * n);
    Qf << Q.Qpp, Q.Qpm, Q.Qmp, Q.Qmm;
    Qfb << Qbar.Qpp, Qbar.Qpm, Qbar.Qmp, Qbar.Qmm;
    const MatC Bb = Q.A.partialPivLu().solve(MatC(Qfb.adjoint())).adjoint();
    const MatC Ab = -(Qf.adjoint() * Bb);
    RadialSteps a{rings.radii, {}}, b{rings.radii, {}};
    for (cplx e : rings.eps) {
      a.value.push_back(0.5 * (ei / e - 1.0));
      b.value.push_back(0.5 * (e - ei));
    }
    std::vector<cplx> ab(nr, 0.0), bb(nr, 0.0);
    auto blk = [&](const MatC& X, int r, int c) { return MatC(X.block(r * n, c * n, n, n)); };
    galerkinBar(grid, p, p, a, blk(Ab, 0, 0), ab, out.radii);
    galerkinBar(grid, p, q, a, blk(Ab, 0, 1), ab, out.radii);
    galerkinBar(grid, q, p, a, blk(Ab, 1, 0), ab, out.radii);
    galerkinBar(grid, q, q, a, blk(Ab, 1, 1), ab, out.radii);
    out.epsInf += blk(Bb, 0, 0).trace() + blk(Bb, 1, 1).trace();
    galerkinBar(grid, p, p, b, blk(Bb, 0, 0), bb, out.radii);
    galerkinBar(grid, p, q, b, MatC(-blk(Bb, 0, 1)), bb, out.radii);
    galerkinBar(grid, q, p, b, MatC(-blk(Bb, 1, 0)), bb, out.radii);
    galerkinBar(grid, q, q, b, blk(Bb, 1, 1), bb, out.radii);
    // a_i = (eps_inf / e_i - 1) / 2, b_i = (e_i - eps_inf) / 2 (holomorphic: zbar += conj(dw/dz) wbar)
    for (std::size_t i = 0; i < nr; ++i) {
      const cplx e = rings.eps[i];
      out.eps[i] += std::conj(-0.5 * ei / (e * e)) * ab[i] + 0.5 * bb[i];
      out.epsInf += std::conj(0.5 / e) * ab[i] - 0.5 * bb[i];
    }
  }
  for (std::size_t i = 0; i < nr; ++i) {
    out.eps[i] += db[i];
    out.epsInf -= db[i];
  }
  return out;
}

}  // namespace fg::cyl
