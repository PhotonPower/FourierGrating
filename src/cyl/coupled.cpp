// FourierGrating - coupled harmonics of a P-fold symmetric cross section (M21, Plan eq:cylP_Q).
#include "fg/cyl/coupled.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/cyl/angular.hpp"
#include "fg/cyl/bessel.hpp"

namespace fg::cyl {

namespace {

/// Radial rule on [0, breaks.back()]: panels between the breaks, sub-panels of width <= hmax; a sub-panel that
/// starts (ends) at a break gets r = u_a + h (1 - cos t) (r = u_a + h sin t), t in [0, pi/2], so that dr/dt
/// vanishes like sqrt|r - r_b| there; a panel between two breaks has at least two sub-panels.
void radialRule(const std::vector<double>& breaks, double hmax, int nodes, std::vector<double>& r,
                std::vector<double>& W) {
  const auto& gl = quadrature::gaussLegendre(nodes);
  for (std::size_t e = 0; e + 1 < breaks.size(); ++e) {
    const double a = breaks[e], b = breaks[e + 1];
    if (b <= a) continue;
    const int nsub = std::max(2, static_cast<int>(std::ceil((b - a) / hmax - 1e-12)));
    const double h = (b - a) / nsub;
    for (int s = 0; s < nsub; ++s) {
      const double ua = a + s * h;
      const bool lo = s == 0, hi = s == nsub - 1;
      const double t0 = (lo || hi) ? 0.0 : -1.0, t1 = (lo || hi) ? PI / 2 : 1.0;
      for (std::size_t q = 0; q < gl.x.size(); ++q) {
        const double t = 0.5 * (t1 - t0) * gl.x[q] + 0.5 * (t1 + t0), wt = 0.5 * (t1 - t0) * gl.w[q];
        if (lo) {
          r.push_back(ua + h * (1 - std::cos(t)));
          W.push_back(wt * h * std::sin(t));
        } else if (hi) {
          r.push_back(ua + h * std::sin(t));
          W.push_back(wt * h * std::cos(t));
        } else {
          r.push_back(ua + 0.5 * h * (t + 1));
          W.push_back(wt * 0.5 * h);
        }
      }
    }
  }
}

}  // namespace

VecC tiledK(const VecC& k, int blocks) {
  FG_REQUIRE(blocks >= 1, "tiledK: blocks >= 1");
  VecC t(k.size() * blocks);
  for (int b = 0; b < blocks; ++b) t.segment(b * k.size(), k.size()) = k;
  return t;
}

namespace {

/// The radial rule of a class and d_q(r_i) r_i W_i = sum_b (eps_b - epsBg) chi_q^(b)(r_i) r_i W_i for
/// q = -qmax..qmax (row qmax + q; chi_{-q} = conj chi_q). The panels are as wide as allowed by kabs.
struct RadialTable {
  std::vector<double> r, rw;  ///< nodes, r W
  MatC d;    ///< contrast eps - epsBg times r W
  MatC dIm;  ///< Im eps (with Im epsBg in q = 0) times r W
  int qmax = 0;
};

RadialTable radialTable(const HarmonicMap& map, const AngularSection& section, double kabs, int nodesPerPanel) {
  FG_REQUIRE(map.P >= 1 && map.L >= 0, "coupledQ: P >= 1, L >= 0");
  FG_REQUIRE(section.regions.size() == section.eps.size(), "coupledQ: one permittivity per region");
  RadialTable t;
  t.qmax = 2 * map.L;
  // critical radii of all regions: chi_q is smooth between them
  std::vector<double> breaks{0.0};
  for (const auto& reg : section.regions) {
    const auto c = criticalRadii(reg);
    breaks.insert(breaks.end(), c.begin(), c.end());
  }
  std::sort(breaks.begin(), breaks.end());
  const double tol = 1e-12 * std::max(breaks.back(), 1e-300);
  breaks.erase(std::unique(breaks.begin(), breaks.end(), [&](double a, double b) { return b - a <= tol; }),
               breaks.end());
  std::vector<double> W;
  radialRule(breaks, PI / (2 * std::max(kabs, 1e-300)), nodesPerPanel, t.r, W);
  const auto nr = static_cast<Eigen::Index>(t.r.size());
  const int qmax = t.qmax;
  t.d = MatC::Zero(2 * qmax + 1, nr);
  t.dIm = MatC::Zero(2 * qmax + 1, nr);
  for (Eigen::Index i = 0; i < nr; ++i) {
    const double rw = t.r[static_cast<std::size_t>(i)] * W[static_cast<std::size_t>(i)];
    t.rw.push_back(rw);
    t.dIm(qmax, i) += section.epsBg.imag() * rw;
    for (std::size_t b = 0; b < section.regions.size(); ++b) {
      const cplx de = section.eps[b] - section.epsBg;
      const VecC chi = angularCoefficients(section.regions[b], t.r[static_cast<std::size_t>(i)], map.P, qmax);
      for (int q = 0; q <= qmax; ++q) {
        t.d(qmax + q, i) += de * chi(q) * rw;
        t.dIm(qmax + q, i) += de.imag() * chi(q) * rw;
        if (q) {
          t.d(qmax - q, i) += de * std::conj(chi(q)) * rw;
          t.dIm(qmax - q, i) += de.imag() * std::conj(chi(q)) * rw;
        }
      }
    }
  }
  return t;
}

/// Bessel tables T_nu(j, i) = J_nu(k_j r_i), nu = nuMin..nuMax.
std::vector<MatC> besselTables(const VecC& k, const std::vector<double>& r, int nuMin, int nuMax) {
  const Eigen::Index nk = k.size(), nr = static_cast<Eigen::Index>(r.size());
  const int aMax = std::max(std::abs(nuMin), std::abs(nuMax));
  std::vector<MatC> T(static_cast<std::size_t>(nuMax - nuMin + 1), MatC(nk, nr));
  for (Eigen::Index j = 0; j < nk; ++j)
    for (Eigen::Index i = 0; i < nr; ++i) {
      const auto v = besselJAll(aMax, k(j) * r[static_cast<std::size_t>(i)]);
      for (int nu = nuMin; nu <= nuMax; ++nu) {
        const int a = std::abs(nu);
        const double sg = (nu < 0 && a % 2) ? -1.0 : 1.0;  // J_{-n} = (-1)^n J_n
        T[static_cast<std::size_t>(nu - nuMin)](j, i) = sg * v[static_cast<std::size_t>(a)];
      }
    }
  return T;
}

}  // namespace

CylQ coupledQ(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section, int nodesPerPanel) {
  const Eigen::Index nk = grid.k.size();
  const int nb = map.blocks();
  const RadialTable rt = radialTable(map, section, grid.k.cwiseAbs().maxCoeff(), nodesPerPanel);
  const int qmax = rt.qmax, nuMin = map.m(0) - 1, nuMax = map.m(nb - 1) + 1;
  const std::vector<MatC> T = besselTables(grid.k, rt.r, nuMin, nuMax);
  auto table = [&](int nu) -> const MatC& { return T[static_cast<std::size_t>(nu - nuMin)]; };
  const VecC wk = grid.w.cwiseProduct(grid.k);
  const Eigen::Index nf = nb * nk;
  auto assemble = [&](int s) {
    MatC Q = MatC::Zero(nf, nf);
    for (int b = 0; b < nb; ++b)
      for (int bp = 0; bp < nb; ++bp) {
        const VecC dq = rt.d.row(qmax + (b - bp)).transpose();
        Q.block(b * nk, bp * nk, nk, nk) =
            (table(map.m(b) + s) * dq.asDiagonal() * table(map.m(bp) + s).transpose()) * wk.asDiagonal();
        if (b == bp) Q.block(b * nk, b * nk, nk, nk).diagonal().array() += section.epsBg;
      }
    return Q;
  };
  CylQ Q;
  Q.Qpp = assemble(+1);
  Q.Qmm = assemble(-1);
  Q.Qyy = assemble(0);
  Q.Qpm = MatC::Zero(nf, nf);
  Q.Qmp = MatC::Zero(nf, nf);
  return Q;
}

PlaneWaveSource coupledSource(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                              double kpar, double Psi, int nodesPerPanel) {
  const Eigen::Index nk = grid.k.size();
  const int nb = map.blocks();
  const RadialTable rt =
      radialTable(map, section, std::max(grid.k.cwiseAbs().maxCoeff(), std::abs(kpar)), nodesPerPanel);
  const int qmax = rt.qmax, nuMin = map.m(0) - 1, nuMax = map.m(nb - 1) + 1;
  const std::vector<MatC> T = besselTables(grid.k, rt.r, nuMin, nuMax);
  const std::vector<MatC> Tk = besselTables(VecC::Constant(1, kpar), rt.r, nuMin, nuMax);
  auto idx = [&](int nu) { return static_cast<std::size_t>(nu - nuMin); };
  // i^nu e^{-i nu Psi} (Plan eq:cyl_jacobi)
  auto phase = [Psi](int nu) { return std::polar(1.0, nu * (0.5 * PI - Psi)); };
  const Eigen::Index nf = nb * nk;
  // block l: sum_l' phase(nu_l') int d_{l-l'} J_{nu_l}(k_j r) J_{nu_l'}(kpar r) r dr, nu_l = m_l + s
  auto column = [&](int s) {
    VecC c = VecC::Zero(nf);
    for (int b = 0; b < nb; ++b) {
      VecC g = VecC::Zero(static_cast<Eigen::Index>(rt.r.size()));
      for (int bp = 0; bp < nb; ++bp) {
        const int nup = map.m(bp) + s;
        g += phase(nup) * rt.d.row(qmax + (b - bp)).transpose().cwiseProduct(Tk[idx(nup)].row(0).transpose());
      }
      c.segment(b * nk, nk) = T[idx(map.m(b) + s)] * g;
    }
    return c;
  };
  PlaneWaveSource B;
  B.By = column(0);
  B.Bpm = MatC::Zero(2 * nf, 2);
  B.Bpm.col(0).head(nf) = column(+1);
  B.Bpm.col(1).tail(nf) = column(-1);
  return B;
}

SectionVolume sectionVolume(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                            int nodesPerPanel) {
  RadialTable t = radialTable(map, section, grid.k.cwiseAbs().maxCoeff(), nodesPerPanel);
  SectionVolume v;
  v.r = std::move(t.r);
  v.imEps = std::move(t.dIm);
  v.contrast = std::move(t.d);
  v.qmax = t.qmax;
  return v;
}

CoupledQBar coupledQBar(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                        const CylQ& Qbar, int nodesPerPanel) {
  const Eigen::Index nk = grid.k.size();
  const int nb = map.blocks();
  FG_REQUIRE(Qbar.Qpp.rows() == nb * nk && Qbar.Qmm.rows() == nb * nk && Qbar.Qyy.rows() == nb * nk,
             "coupledQBar: sizes");
  RadialTable rt = radialTable(map, section, grid.k.cwiseAbs().maxCoeff(), nodesPerPanel);
  const int qmax = rt.qmax, nuMin = map.m(0) - 1, nuMax = map.m(nb - 1) + 1;
  const auto nr = static_cast<Eigen::Index>(rt.r.size());
  const std::vector<MatC> T = besselTables(grid.k, rt.r, nuMin, nuMax);
  auto table = [&](int nu) -> const MatC& { return T[static_cast<std::size_t>(nu - nuMin)]; };
  const VecC wk = grid.w.cwiseProduct(grid.k);
  CoupledQBar out;
  out.qmax = qmax;
  out.r = rt.r;
  out.rw = rt.rw;
  // block (b, b') of Q_ss = T_a diag(d_{b-b'} r W) T_b'^T diag(w k): tr(X^H T_a D B) = sum_i D_ii (B X^H T_a)_ii
  out.G = MatC::Zero(2 * qmax + 1, nr);
  for (const auto& [Qb, s] : {std::pair<const MatC*, int>{&Qbar.Qpp, 1}, {&Qbar.Qmm, -1}, {&Qbar.Qyy, 0}}) {
    for (int b = 0; b < nb; ++b)
      for (int bp = 0; bp < nb; ++bp) {
        const MatC Y = Qb->block(b * nk, bp * nk, nk, nk).adjoint() * table(map.m(b) + s);  // nk x nr
        const MatC& Tb = table(map.m(bp) + s);
        for (Eigen::Index i = 0; i < nr; ++i)
          out.G(qmax + (b - bp), i) += Tb.col(i).cwiseProduct(wk).cwiseProduct(Y.col(i)).sum();
        if (b == bp) out.epsBg += Qb->block(b * nk, b * nk, nk, nk).trace();
      }
  }
  // d_q(r_i) r_i W_i = sum_b (eps_b - epsBg) chi_q^(b)(r_i) r_i W_i: dL = Re sum G d(d); bar = conj(sum G chi r W)
  for (std::size_t reg = 0; reg < section.regions.size(); ++reg) {
    cplx g = 0;
    for (Eigen::Index i = 0; i < nr; ++i) {
      const VecC chi = angularCoefficients(section.regions[reg], rt.r[static_cast<std::size_t>(i)], map.P, qmax);
      cplx gi = out.G(qmax, i) * chi(0);
      for (int q = 1; q <= qmax; ++q) gi += out.G(qmax + q, i) * chi(q) + out.G(qmax - q, i) * std::conj(chi(q));
      g += gi * rt.rw[static_cast<std::size_t>(i)];
    }
    out.eps.push_back(std::conj(g));
    out.epsBg -= std::conj(g);
  }
  return out;
}

double coupledQShape(const CoupledQBar& bar, const HarmonicMap& map, const AngularSection& section,
                     const std::vector<BoundaryVelocity>& velocity) {
  FG_REQUIRE(velocity.size() == section.regions.size(), "coupledQShape: one velocity per region");
  const int qmax = bar.qmax;
  double dL = 0;
  for (std::size_t reg = 0; reg < section.regions.size(); ++reg) {
    if (!velocity[reg]) continue;
    const cplx de = section.eps[reg] - section.epsBg;
    cplx g = 0;
    for (std::size_t i = 0; i < bar.r.size(); ++i) {
      const auto ii = static_cast<Eigen::Index>(i);
      const VecC dchi = angularCoefficientsShape(section.regions[reg], bar.r[i], map.P, qmax, velocity[reg]);
      cplx gi = bar.G(qmax, ii) * dchi(0);
      for (int q = 1; q <= qmax; ++q) gi += bar.G(qmax + q, ii) * dchi(q) + bar.G(qmax - q, ii) * std::conj(dchi(q));
      g += gi * bar.rw[i];
    }
    dL += std::real(de * g);
  }
  return dL;
}

}  // namespace fg::cyl
