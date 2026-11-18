// FourierGrating - real-space fields of the Hankel basis (M20, Plan section "M20: S-Matrix, Loeser und
// Nahfeld"): synthesis of the spectra, radial quadrature for volume integrals, FFF-consistent
// reconstruction from the continuous auxiliary fields G_pm (Plan eq:cyl_G).
#pragma once

#include <cstddef>
#include <vector>

#include "fg/core/types.hpp"
#include "fg/cyl/spectrum.hpp"

namespace fg::cyl {

/// Synthesis matrices of harmonic m on the radii r_s (Nystroem, Plan sec:cylgrid):
/// T_nu(s, j) = w_j k_j J_nu(k_j r_s), so that f(r_s) = (T_nu f^)(s), for nu = m - 1, m, m + 1.
struct Synthesis {
  MatC Tm, T0, Tp;  ///< orders m - 1 (E_-, H_-), m (E_y, H_y), m + 1 (E_+, H_+)
};
Synthesis synthesis(const SpectralGrid& grid, int m, const std::vector<double>& r);

/// J_nu(k_j r) for nu = 0 .. nmax and all nodes at one radius: (nmax + 1) x N_k.
MatC besselTable(const VecC& k, int nmax, double r);

/// J_nu(k_j r) of any integer order from a besselTable row set (J_{-n} = (-1)^n J_n).
inline cplx tableJ(const MatC& table, int nu, Eigen::Index j) {
  const int n = nu < 0 ? -nu : nu;
  const cplx v = table(n, j);
  return (nu < 0 && n % 2 == 1) ? -v : v;
}

/// Composite Gauss-Legendre rule on [0, radii.back()] with breaks at the radii and panels
/// <= pi / (2 kmax) (Plan sec:cylgal), weights including r: int_0^R f(r) r dr ~ sum_s w_s f(r_s).
/// ring[s] is the ring (index into radii) that contains r_s.
struct RadialRule {
  std::vector<double> r, w;
  std::vector<std::size_t> ring;
};
RadialRule radialRule(const std::vector<double>& radii, double kmax, int nodes = 16);

/// E_pm from the continuous G_pm (Plan eq:cyl_G) with rho = eps_inf / eps at the point.
inline void fieldsFromG(cplx eps, cplx epsInf, cplx gp, cplx gm, cplx& ep, cplx& em) {
  const cplx rho = epsInf / eps;
  ep = 0.5 * ((rho + 1.0) * gp + (rho - 1.0) * gm);
  em = 0.5 * ((rho - 1.0) * gp + (rho + 1.0) * gm);
}

}  // namespace fg::cyl
