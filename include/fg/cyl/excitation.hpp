// FourierGrating - plane-wave excitation of the Hankel basis (M20, Plan section "M20: Anregung",
// paragraph "Ebene Welle (Streufeld)").
#pragma once

#include <Eigen/LU>

#include "fg/core/types.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/spectrum.hpp"

namespace fg::cyl {

/// int_0^inf g(r) J_nu(k_j r) J_nu'(kp r) r dr of a compact piecewise-constant g: the column of Plan
/// eq:cyl_gal at the extra node k' = kp without weight (Lommel integrals, Plan eq:lommel, for
/// |nu| = |nu'|, Gauss quadrature otherwise).
VecC besselProjection(const SpectralGrid& grid, int nu, int nup, const RadialSteps& g, cplx kp);

/// Source spectra of harmonic m for a background plane wave with lateral wave number kpar along the
/// azimuth Psi (Plan eq:cyl_jacobi): (p_+; p_-) = Bpm (c_+; c_-), p_y = e_y By (c_pm = e_x +- i e_z,
/// e_y: background field on the axis). Laurent: P = (eps - eps_inf) E_bg; Li in r: the source
/// consistent with the factorization, P = Q E_bg - eps_inf E_bg = b^ - Q a^ (Q = B A^{-1} of Plan
/// eq:cyl_li, a^, b^ the columns of [[a]], [[b]] at k' = kpar; the delta parts cancel since
/// eps_bg = eps_inf).
struct PlaneWaveSource {
  MatC Bpm;  ///< 2 N_k x 2
  VecC By;
};
PlaneWaveSource planeWaveSource(const SpectralGrid& grid, int m, const Rings& rings, double kpar, double Psi,
                                CylFactorization f, const CylQ& Q);

/// sigma = (kappa E p_y, -kappa E p_y, p_+, -p_-) of Plan eq:cyl_scat with E = Q_yy^{-1} (LU of Q_yy).
VecC scatteredSource(const VecC& k, const Eigen::PartialPivLU<MatC>& Qyy, const PlaneWaveSource& B, cplx cp,
                     cplx cm, cplx ey);

/// Harmonic cut-off m_cut = ceil(x + 4 x^{1/3} + 2) (Wiscombe; Plan eq:cyl_scat).
int harmonicCutoff(double x);

}  // namespace fg::cyl
