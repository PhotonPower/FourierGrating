// FourierGrating - coupled azimuthal harmonics of a P-fold symmetric cross section (M21, Plan sections
// "Winkelentwicklung und Kopplung", "Gekoppelte Propagationsgleichung", eq:cylP_eps, eq:cylP_Q).
#pragma once

#include <vector>

#include "fg/core/types.hpp"
#include "fg/cyl/angular.hpp"
#include "fg/cyl/excitation.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/geometry/contour.hpp"

namespace fg::cyl {

/// Harmonics m_l = m0 + l P, |l| <= L, of one symmetry class (Plan sec:cylP). Block b = l + L; the state index of
/// (block b, node j) is b N_k + j (counterpart of IndexMap).
struct HarmonicMap {
  int m0 = 0, P = 1, L = 0;
  int blocks() const { return 2 * L + 1; }
  int m(int b) const { return m0 + (b - L) * P; }
};

/// Cross section with angular dependence, normalized coordinates (k0 x, k0 z), centred on the axis: regions with
/// permittivities eps on the background epsBg, which continues outside (Plan sec:cylP: eps_q is compact for q != 0,
/// eps_0 = background + compact rest).
struct AngularSection {
  std::vector<geometry::Region> regions;  ///< normalized (Region::normalize)
  std::vector<cplx> eps;                  ///< one per region
  cplx epsBg = 1.0;
};

/// Laurent blocks of Plan eq:cylP_Q for one class, n_f = (2L+1) N_k: (Q_ss)_{ll'} = epsBg delta_ll' 1 +
/// sum_b (eps_b - epsBg) [[chi^(b)_{l-l'}]]_{m_l+s <- m_l'+s} for s = +1 (Qpp), -1 (Qmm), 0 (Qyy); Qpm = Qmp = 0.
/// The Galerkin integrals run over the panels between the critical radii of all regions (width <= pi/(2 max|k|),
/// 16 Gauss nodes, square-root substitution at panel ends that are critical radii: chi ~ sqrt|r - r_t| at tangent
/// radii, Plan sec:cylP). nodesPerPanel overrides the 16 (convergence checks).
CylQ coupledQ(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section, int nodesPerPanel = 16);

/// Source spectra of one class for a background plane wave with lateral wave number kpar along the azimuth Psi (Plan
/// sec:cylP "Anregung und Symmetrieklassen", Laurent): block l of (p_+; p_-) = Bpm (c_+; c_-) and p_y = e_y By is
/// sum_{|l'| <= L} [[eps_{l-l'} - epsBg delta]]_{nu_l <- nu_l'} at the extra node k' = kpar (without weight) times
/// the Jacobi-Anger phase i^nu' e^{-i nu' Psi}, nu = m + 1, m - 1, m. The background harmonics of the class are
/// truncated like the field (|l'| <= L), consistent with coupledQ. Radial rule as coupledQ.
PlaneWaveSource coupledSource(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                              double kpar, double Psi, int nodesPerPanel = 16);

/// Volume integrals over an angular layer (M21, Plan sec:cylP "Volumenintegrale"): the radial rule of coupledQ (same
/// panels and nodes, so that the integrals are consistent with the Galerkin blocks) and per node r_i the angular
/// coefficients of Im eps (q = 0 includes Im epsBg on [0, largest critical radius]) and of the contrast eps - epsBg,
/// both times r_i W_i, row qmax + q, q = -qmax..qmax, qmax = 2L. With E = sum_l E_l(r) e^{i m_l phi}:
/// int_0^{2 pi} int Im eps |E|^2 r dr dphi = 2 pi sum_i sum_{l,l'} conj(E_l(r_i)) imEps(qmax + l - l', i) E_l'(r_i).
struct SectionVolume {
  std::vector<double> r;
  MatC imEps, contrast;
  int qmax = 0;
};
SectionVolume sectionVolume(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                            int nodesPerPanel = 16);

/// Reverse of coupledQ (M21 gradients, V22g) for the cotangent Qbar (Qpp, Qmm, Qyy; Qpm = Qmp = 0 in the forward
/// map), dL = Re sum conj(Qbar) dQ: the cotangents of the region permittivities and of epsBg (dL = Re conj(bar) d eps)
/// and, for coupledQShape, G(qmax + q, i) = dL / d(d_q(r_i) r_i W_i) on the radial rule of coupledQ (r, r W).
struct CoupledQBar {
  std::vector<cplx> eps;
  cplx epsBg = 0.0;
  std::vector<double> r, rw;
  MatC G;
  int qmax = 0;
};
CoupledQBar coupledQBar(const SpectralGrid& grid, const HarmonicMap& map, const AngularSection& section,
                        const CylQ& Qbar, int nodesPerPanel = 16);

/// Shape derivative dL/dtheta of coupledQ for regions moving with velocity[b] (empty: region b fixed): sum_b
/// (eps_b - epsBg) d chi_q^(b) / d theta (angularCoefficientsShape) on the radial rule of bar. chi_q is continuous in
/// r, so the moving critical radii add no terms; the derivative is that of the exact integrals (the quadrature of
/// coupledQ resolves them to about 1e-12, the square-root substitution covers d chi ~ |r - r_t|^{-1/2}).
double coupledQShape(const CoupledQBar& bar, const HarmonicMap& map, const AngularSection& section,
                     const std::vector<BoundaryVelocity>& velocity);

/// 1_{2L+1} (x) kappa: the node values repeated over the blocks (the derivative parts of Plan eq:cyl_M do not depend
/// on m); with coupledQ the argument of propagationMatrix.
VecC tiledK(const VecC& k, int blocks);

}  // namespace fg::cyl
