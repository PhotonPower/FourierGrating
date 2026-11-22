// FourierGrating - coupled azimuthal harmonics of a P-fold symmetric cross section (M21, Plan sections
// "Winkelentwicklung und Kopplung", "Gekoppelte Propagationsgleichung", eq:cylP_eps, eq:cylP_Q).
#pragma once

#include <vector>

#include "fg/core/types.hpp"
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

/// 1_{2L+1} (x) kappa: the node values repeated over the blocks (the derivative parts of Plan eq:cyl_M do not depend
/// on m); with coupledQ the argument of propagationMatrix.
VecC tiledK(const VecC& k, int blocks);

}  // namespace fg::cyl
