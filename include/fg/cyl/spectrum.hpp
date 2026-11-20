// FourierGrating - discrete Hankel spectrum (M20, Plan section "Diskretes Spektrum").
#pragma once

#include <vector>

#include "fg/core/types.hpp"

namespace fg::cyl {

/// Options of the spectral grid (Plan sec:cylgrid); normalized units.
struct SpectrumOptions {
  double kmax = 0.0;                    ///< cut-off k_max = kappa_max max Re sqrt(eps)
  double Rs = 0.0;                      ///< outer radius of the structure (sets the panel width)
  double dkFactor = 1.0;                ///< panel width pi / (2 Rs dkFactor)
  int nodesPerPanel = 8;                ///< Gauss-Legendre nodes per panel
  std::vector<double> branchPoints;     ///< Re sqrt(eps_h) of the half spaces (beta_h = 0)
  double deformA = 0.0, deformB = 0.0;  ///< deformed section [k_a, k_b] (empty: real path)
  double delta = 0.0;                   ///< depth of the deformation (0: real path)
  /// poles p (Im p >= 0) close to the real section of the path, e.g. plasmons that leak between the branch
  /// points (V21h): graded panels with edges Re p +- w 2^j, w = max(Im p, 1e-3 h_max), up to the panel width
  std::vector<cplx> resonances;
  /// real poles of the guided modes of a lossless planar background on the deformed section (V21d, Plan
  /// v2.38): graded panels with edges p +- w 2^j, w = max(distance of p to the path, 1e-3 h_max)
  std::vector<double> boundModes;
};

/// Nodes k_j (complex on the deformed section) and weights w_j including dk/dt, so that
/// int_path f(k) dk ~ sum_j w_j f(k_j) (Nystroem, Plan sec:cylgrid).
struct SpectralGrid {
  VecC k;
  VecC w;
  std::vector<double> edges;  ///< panel boundaries in Re k (branch points, deformation ends, k_max)
};

/// Composite Gauss-Legendre grid on [0, k_max] with square-root substitutions at the branch points
/// and the path k(u) = u - i delta sin(pi (u - k_a)/(k_b - k_a)) on [k_a, k_b]; graded panels around the
/// resonances on the real section (Plan sec:cylgrid).
SpectralGrid spectralGrid(const SpectrumOptions& opt);

/// Guided modes (TE and TM) of a lossless planar stack (layers bottom to top, normalized thicknesses k0 h)
/// between the half spaces: the real zeros k in (k_a, k_b) of the characteristic function of the layer
/// matrices (Plan sec:cylgrid, v2.38), sorted; empty if a material is lossy.
std::vector<double> planarBoundModes(cplx epsSub, cplx epsSup, const std::vector<cplx>& epsLayers,
                                     const std::vector<double>& thickness, double ka, double kb);

}  // namespace fg::cyl
