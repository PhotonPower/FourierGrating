// FourierGrating - what the stack assembly needs from a field representation (Plan section "M20:
// Architektur und Aufwand"): the Cartesian Fourier basis and the cylindrical Hankel basis provide
// the propagation matrices, layer modes and half-space bases; RK, RCWA, S-matrices and the
// star-product chain are shared.
#pragma once

#include <optional>
#include <vector>

#include "fg/core/types.hpp"
#include "fg/physics/modes.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::detail {

/// One layer of the stack as seen by the assembly (micrometres, y relative to the layer bottom).
class StackLayer {
 public:
  virtual ~StackLayer() = default;
  virtual double thickness_um() const = 0;
  virtual bool homogeneous() const { return false; }
  virtual bool yInvariant() const { return true; }
  /// Number of RCWA staircase slices (y-invariant layers: 1).
  int rcwaSlices = 1;
  /// Interior heights 0 < y < thickness (micrometres, ascending) where Q(y) is not smooth: RK steps and
  /// Chebyshev interpolants end there (Plan section "M24: Kritische Höhen und Bänder"). Empty for most layers.
  virtual std::vector<double> yBreaks_um() const { return {}; }
};

/// Field representation of the assembly: F = Psi V in the reference medium (superstrate) and the
/// substrate (Plan eq:Stot), dF/dy = M(y) F in the layers.
class ModalBasis {
 public:
  virtual ~ModalBasis() = default;
  virtual const ModeBasis& ref() const = 0;
  virtual const ModeBasis& sub() const = 0;
  /// Propagation matrix of a layer at height y_um (relative to its bottom).
  virtual MatC M(const StackLayer& layer, double y_um) const = 0;
  /// Eigen-decomposition of a slice matrix M (RCWA).
  virtual LayerModes modes(const MatC& M) const = 0;
  /// Modes of a homogeneous layer (analytic where available).
  virtual LayerModes homogeneousModes(const StackLayer& layer) const = 0;
  /// Estimate of the spectral radius of M in a layer (slice criterion, Plan eq:slicecrit).
  virtual double rhoEstimate(const StackLayer& layer) const = 0;
  /// Reduced 2 n_f system of one polarization (Cartesian 1D only, Plan eq:TE / eq:TM).
  virtual std::optional<Polarization> polarization() const { return std::nullopt; }
  /// G > 0: the half-space bases and the modes of homogeneous layers are block-diagonal over G groups (state rows
  /// {g + k G}, k = 0..3; mode columns {g, G + g}), so that their S-matrices are grouped (groupedLayerSMatrix,
  /// groupedFromTransfer, groupedStar). Hankel basis: one group per spectral node and harmonic (M21 acceleration).
  virtual Eigen::Index groups() const { return 0; }
};

}  // namespace fg::detail
