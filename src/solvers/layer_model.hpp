// FourierGrating - internal layer models shared by the solvers and the near-field reconstruction.
// Units: x, y in micrometres (y relative to the layer bottom); Q blocks in normalized form.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/geometry/stack.hpp"
#include "fg/physics/qmatrix.hpp"

namespace fg::detail {

class LayerModel {
 public:
  virtual ~LayerModel() = default;
  virtual double thickness_um() const = 0;
  virtual bool homogeneous() const { return false; }
  virtual bool yInvariant() const { return true; }
  /// Homogeneous layers: the material.
  virtual const Material& material() const;
  /// Factorized permittivity blocks at height y (Plan section "Faktorisierung").
  virtual QBlocks Q(double y_um) const = 0;
  /// Pointwise permittivity tensor at (x, y).
  virtual Mat3C epsAt(double x_um, double y_um) const = 0;
  /// Normal vector field used by the factorization (FFF-consistent reconstruction).
  virtual Vec3D normalAt(double x_um) const = 0;
  /// Fourier coefficients of N_i N_j (FFF-consistent reconstruction).
  virtual const NormalFieldCoefficients& NN() const = 0;
  /// Permittivities occurring in the layer (spectral radius estimate).
  std::vector<cplx> epsList;
  /// Number of RCWA staircase slices (y-invariant layers: 1).
  int rcwaSlices = 1;
  std::vector<std::string> warnings;
};

/// Builds the models of all layers of a stack (bottom to top) for a given truncation.
std::vector<std::unique_ptr<LayerModel>> buildLayerModels(const Stack& stack, const IndexMap& map,
                                                          int defaultSlices);

}  // namespace fg::detail
