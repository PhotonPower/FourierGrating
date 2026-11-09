// FourierGrating - complete structure: substrate (bottom), layers, superstrate (top).
#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "fg/geometry/aperiodic.hpp"
#include "fg/geometry/coordinate_transform.hpp"
#include "fg/geometry/cylindrical.hpp"
#include "fg/geometry/layer.hpp"

namespace fg {

struct Stack {
  double dx_um = 1.0;
  double dz_um = 0.0;                       ///< 0 for 1D gratings
  std::map<std::string, Material> materials;
  std::string substrate = "substrate";
  std::string superstrate = "superstrate";
  std::vector<LayerEntry> layers;           ///< ordered from bottom (substrate side) to top
  CoordinateTransformSpec transform;        ///< coordinate transformation x = f(X), z = g(Z) (M17)
  AperiodicSpec aperiodic;                  ///< isolated structure with PML zones (M18)
  /// Structured half spaces (M18, Plan "Einfall (b)"): y-invariant cross-sections (thickness ignored),
  /// e.g. a waveguide propagating along y; superstrate/substrate then name the cladding (background).
  std::optional<LayerEntry> superstrateStructure, substrateStructure;
  /// Rotationally symmetric structure about the y axis (M20); dx_um, dz_um are then unused.
  CylindricalSpec cylindrical;

  bool is2D() const { return dz_um > 0.0; }
  const Material& material(const std::string& name) const;
};

}  // namespace fg
