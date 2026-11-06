// FourierGrating - complete structure: substrate (bottom), layers, superstrate (top).
#pragma once

#include <map>
#include <string>
#include <vector>

#include "fg/geometry/coordinate_transform.hpp"
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

  bool is2D() const { return dz_um > 0.0; }
  const Material& material(const std::string& name) const;
};

}  // namespace fg
