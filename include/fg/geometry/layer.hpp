// FourierGrating - layer models (Plan section "Schichtmodelle und Scheibenzerlegung").
// All lengths here are in micrometres; conversion to normalized units happens in the solver setup.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <variant>
#include <vector>

#include "fg/geometry/material.hpp"

namespace fg {

enum class Factorization { Auto, Laurent, LiLamellar, FFF_iso, FFF_A, FFF_B };

struct HomogeneousLayer {
  double thickness_um = 0;
  std::string material;
};

struct Lamellar1DLayer {            ///< y-invariant, piecewise constant in x
  double thickness_um = 0;
  std::vector<double> breakpoints_um;  ///< 0 = x0 < x1 < ... < xK = dx
  std::vector<std::string> materials;  ///< K entries
};

struct Profile1DLayer {             ///< region below y = g(x) is 'below', above is 'above'
  double depth_um = 0;
  std::function<double(double)> g;     ///< g(x) in [0, depth], period dx
  std::function<double(double)> dg;    ///< g'(x) (normal field), optional
  std::string below, above;
  int slices = 0;                      ///< 0 = automatic (Plan eq:slicecrit)
};

struct Polygon2DLayer {             ///< y-invariant cross-section in the (x, z) cell
  double thickness_um = 0;
  std::string background;
  struct Inclusion { std::vector<std::array<double, 2>> vertices; std::string material; };
  std::vector<Inclusion> inclusions;   ///< non-overlapping, counter-clockwise
};

struct Profile2DLayer {
  double depth_um = 0;
  std::function<double(double, double)> g;
  std::string below, above;
  int slices = 0;
};

using LayerSpec = std::variant<HomogeneousLayer, Lamellar1DLayer, Profile1DLayer, Polygon2DLayer,
                               Profile2DLayer>;

struct LayerEntry {
  LayerSpec spec;
  Factorization factorization = Factorization::Auto;
};

}  // namespace fg
