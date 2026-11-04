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

/// FFF_iso: book (IV.32), Q = [[eps]] - Delta [[NN]]; FFF_sym: reciprocity-preserving
/// Q = [[eps]] - (Delta [[NN]] + [[NN]] Delta)/2 (Plan eq:Qsym, default for isotropic profiles).
enum class Factorization { Auto, Laurent, LiLamellar, FFF_iso, FFF_sym, FFF_rig, FFF_A, FFF_B };

struct HomogeneousLayer {
  double thickness_um = 0;
  std::string material;
};

struct Lamellar1DLayer {            ///< y-invariant, piecewise constant in x
  double thickness_um = 0;
  std::vector<double> breakpoints_um;  ///< 0 = x0 < x1 < ... < xK = dx
  std::vector<std::string> materials;  ///< K entries
};

struct Profile1DLayer {             ///< region y < g(x) is 'below', y > g(x) is 'above'
  enum class Shape { Sinus, PiecewiseLinear, Function };
  double depth_um = 0;                 ///< a = max g - min g; g in [0, a]
  Shape shape = Shape::Sinus;
  double x0_um = 0;                    ///< Sinus: g = a/2 (1 + cos(K (x - x0)))
  std::vector<double> xs_um, ys_um;    ///< PiecewiseLinear: periodic polygon, xs ascending in [0, dx)
  std::function<double(double)> g;     ///< Function: g(x) in [0, a], period dx
  std::function<double(double)> dg;    ///< Function: g'(x) (optional, else finite differences)
  std::string below, above;
  int slices = 0;                      ///< RCWA staircase slices (0 = solver default)
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
