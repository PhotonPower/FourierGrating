// FourierGrating - gradients of efficiency objectives by the discrete adjoint (Plan section "M19").
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "fg/solvers/solver.hpp"

namespace fg {

/// An efficiency of the result: side 'R' (reflected) or 'T' (transmitted), order (n, m).
struct OrderRef {
  char side = 'R';
  int n = 0, m = 0;
};

/// L = sum of weight * efficiency (Plan "Schnittstellen"); as cotangents of the efficiencies the same
/// object defines a vector-Jacobian product.
struct Objective {
  std::vector<std::pair<OrderRef, double>> weights;
};

struct GradientResult {
  double value = 0;                ///< L
  std::vector<double> grad;        ///< dL / d(parameter), per parameter, in the units of the configuration
  Result result;                   ///< the forward result
  std::vector<std::string> notes;  ///< fixed choices (Plan "Methode": slices, frozen normal field)
};

/// Gradient of L with respect to numbers of the configuration (JSON pointers as in sweeps), first scan
/// point: "/layers/k/thickness_um", "/layers/k/breakpoints_um/j", "/layers/k/depth_um" and "/x0_um"
/// (sinus profile), "/layers/k/regions/..." and "/layers/k/inclusions/..." (any number of a contour2d or
/// polygon2d geometry: control points, weights, circle centre and radius, vertices), "/materials/name/n",
/// "/materials/name/n/0|1" (real/imaginary part), "/materials/name/eps", "/eps/0|1". RCWA and RK, periodic
/// structures without coordinate transformation; the slices, the RK step sequence and the raster normal
/// field of 2D cross-sections are fixed (profile normal fields are differentiated).
GradientResult gradient(const std::string& configJson, const Objective& objective, const std::vector<std::string>& params);
GradientResult gradient(const std::string& configJson, const Wavelength& wl, const Objective& objective,
                        const std::vector<std::string>& params);

}  // namespace fg
