// FourierGrating - gradients of efficiency objectives by the discrete adjoint (Plan section "M19").
#pragma once

#include <functional>
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

/// A geometric parameter given by the structure at theta + delta (C++ API): the regions of contour2d,
/// polygon2d and sectioned (CAD) layers move with V = (r(theta + delta) - r(theta - delta)) / (2 delta),
/// Plan eq:shape2d per section height (eq:shape3d for CAD bodies: lateral moves only, no critical height
/// may move, and the segment structure of the sections must stay). at(delta) returns the whole stack.
struct StackParameter {
  std::string name;
  std::function<Stack(double delta)> at;
  double delta = 1e-6;
};
GradientResult gradient(const SolverOptions& opt, const Stack& stack, const Wavelength& wl, const Objective& objective,
                        const std::vector<StackParameter>& params);
GradientResult gradient(const std::string& configJson, const Wavelength& wl, const Objective& objective,
                        const std::vector<std::string>& params);

}  // namespace fg
