// FourierGrating - 1D surface-relief profiles y = g(x) (Plan section "Schichtmodelle").
#pragma once

#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/geometry/layer.hpp"
#include "fg/geometry/normal_field.hpp"

namespace fg {

/// g(x) in micrometres (x in [0, dx) or any real x, periodic).
double profileValue(const Profile1DLayer& p, double x, double dx);
/// g'(x).
double profileSlope(const Profile1DLayer& p, double x, double dx);
/// Sorted crossings x in [0, dx) with g(x) = y (y in micrometres, 0 < y < depth).
std::vector<double> profileCrossings(const Profile1DLayer& p, double y, double dx);
/// Breakpoints (0 = x0 < ... < xK = dx) and flags 'below' for each interval at height y.
void profileIntervals(const Profile1DLayer& p, double y, double dx, std::vector<double>& breaks,
                      std::vector<bool>& below);

}  // namespace fg
