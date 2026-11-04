// FourierGrating - normal vector fields for FFF (Plan section "Normalen- und Tangentenfelder").
#pragma once

#include "fg/core/index_map.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/geometry/layer.hpp"

namespace fg {

/// Fourier coefficient grids of N_i N_j (i, j in x, y, z), |p| <= 2Nx, |q| <= 2Nz.
struct NormalFieldCoefficients {
  std::array<std::array<fourier::CoefficientGrid, 3>, 3> NN;
};

/// 1D profile: N = (-g', 1, 0)/sqrt(1+g'^2) (Plan eq:N1d). TODO(M3).
NormalFieldCoefficients normalField(const Profile1DLayer& layer, double dx_norm, const IndexMap& map);
/// y-invariant 2D cross-section: Algorithm NF (signed distance, smoothing). TODO(M7).
NormalFieldCoefficients normalField(const Polygon2DLayer& layer, double dx_norm, double dz_norm,
                                    const IndexMap& map);

}  // namespace fg
