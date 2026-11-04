// FourierGrating - normal vector fields for FFF (Plan section "Normalen- und Tangentenfelder").
#pragma once

#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/geometry/layer.hpp"

namespace fg {

/// Fourier coefficient grids of N_i N_j (i, j in x, y, z), |p| <= 2Nx, |q| <= 2Nz.
struct NormalFieldCoefficients {
  std::array<std::array<fourier::CoefficientGrid, 3>, 3> NN;
  std::array<fourier::CoefficientGrid, 3> N;  ///< signed components (profile layers only)
  bool hasSigned = false;
};

/// 1D profile: N = (-g', 1, 0)/sqrt(1+g'^2) (Plan eq:N1d); exact for piecewise-linear profiles.
NormalFieldCoefficients normalField(const Profile1DLayer& layer, double dx_um, const IndexMap& map);

/// Raster of Algorithm NF on Gx x Gz points x_i = dx i/Gx, z_k = dz k/Gz (row-major, index
/// i Gz + k). NN = (N_x^2, N_x N_z, N_z^2); for the vector method also the signed (N_x, N_z).
struct NormalFieldRaster {
  int Gx = 0, Gz = 0;
  std::vector<double> nxx, nxz, nzz;
  std::vector<double> nx, nz;  ///< Vector method only
  std::vector<int> region;     ///< -1 background, k = inclusion k (periodic lookup)
};
/// Algorithm NF steps 1-4 (Plan alg:nfield): exact distances to the material boundaries of the
/// polygons (periodic images, internal edges removed), smoothing with a periodic Gaussian.
NormalFieldRaster normalFieldRaster(const Polygon2DLayer& layer, double dx_um, double dz_um, int Gx,
                                    int Gz);
/// y-invariant 2D cross-section: Algorithm NF incl. step 5 (DFT of the smooth raster).
NormalFieldCoefficients normalField(const Polygon2DLayer& layer, double dx_um, double dz_um,
                                    const IndexMap& map);

/// Index of the inclusion containing (x, z) (periodic), -1 for the background.
int polygonRegion(const Polygon2DLayer& layer, double dx_um, double dz_um, double x, double z);

/// 2D profile y = g(x, z): N = (-g_x, 1, -g_z)/sqrt(1 + g_x^2 + g_z^2) (Plan eq:N2d), sampled
/// (smooth function); SinusX reproduces the 1D sinus coefficients exactly (V9).
NormalFieldCoefficients normalField(const Profile2DLayer& layer, double dx_um, double dz_um,
                                    const IndexMap& map);

}  // namespace fg
