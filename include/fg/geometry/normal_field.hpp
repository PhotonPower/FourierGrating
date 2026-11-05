// FourierGrating - normal vector fields for FFF (Plan section "Normalen- und Tangentenfelder").
#pragma once

#include <functional>
#include <string>
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
/// i Gz + k). NN = (N_x^2, N_x N_z, N_z^2) according to the method.
/// (nx, nz): unit rank-1 field for the projector variant B (M8, sign arbitrary): normalized
/// smoothed signed-distance gradient, principal axis of the smoothed tensor where it vanishes.
struct NormalFieldRaster {
  int Gx = 0, Gz = 0;
  std::vector<double> nxx, nxz, nzz;
  std::vector<double> nx, nz;
  std::vector<int> region;     ///< -1 background, k = inclusion k (periodic lookup)
  /// 3D normals (boundary samples with n_y != 0, M13): the y components; empty otherwise
  std::vector<double> nxy, nyy, nyz, ny;
  bool has3D() const { return !nyy.empty(); }
};

/// Boundary sample of Algorithm NF (Plan section "M13: Normalenfeld mit Randstichproben"): point
/// (x, z) on a material boundary and unit normal (3D; y-invariant layers: n_y = 0), oriented out of
/// the region it belongs to. exactChord: the chord to the next sample is the boundary itself
/// (polygon edge, normal of this sample along the whole chord); otherwise it approximates a curve
/// with sagitta <= 1e-3 raster spacings (normals interpolated between the two samples).
struct BoundarySample {
  double x = 0, z = 0;
  Vec3D n = Vec3D(1, 0, 0);
  bool exactChord = false;
};
/// Consecutive samples of one boundary curve (closed: chord from the last to the first sample).
struct BoundaryChain {
  std::vector<BoundarySample> samples;
  bool closed = true;
};

/// Algorithm NF steps 1-4 from boundary samples (M13): nearest chord per raster point via a bucket
/// grid (periodic images), foot point on the chord; in-plane direction as in Algorithm NF, n_y
/// interpolated linearly along the chord. Only chords with different materials on both sides count
/// (regionAt: index of the region at (x, z), -1 background; materialOf(k): its material name).
/// For polygons (one exact chord per edge) the result is bit-identical to normalFieldRaster(Polygon2D).
NormalFieldRaster normalFieldRaster(const std::vector<BoundaryChain>& chains,
                                    const std::function<int(double, double)>& regionAt,
                                    const std::function<std::string(int)>& materialOf, const NormalFieldOptions& opt,
                                    double dx_um, double dz_um, int Gx, int Gz);
/// One exact chord per polygon edge, normal (t_z, 0, -t_x)/|t| (counter-clockwise: outward).
std::vector<BoundaryChain> polygonChains(const std::vector<std::array<double, 2>>& vertices);
/// Step 5 (DFT of the smooth raster): coefficients of the projector field, |p| <= 2Nx, |q| <= 2Nz.
NormalFieldCoefficients rasterCoefficients(const NormalFieldRaster& R, const IndexMap& map);
/// Raster size used for N Fourier orders: max(grid, 8 (2N + 1)), rounded up to a power of two.
int normalFieldGrid(const NormalFieldOptions& o, int N);
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

/// Gradient (dg/dx, dg/dz) and unit normal (Plan eq:N2d) of a 2D profile.
std::array<double, 2> profile2DGradient(const Profile2DLayer& l, double dx_um, double dz_um, double x, double z);
Vec3D profile2DNormal(const Profile2DLayer& l, double dx_um, double dz_um, double x, double z);

}  // namespace fg
