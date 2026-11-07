// FourierGrating - internal layer models shared by the solvers and the near-field reconstruction.
// Units: x, y in micrometres (y relative to the layer bottom); Q blocks in normalized form.
#pragma once

#include <functional>

#include <memory>
#include <string>
#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/geometry/stack.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/physics/qmatrix.hpp"
#include "fg/adjoint/gradient.hpp"
#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

namespace fg::detail {

class LayerModel;

/// Adjoint data of one slice (M19): the factorization and its coefficient grids, and the grids as a
/// linear combination eps = sum_k eps(materials[k]) chi[k] (1/eps likewise) of indicator coefficients;
/// 1D interval form: chi[k] is the indicator of [breaks[k], breaks[k + 1]]. The normal field in sc is
/// frozen (Plan section "M19: Methode", 2.).
struct AdjointSlice {
  Factorization fac = Factorization::Laurent;
  SliceCoefficients sc;
  std::vector<std::string> materials;
  std::vector<fourier::CoefficientGrid> chi;
  std::vector<double> breaks;
  std::shared_ptr<const LayerModel> keep;  ///< keeps the section model of sc.normal alive (sectioned layers)
  // anisotropic slices (M19 Phase 2): grids of the coefficient functions in the order of adjoint::qAnisotropic
  bool aniso = false;
  std::vector<fourier::CoefficientGrid> grids;
  std::vector<Mat3C> tensors;  ///< per region (order of materials)
  /// grids for other region tensors (same geometry)
  std::function<std::vector<fourier::CoefficientGrid>(const std::vector<Mat3C>&)> gridsFor;
  /// 1D layers: coefficient functions of region k at x (moving jumps)
  std::function<void(std::size_t, double, cplx*)> evalAt;
  /// profiles: grids for the depth a + delta at the same crossings (only the normal field moves)
  std::function<std::vector<fourier::CoefficientGrid>(double)> gridsDepth;
};

class LayerModel {
 public:
  virtual ~LayerModel() = default;
  virtual double thickness_um() const = 0;
  virtual bool homogeneous() const { return false; }
  virtual bool yInvariant() const { return true; }
  /// Homogeneous layers: the material.
  virtual const Material& material() const;
  /// Factorized permittivity blocks at height y (Plan section "Faktorisierung").
  virtual QBlocks Q(double y_um) const = 0;
  /// Pointwise permittivity tensor at (x, y, z) (y relative to the layer bottom).
  virtual Mat3C epsAt(double x_um, double y_um, double z_um) const = 0;
  /// Unit normal (rank-1 field, sign arbitrary) at (x, y, z) (y relative to the layer bottom) used by
  /// the factorization and by the anisotropic FFF-consistent reconstruction (Plan eq:recon). Layers
  /// whose cross-section does not change with y ignore y (Plan section "M13: Hoehenabhaengige Schichtmodelle").
  virtual Vec3D normalAt(double x_um, double y_um, double z_um) const = 0;
  /// Projector field P(x, y, z) whose coefficients are NN(y) (default N N^T). Isotropic 2D cross-
  /// sections use the tensor field of Algorithm NF (not a projector, rank1Projector() = false).
  virtual Eigen::Matrix3d projectorAt(double x_um, double y_um, double z_um) const {
    const Vec3D N = normalAt(x_um, y_um, z_um);
    return N * N.transpose();
  }
  virtual bool rank1Projector() const { return true; }
  /// Fourier coefficients of P_ij = N_i N_j at height y (FFF-consistent reconstruction).
  virtual const NormalFieldCoefficients& NN(double y_um) const = 0;
  /// Adjoint data at height y (isotropic layers without coordinate transformation); false: not available.
  virtual bool adjointSlice(double /*y_um*/, AdjointSlice& /*out*/) const { return false; }
  /// Sectioned layers: the cross-section regions at height y (shape derivatives); empty otherwise.
  virtual std::vector<geometry::Region> sectionRegions(double /*y_um*/) const { return {}; }
  /// Permittivities occurring in the layer (spectral radius estimate).
  std::vector<cplx> epsList;
  /// Number of RCWA staircase slices (y-invariant layers: 1).
  int rcwaSlices = 1;
  mutable std::vector<std::string> warnings;  ///< may grow lazily (e.g. y interpolation)
};

using LayerModels = std::vector<std::unique_ptr<LayerModel>>;

/// Builds the models of all layers of a stack (bottom to top) for a given truncation.
LayerModels buildLayerModels(const Stack& stack, const IndexMap& map, int defaultSlices);
/// Model of layer k only (Simulation::setParameter rebuilds single layers, M12).
std::unique_ptr<LayerModel> buildLayerModel(const Stack& stack, std::size_t k, const IndexMap& map, int defaultSlices);

/// Thread-local settings of a solver run (backend dispatch, FFT sampling); also needed before
/// building layer models outside a solver.
void applyThreadSettings(const SolverOptions& opt);
/// Working precision of a solve: "precision" (FG_PRECISION overrides); auto = FP64 on the CPU, the
/// device policy on the GPU (Plan section "M16: Standardwahl").
linalg::Precision resolvePrecision(const SolverOptions& opt);
/// fg::gradient with prebuilt layer models of cfg (Simulation); nullptr: built here.
GradientResult gradientWithModels(const std::string& configJson, const io::RunConfig& cfg, const Wavelength& wl,
                                  const LayerModels* prebuilt, const Objective& objective,
                                  const std::vector<std::string>& params);
/// Identical to makeSolver(opt)->solve(stack, wl) with the models of buildLayerModels.
Result solveWithModels(const SolverOptions& opt, const Stack& stack, const Wavelength& wl, const LayerModels& models);

}  // namespace fg::detail
