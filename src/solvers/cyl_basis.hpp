// FourierGrating - cylindrical field representation of the stack assembly (M20, Plan section
// "M20: S-Matrix, Loeser und Nahfeld"): one azimuthal harmonic m on the discrete Hankel spectrum.
#pragma once

#include <memory>
#include <vector>

#include "fg/cyl/excitation.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/geometry/stack.hpp"
#include "modal_basis.hpp"

namespace fg::detail {

/// Homogeneous or ring layer of a cylindrical stack (radii normalized).
class CylLayer final : public StackLayer {
 public:
  double thickness_um() const override { return thickness_; }
  bool homogeneous() const override { return homogeneous_; }
  double thickness_ = 0;
  bool homogeneous_ = true;
  cplx eps = 1.0;    ///< homogeneous layers
  cyl::Rings rings;  ///< ring layers (normalized radii)
  std::vector<cplx> epsList;
};

/// Layers of a cylindrical stack (bottom to top); rings radii times k0.
std::vector<std::unique_ptr<CylLayer>> buildCylLayers(const Stack& stack, double k0);

/// Hankel basis of one harmonic: analytic bases of superstrate and substrate (Plan eq:cyl_Psi), M
/// of the ring layers from Plan eq:cyl_M with Laurent or Li in r (precomputed, the layers are
/// y-invariant).
class CylBasis final : public ModalBasis {
 public:
  CylBasis(const cyl::SpectralGrid& grid, int m, cplx epsSup, cplx epsSub,
           const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f);
  const ModeBasis& ref() const override { return ref_; }
  const ModeBasis& sub() const override { return sub_; }
  MatC M(const StackLayer& layer, double y_um) const override;
  LayerModes modes(const MatC& M) const override;
  LayerModes homogeneousModes(const StackLayer& layer) const override;
  double rhoEstimate(const StackLayer& layer) const override;

  /// Plane-wave background with lateral wave number kpar along the azimuth Psi (Plan eq:cyl_jacobi):
  /// prepares the source spectra of the ring layers.
  void setPlaneWave(double kpar, double Psi);
  /// sigma of Plan eq:cyl_scat in a ring layer for the background (c_+, c_-, e_y) on the axis;
  /// zero vector for homogeneous layers (their material is the background).
  VecC source(const StackLayer& layer, cplx cp, cplx cm, cplx ey) const;

 private:
  std::size_t index(const StackLayer& layer) const;
  const cyl::SpectralGrid& grid_;
  ModeBasis ref_, sub_;
  std::vector<const StackLayer*> layers_;
  int m_;
  cyl::CylFactorization f_;
  std::vector<cyl::CylQ> Q_;
  std::vector<MatC> M_;
  std::vector<Eigen::PartialPivLU<MatC>> Qyy_;
  std::vector<cyl::PlaneWaveSource> src_;
};

}  // namespace fg::detail
