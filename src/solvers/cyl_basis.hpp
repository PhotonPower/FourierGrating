// FourierGrating - cylindrical field representation of the stack assembly (M20, Plan section
// "M20: S-Matrix, Loeser und Nahfeld"): one azimuthal harmonic m on the discrete Hankel spectrum.
#pragma once

#include <memory>
#include <vector>

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

 private:
  const cyl::SpectralGrid& grid_;
  ModeBasis ref_, sub_;
  std::vector<const StackLayer*> layers_;
  std::vector<MatC> M_;
};

}  // namespace fg::detail
