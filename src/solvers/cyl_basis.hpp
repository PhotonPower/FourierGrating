// FourierGrating - cylindrical field representation of the stack assembly (M20, Plan section
// "M20: S-Matrix, Loeser und Nahfeld"): one azimuthal harmonic m on the discrete Hankel spectrum.
#pragma once

#include <memory>
#include <vector>

#include "fg/cyl/coupled.hpp"
#include "fg/cyl/excitation.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/geometry/stack.hpp"
#include "fg/solvers/solver.hpp"
#include "modal_basis.hpp"

namespace fg::detail {

/// Homogeneous, ring or angular (polygon2d / contour2d, M21) layer of a cylindrical stack (lengths normalized).
class CylLayer final : public StackLayer {
 public:
  double thickness_um() const override { return thickness_; }
  bool homogeneous() const override { return homogeneous_; }
  double thickness_ = 0;
  bool homogeneous_ = true;
  cplx eps = 1.0;    ///< homogeneous layers
  cyl::Rings rings;  ///< ring layers (normalized radii)
  /// polygon2d / contour2d (M21): regions with angular dependence on the background (epsList.back())
  bool angular = false;
  cyl::AngularSection section;
  std::vector<cplx> epsList;  ///< every material of the layer, its background (outer material) last
};

/// Layers of a cylindrical stack (bottom to top); rings radii and cross sections times k0. polygon2d / contour2d need
/// "symmetry": {"rotation": P} (M21) and are checked for the P-fold symmetry (tolerance 1e-9 of their largest radius).
std::vector<std::unique_ptr<CylLayer>> buildCylLayers(const Stack& stack, double k0);

/// The block of harmonic b (N_k nodes) of a block-diagonal class matrix (ring layers of a class; frak A if present).
cyl::CylQ diagonalBlock(const cyl::CylQ& Q, Eigen::Index b, Eigen::Index n);

/// The stack with every rings_profile layer replaced by its staircase of rings layers (ringsStaircase; the
/// factorization entry is kept).
Stack expandRingsProfiles(const Stack& stack, int defaultSlices);

/// Hankel basis of one harmonic or of one symmetry class (M21, Plan sec:cylP): analytic bases of superstrate and
/// substrate (Plan eq:cyl_Psi with 1 (x) kappa), M of the structured layers from Plan eq:cyl_M (precomputed, the
/// layers are y-invariant): ring layers with Laurent or Li in r, block-diagonal over the harmonics of a class;
/// angular layers coupled (coupledQ, Laurent). State index of (block b, node j): b N_k + j in every component.
class CylBasis final : public ModalBasis {
 public:
  CylBasis(const cyl::SpectralGrid& grid, int m, cplx epsSup, cplx epsSub,
           const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f);
  CylBasis(const cyl::SpectralGrid& grid, const cyl::HarmonicMap& map, cplx epsSup, cplx epsSub,
           const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f);
  const ModeBasis& ref() const override { return ref_; }
  const ModeBasis& sub() const override { return sub_; }
  MatC M(const StackLayer& layer, double y_um) const override;
  LayerModes modes(const MatC& M) const override;
  LayerModes homogeneousModes(const StackLayer& layer) const override;
  double rhoEstimate(const StackLayer& layer) const override;
  /// One group per spectral node and harmonic: superstrate, substrate and homogeneous layers decouple there.
  Eigen::Index groups() const override { return k_.size(); }

  /// Plane-wave background with lateral wave number kpar along the azimuth Psi (Plan eq:cyl_jacobi):
  /// prepares the source spectra of the ring layers.
  void setPlaneWave(double kpar, double Psi);
  /// sigma of Plan eq:cyl_scat in a ring layer for the background (c_+, c_-, e_y) on the axis;
  /// zero vector for homogeneous layers (their material is the background).
  VecC source(const StackLayer& layer, cplx cp, cplx cm, cplx ey) const;

  // ---- near fields (Plan "M20: S-Matrix, Loeser und Nahfeld")
  int harmonic() const { return map_.m0; }
  const cyl::HarmonicMap& harmonics() const { return map_; }
  cyl::CylFactorization factorization() const { return f_; }
  std::size_t index(const StackLayer& layer) const;
  /// The Galerkin matrices of layer i (empty for homogeneous layers; gradients).
  const cyl::CylQ& Q(std::size_t i) const { return Q_.at(i); }
  /// The plane-wave source spectra of layer i (after setPlaneWave; empty for homogeneous layers; gradients).
  const cyl::PlaneWaveSource& planeWaveSource(std::size_t i) const { return src_.at(i); }
  /// The source spectra of layer i for the incidence (kpar, Psi) without storing them (setPlaneWave per layer; the
  /// superstrate gradient of plane waves differentiates them in kpar, Plan v2.51).
  cyl::PlaneWaveSource planeWaveSourceAt(std::size_t i, double kpar, double Psi) const;
  /// e_y = Q_yy^{-1} (d_y - p_y) with d_y = -1/2 kappa (h_+ + h_-) (Plan eq:cyl_spec, eq:cyl_scat); p_y of
  /// the plane-wave background e_y on the axis (zero for beams). Homogeneous layers: d_y / eps.
  VecC ey(std::size_t layer, const VecC& hp, const VecC& hm, cplx eyBg = 0.0) const;
  /// Spectra of G_pm = A^{-1} ((e_+; e_-) - Apm (c_+; c_-)) in a ring layer with Li in r (Plan eq:cyl_li,
  /// eq:cyl_G); the plane-wave background itself is added in real space. Empty for Laurent and angular layers.
  VecC g(std::size_t layer, const VecC& epm, cplx cp = 0.0, cplx cm = 0.0) const;

 private:
  const cyl::SpectralGrid& grid_;
  cyl::HarmonicMap map_;
  VecC k_;  ///< 1 (x) kappa over the blocks
  ModeBasis ref_, sub_;
  std::vector<const StackLayer*> layers_;
  cyl::CylFactorization f_;
  std::vector<cyl::CylQ> Q_;
  std::vector<MatC> M_;
  std::vector<Eigen::PartialPivLU<MatC>> Qyy_;
  std::vector<cyl::PlaneWaveSource> src_;
  mutable std::vector<std::unique_ptr<Eigen::PartialPivLU<MatC>>> Alu_;  ///< LU of frak A, on demand
};

/// The cylindrical forward problem as the solver sets it up (gradients, V21i): the stack with rings_profile
/// expanded, the layers, the spectral grid, the harmonics and the incident beam amplitudes per harmonic (empty
/// for plane waves). Normalized lengths (times k0).
struct CylProblem {
  Stack stack;
  cyl::SpectralGrid grid;
  std::vector<std::unique_ptr<CylLayer>> layers;
  cyl::CylFactorization fac = cyl::CylFactorization::LiR;
  cplx epsSup = 1.0, epsSub = 1.0;
  double k0 = 0, yTop = 0;
  std::vector<int> ms;
  /// The units of the computation: one harmonic {m, 1, 0} per entry of ms (M20) or the symmetry classes (M21,
  /// symmetry.rotation); aIn per unit (2 n_f, empty for plane waves).
  std::vector<cyl::HarmonicMap> classes;
  std::vector<VecC> aIn;
  bool planeWave = false, beamScattered = false;
  /// Plane wave (scattered field): lateral wave number, azimuth, normalized power -> um^2, and per layer the
  /// background exponentials (kappa, amplitude, reference height, c_+, c_-, e_y on the axis; Plan eq:cyl_scat).
  double kpar = 0, Psi = 0, toArea = 0;
  struct BgTerm {
    cplx kappa, amp, cp, cm, ey;
    double yr = 0;
  };
  std::vector<std::vector<BgTerm>> bg;
  std::vector<double> y0;  ///< normalized bottom of every layer
};
CylProblem cylProblem(const SolverOptions& opt, const Stack& stack, const Wavelength& wl);

}  // namespace fg::detail
