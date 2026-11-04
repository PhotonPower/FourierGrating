// FourierGrating - factorized permittivity matrices Q_eps (Plan section "Faktorisierung (FFF)").
#pragma once

#include "fg/core/index_map.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/physics/propagation_matrix.hpp"

namespace fg {

/// Fourier data of one height y inside a layer. Grids must cover |p| <= 2Nx, |q| <= 2Nz.
struct SliceCoefficients {
  bool isotropic = true;
  fourier::CoefficientGrid eps;      ///< isotropic: eps(x, z)
  fourier::CoefficientGrid invEps;   ///< isotropic: 1/eps(x, z)
  std::array<std::array<fourier::CoefficientGrid, 3>, 3> epsTensor{};  ///< anisotropic: eps_ij
  const NormalFieldCoefficients* normal = nullptr;  ///< N_i N_j (FFF)
};

/// Laurent: Q_ij = [[eps_ij]] (isotropic: delta_ij [[eps]]).
QBlocks qLaurent(const SliceCoefficients& c, const IndexMap& map);
/// Li, lamellar 1D: Q_xx = [[1/eps]]^-1, Q_yy = Q_zz = [[eps]] (book IV.5/IV.6).
QBlocks qLiLamellar(const SliceCoefficients& c, const IndexMap& map);
/// Isotropic FFF: Q_ij = delta_ij [[eps]] - Delta [[N_i N_j]], Delta = [[eps]] - [[1/eps]]^-1
/// (Plan eq:Qiso, book IV.32). symmetric = true: Q_ij = delta_ij [[eps]] - (Delta [[N_iN_j]] +
/// [[N_iN_j]] Delta)/2 (Plan eq:Qsym), which keeps the +-lambda pairing of reciprocal slices.
QBlocks qFFFIsotropic(const SliceCoefficients& c, const IndexMap& map, bool symmetric = false);
/// Isotropic FFF without the commutation approximation of the book (Plan eq:Qrig):
/// Q_ij = delta_ij [[eps]] - [[eps]] [[N_i]] [[N_j]] + [[N_i]] [[1/eps]]^-1 [[N_j]].
/// Requires the signed normal components (c.normal->hasSigned).
QBlocks qFFFIsotropicRigorous(const SliceCoefficients& c, const IndexMap& map);
/// Anisotropic FFF variant A: [[eps C]] [[C]]^-1 (Plan eq:QaniA).      TODO(M5)
QBlocks qFFFAnisotropicA(const SliceCoefficients& c, const IndexMap& map);
/// Anisotropic FFF variant B (projectors, Plan eq:QaniB).              TODO(M5)
QBlocks qFFFAnisotropicB(const SliceCoefficients& c, const IndexMap& map);

}  // namespace fg
