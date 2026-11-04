// FourierGrating - factorized permittivity matrices Q_eps (Plan section "Faktorisierung (FFF)").
#pragma once

#include "fg/core/index_map.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/physics/propagation_matrix.hpp"

namespace fg {

/// Input for one height y inside a layer: Fourier coefficients of eps (or eps_ij) and N N^T.
struct SliceCoefficients {
  bool isotropic = true;
  fourier::CoefficientGrid eps{0, 0};          ///< isotropic: eps(x, z)
  fourier::CoefficientGrid invEps{0, 0};       ///< isotropic: 1/eps(x, z)
  std::array<std::array<fourier::CoefficientGrid, 3>, 3> epsTensor{};  ///< anisotropic
  const NormalFieldCoefficients* normal = nullptr;
};

/// Laurent: Q_ij = [[eps_ij]].                                         TODO(M2)
QBlocks qLaurent(const SliceCoefficients& c, const IndexMap& map);
/// Li lamellar 1D: Q_xx = [[1/eps]]^-1, Q_yy = Q_zz = [[eps]].          TODO(M2)
QBlocks qLiLamellar(const SliceCoefficients& c, const IndexMap& map);
/// Isotropic FFF: Q_ij = delta_ij [[eps]] - Delta [[N_i N_j]] (Plan eq:Qiso). TODO(M3)
QBlocks qFFFIsotropic(const SliceCoefficients& c, const IndexMap& map);
/// Anisotropic FFF variant A: [[eps C]] [[C]]^-1 (Plan eq:QaniA).      TODO(M5)
QBlocks qFFFAnisotropicA(const SliceCoefficients& c, const IndexMap& map);
/// Anisotropic FFF variant B (projectors, Plan eq:QaniB).              TODO(M5)
QBlocks qFFFAnisotropicB(const SliceCoefficients& c, const IndexMap& map);

}  // namespace fg
