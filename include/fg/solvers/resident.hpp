// FourierGrating - resident RCWA chain (Plan section "M15: Profiling und Parallelisierung", step
// 5). Slice modes, slice S-matrices and star products on DeviceMatrix handles of one backend: with
// the CUDA backend everything stays on the GPU and only the total S-matrix is downloaded. Same
// formulas as the host path (Plan eq:rcwa_reduced, eq:classify, eq:rcwa_S, eq:star); with the
// default (host) implementation of the resident operations the results agree with it to rounding.
#pragma once

#include <optional>

#include "fg/linalg/backend.hpp"
#include "fg/physics/modes.hpp"
#include "fg/smatrix/smatrix.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::resident {

struct SMatrixD {
  linalg::DeviceMatrix S11, S12, S21, S22;
};

SMatrixD upload(linalg::Backend& d, const SMatrix& S);
SMatrix download(linalg::Backend& d, const SMatrixD& S);
/// Redheffer star product (Plan eq:star), as fg::star.
SMatrixD star(linalg::Backend& d, const SMatrixD& lower, const SMatrixD& upper);

struct LayerModesD {
  linalg::DeviceMatrix Wu, Wd;
  VecC lu, ld;
};

LayerModesD upload(linalg::Backend& d, const LayerModes& m);
/// Modes of a slice with propagation matrix M: layerModesAuto (pol empty) or
/// layerModesPolarization(M, *pol), with the eigenproblem on the backend.
LayerModesD layerModes(linalg::Backend& d, const MatC& M, std::optional<Polarization> pol);
/// Stable slice S-matrix S = H G^-1 in the reference basis (Plan eq:rcwa_S), as stableLayerSMatrix.
SMatrixD layerSMatrix(linalg::Backend& d, const LayerModesD& m, double h,
                      const linalg::DeviceMatrix& refPsiInv);

}  // namespace fg::resident
