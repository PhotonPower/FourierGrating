// FourierGrating - RCWA building blocks (Plan section "Loeser B").
#pragma once

#include "fg/physics/modes.hpp"
#include "fg/smatrix/smatrix.hpp"

namespace fg {

/// Eigenmodes of a y-invariant slice: up-going (field ~ exp(lu (y - y_bottom)), Re lu <= 0)
/// and down-going (field ~ exp(ld (y - y_top)), Re ld >= 0). Columns are 4 n_f state vectors.
struct LayerModes {
  MatC Wu, Wd;
  VecC lu, ld;
};

/// True if M has the block anti-diagonal structure of isotropic y-invariant layers.
bool isBlockAntiDiagonal(const MatC& M, double tol = 1e-13);
/// Reduced eigenproblem P R u = mu u (Plan eq:rcwa_reduced).
LayerModes layerModesReduced(const MatC& M, double eta = 1e-10);
/// Full eigen-decomposition with classification (Plan eq:classify).
LayerModes layerModesGeneral(const MatC& M);
/// Analytic modes of a homogeneous isotropic layer.
LayerModes layerModesHomogeneous(const ModeBasis& mb);

/// Data kept per slice for near-field reconstruction (Plan eq:rcwa_bc).
struct RCWASliceRecord {
  LayerModes modes;
  double h = 0;
  MatC G;  ///< system matrix of eq:rcwa_S
};

/// Stable slice S-matrix S = H G^-1 in the reference basis (Plan eq:rcwa_S).
SMatrix stableLayerSMatrix(const LayerModes& m, double h, const ModeBasis& ref,
                           RCWASliceRecord* record = nullptr);

}  // namespace fg
