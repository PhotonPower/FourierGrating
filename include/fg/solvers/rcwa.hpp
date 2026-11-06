// FourierGrating - RCWA building blocks (Plan section "Loeser B").
#pragma once

#include <functional>
#include <vector>

#include "fg/physics/modes.hpp"
#include "fg/smatrix/smatrix.hpp"

namespace fg {

/// Eigenmodes of a y-invariant slice: up-going (field ~ exp(lu (y - y_bottom)), Re lu <= 0)
/// and down-going (field ~ exp(ld (y - y_top)), Re ld >= 0). Columns are 4 n_f state vectors.
struct LayerModes {
  MatC Wu, Wd;
  VecC lu, ld;
};

/// Relative tolerance tau / max|lambda| of the classification: 1e-9 in FP64, 1e-4 when the
/// eigenproblems run in single precision (Plan section "M16: Angepasste Parameter").
double classificationTolerance();
/// Order of the modes for the up/down classification (Plan eq:classify): sorted by Re lambda, modes
/// with |Re lambda| <= tau by the flux sign (flux(k) of eigenvector k); first half up-going.
std::vector<Eigen::Index> modeOrder(const VecC& lam, const std::function<double(Eigen::Index)>& flux);
/// Flux 1/2 Re(Ez Hx* - Ex Hz*) of a 4 n_f state vector; of a TE/TM block (E, H) with its sign.
double modeFlux(const VecC& w);
double blockModeFlux(const VecC& w, double fluxSign);

/// True if M has the block anti-diagonal structure of isotropic y-invariant layers.
bool isBlockAntiDiagonal(const MatC& M, double tol = 1e-13);
/// Reduced eigenproblem P R u = mu u (Plan eq:rcwa_reduced).
LayerModes layerModesReduced(const MatC& M, double eta = 1e-10);
/// Full eigen-decomposition with classification (Plan eq:classify).
LayerModes layerModesGeneral(const MatC& M);
/// Picks the cheapest exact path: reduced (anti-diagonal M), TE/TM-decoupled blocks
/// (1D in-plane: (Ex, Hz) and (Ez, Hx)), or the full eigenproblem.
LayerModes layerModesAuto(const MatC& M);
/// Modes of one polarization of the decoupled 1D in-plane problem (M from
/// buildPropagationMatrix(Q, alpha, p), 2 n_f rows): reduced path or classified eigenproblem.
LayerModes layerModesPolarization(const MatC& M, Polarization p);
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
