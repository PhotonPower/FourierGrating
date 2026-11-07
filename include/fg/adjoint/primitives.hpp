// FourierGrating - reverse-mode (adjoint) rules of the building blocks (Plan section "M19: Konvention").
// Convention: for a real objective L and a complex quantity A, Abar = dL/dRe A + i dL/dIm A, so that
// dL = Re tr(Abar^H dA). Backward functions ADD to the cotangents they receive (null: not needed).
#pragma once

#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/core/types.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/physics/propagation_matrix.hpp"
#include "fg/physics/qmatrix.hpp"
#include "fg/smatrix/smatrix.hpp"
#include "fg/solvers/butcher.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::adjoint {

/// Adds b to a (a empty: a = b).
void accumulate(MatC& a, const MatC& b);
void accumulate(VecC& a, const VecC& b);

/// C = c P X (c scalar): Pbar += conj(c) Cbar X^H, Xbar += conj(c) P^H Cbar.
void product(cplx c, const MatC& P, const MatC& X, const MatC& Cbar, MatC* Pbar, MatC* Xbar);
/// X = A^-1 B: Bhat = A^-H Xbar, Bbar += Bhat, Abar -= Bhat X^H.
void solve(const MatC& A, const MatC& X, const MatC& Xbar, MatC* Abar, MatC* Bbar);
/// C = A^-1: Abar -= C^H Cbar C^H.
void inverse(const MatC& C, const MatC& Cbar, MatC& Abar);

/// Eigen-decomposition A = W diag(lambda) W^-1 (Plan "Eigenzerlegung"): Abar = W^-H (diag(lambdabar) +
/// conj(P) o (W^H Wbar)) W^H, P_ij = 1/(lambda_j - lambda_i). Requires L to be independent of the
/// normalization of the eigenvectors. Pairs within one class (cls[i] == cls[j], e.g. up/down) closer than
/// tau max|lambda| (tau = 1e-10) have P_ij = 0 (only their subspace enters); nearly degenerate pairs
/// get the Lorentz broadening conj(dl) / (|dl|^2 + delta^2), delta = 1e-12 max|lambda|.
MatC eig(const MatC& W, const VecC& lambda, const MatC& Wbar, const VecC& lambdabar, const std::vector<int>& cls);

/// T = toeplitz(map, f): fbar_pq += sum over the diagonals n - n' = p, m - m' = q of Tbar.
void toeplitz(const IndexMap& map, const MatC& Tbar, fourier::CoefficientGrid& fbar);

/// Redheffer star product C = A * B (Plan eq:star).
void star(const SMatrix& A, const SMatrix& B, const SMatrix& Cbar, SMatrix& Abar, SMatrix& Bbar);
/// S = fromTransfer(T) (Plan eq:T2S).
void fromTransfer(const MatC& T, const SMatrix& Sbar, MatC& Tbar);
/// Slice S-matrix S = H G^-1 of stableLayerSMatrix: cotangents of the mode matrices, eigenvalues and the
/// (normalized) slice thickness h.
void layerSMatrix(const LayerModes& m, double h, const ModeBasis& ref, const SMatrix& Sbar, MatC& Wubar, MatC& Wdbar,
                  VecC& lubar, VecC& ldbar, double& hbar, MatC* PsiInvbar = nullptr);
/// mb = isotropicModes(eps, alpha, gamma): cotangents of eps, alpha, gamma from those of Psi, PsiInv and beta
/// (empty: none).
void isotropicModes(cplx eps, const VecC& alpha, const VecC& gamma, const ModeBasis& mb, const MatC& Psibar,
                    const MatC& PsiInvbar, const VecC& betabar, cplx& epsbar, VecC& alphabar, VecC& gammabar);

/// M = buildPropagationMatrix(Q, alpha, gamma) (Plan eq:Mgeneral): Qbar from Mbar (alpha, gamma fixed).
void propagationMatrix(const QBlocks& Q, const VecC& alpha, const VecC& gamma, const MatC& Mbar, QBlocks& Qbar);
/// The same with the cotangents of alpha and gamma (M19 Phase 2: angle, wavelength, superstrate).
void propagationMatrix(const QBlocks& Q, const VecC& alpha, const VecC& gamma, const MatC& Mbar, QBlocks& Qbar,
                       VecC& alphabar, VecC& gammabar);
/// The same with dense A, G (M17).
void propagationMatrix(const QBlocks& Q, const MatC& A, const MatC& G, const MatC& Mbar, QBlocks& Qbar);
/// The same, also the cotangents of the dense A, G (transformed coordinates, M19 Phase 2); Abar, Gbar are added to
/// (empty: started at zero).
void propagationMatrix(const QBlocks& Q, const MatC& A, const MatC& G, const MatC& Mbar, QBlocks& Qbar, MatC& Abar,
                       MatC& Gbar);

/// Cotangents of the coefficient grids of an isotropic slice (Plan "Faktorisierung"); the normal field
/// is frozen (Plan "Methode", 2.).
struct IsoGridBar {
  fourier::CoefficientGrid eps, invEps;
  std::array<std::array<fourier::CoefficientGrid, 3>, 3> NN;  ///< FFF: cotangents of [[N_i N_j]]
  std::array<fourier::CoefficientGrid, 3> N;                  ///< FFF_rig: cotangents of [[N_i]]
};
void qLaurent(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g);
void qLiLamellar(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g);
void qFFFIsotropic(const SliceCoefficients& c, const IndexMap& map, bool symmetric, const QBlocks& Qbar, IsoGridBar& g);
void qFFFIsotropicRigorous(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g);
/// Anisotropic Laurent: cotangents of the nine grids eps_ij.
void qLaurentTensor(const IndexMap& map, const QBlocks& Qbar, std::array<std::array<fourier::CoefficientGrid, 3>, 3>& g);

/// Anisotropic factorizations (M19 Phase 2): Q from the grids g of the coefficient functions in the order of
/// the layer models (Laurent: eps_ij; FFF_A: C, eps C (Plan eq:QaniA); FFF_B: G1, G2, H1, P_T, 1/xi (Plan
/// eq:QaniB)); gbar receives the cotangents of all grids (same sizes as g).
void qAnisotropic(Factorization f, const std::vector<fourier::CoefficientGrid>& g, const IndexMap& map, const QBlocks& Qbar,
                  std::vector<fourier::CoefficientGrid>& gbar);
/// One RK step P = rkStepPropagator(nodes, h, tab): cotangents of the stage matrices (and of h).
void rkStep(const std::vector<MatC>& nodes, double h, const ButcherTableau& tab, const MatC& Pbar, std::vector<MatC>& nodesBar,
            double* hbar = nullptr);

}  // namespace fg::adjoint
