// FourierGrating - mode bases of homogeneous media (Plan section "Homogene Medien").
#pragma once

#include "fg/core/types.hpp"

namespace fg {

/// Branch of beta = sqrt(eps - alpha^2 - gamma^2) with Im >= 0 (Re >= 0 if Im == 0),
/// regularized by + i*eta if |beta| < eta (Rayleigh anomaly), Plan eq:beta_branch.
cplx betaBranch(cplx eps, cplx alpha, cplx gamma, double eta = 1e-10);

/// F = Psi V, V = (a_x, a_z | b_x, b_z) (down | up), Psi columns are eigenvectors of M
/// with eigenvalues lambda (Plan eq:Psi_iso). Layout of F: ([Ex],[Ez],[Hx],[Hz]).
struct ModeBasis {
  MatC Psi;
  MatC PsiInv;
  VecC lambda;  ///< eigenvalue of each column (field ~ exp(lambda y))
  VecC beta;    ///< isotropic media: beta of each order (empty otherwise)
};

/// Analytic tangential-E basis of a homogeneous isotropic medium.
ModeBasis isotropicModes(cplx eps, const VecC& alpha, const VecC& gamma, double eta = 1e-10);

/// Per-order 4x4 eigen-decomposition with up/down classification (Plan eq:classify).
/// TODO(M5): implement.
ModeBasis anisotropicModes(const Mat3C& eps, const VecC& alpha, const VecC& gamma);

/// Polarization of the decoupled 1D in-plane systems (Plan eq:TE, eq:TM):
/// TE: F = ([Ez],[Hx]), amplitudes a_z | b_z; TM: F = ([Ex],[Hz]), amplitudes a_x | b_x.
enum class Polarization { TE, TM };

/// Field-block indices (into ([Ex],[Ez],[Hx],[Hz])) and amplitude-block indices (into
/// (a_x, a_z | b_x, b_z)) of one polarization, each as (first, second) block of n_f.
struct PolarizationBlocks {
  int field[2];
  int amp[2];
};
PolarizationBlocks polarizationBlocks(Polarization p);

/// 2 n_f basis of one polarization, cut out of a 4 n_f isotropic basis with gamma = 0 (the
/// basis is block diagonal there, so the blocks of PsiInv are the inverse of the blocks of Psi).
ModeBasis polarizationBasis(const ModeBasis& full, Polarization p);

/// Time-averaged Poynting flux S_y = 1/2 Re(Ez Hx* - Ex Hz*) summed over all orders of F.
double poyntingY(const VecC& F);

}  // namespace fg
