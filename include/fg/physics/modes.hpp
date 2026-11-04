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

/// Time-averaged Poynting flux S_y = 1/2 Re(Ez Hx* - Ex Hz*) summed over all orders of F.
double poyntingY(const VecC& F);

}  // namespace fg
