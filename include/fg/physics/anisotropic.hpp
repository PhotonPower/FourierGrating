// FourierGrating - pointwise matrices of the anisotropic FFF (Plan eq:Aeps, eq:Ceps, eq:projector).
#pragma once

#include "fg/core/types.hpp"

namespace fg {

/// Variant A: C_eps = A_eps^-1 with A_eps = (T1; N^T eps; T2) (rows), Plan eq:Ceps.
/// Columns: (v x T2, N, -v x T1) / xi, v = N^T eps, xi = N^T eps N. Requires T2 x T1 = N.
Mat3C cMatrix(const Mat3C& eps, const Vec3D& N, const Vec3D& T1, const Vec3D& T2);

/// Variant B (Plan eq:projector): with P_N = N N^T, P_T = 1 - P_N, xi = N^T eps N:
/// G1 = eps (1 - P_N eps / xi), G2 = eps / xi, H1 = 1 - P_N eps / xi.
struct ProjectorSet {
  Mat3C G1, G2, H1, PT;
  cplx invXi;
};
ProjectorSet projectorMatrices(const Mat3C& eps, const Vec3D& N);

/// Tangent vectors for a 1D profile normal N = (Nx, Ny, 0): T1 = (Ny, -Nx, 0), T2 = z (Plan eq:T_profile).
void profileTangents(const Vec3D& N, Vec3D& T1, Vec3D& T2);

/// Rotation of a principal tensor diag(e1, e2, e3) by z-x-z Euler angles (radians): R D R^T.
Mat3C rotatedTensor(cplx e1, cplx e2, cplx e3, double phi, double theta, double psi);

}  // namespace fg
