// FourierGrating - propagation matrix, mode bases and flux of one azimuthal harmonic (M20, Plan
// sections "Hankel-Darstellung und Propagationsgleichung", "Homogene Medien, Fluss und
// Leistungen").
#pragma once

#include "fg/core/types.hpp"
#include "fg/cyl/galerkin.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/physics/modes.hpp"

namespace fg::cyl {

/// M of dF/dy = M F with F = (e_+, e_-, h_+, h_-), n_f = N_k (Plan eq:cyl_M); E = Q_yy^{-1} by LU.
MatC propagationMatrix(const VecC& k, const CylQ& Q);

/// Analytic basis of a homogeneous isotropic medium (Plan eq:cyl_Psi): F = Psi V,
/// V = (a_+, a_- | b_+, b_-) (down: lambda = -i beta, up: + i beta), beta by Plan eq:beta_branch.
ModeBasis homogeneousModes(cplx eps, const VecC& k, double eta = 1e-10);

/// Poynting flux P_y of one harmonic through y = const (Plan eq:cyl_flux); only nodes on the real
/// axis contribute (deformed nodes carry no flux, Plan sec:cylgrid).
double flux(const SpectralGrid& grid, const VecC& F);

}  // namespace fg::cyl
