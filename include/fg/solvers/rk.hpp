// FourierGrating - implicit Runge-Kutta integration of dF/dy = M(y) F (Plan section "Loeser A").
#pragma once

#include <functional>

#include "fg/physics/modes.hpp"
#include "fg/smatrix/smatrix.hpp"
#include "fg/solvers/butcher.hpp"

namespace fg {

struct RKStats {
  long steps = 0;
  long rejected = 0;
};

/// Field transfer matrix over a slice of thickness H with constant M: P_{H/k}^k, k doubled
/// until the step-doubling error is below tol (Plan eq:stepcrit). rho: spectral radius estimate.
MatC rkTransferConstant(const MatC& M, double H, const ButcherTableau& tab, double tol,
                        double rho, RKStats& stats);

/// Field transfer matrix over [y0, y1] for y-dependent M with adaptive steps (step doubling).
MatC rkTransfer(const std::function<MatC(double)>& M, double y0, double y1,
                const ButcherTableau& tab, double tol, double rho, RKStats& stats);

/// Slice S-matrix from the field transfer matrix: T = Psi_ref^-1 T_F Psi_ref (Plan eq:rkslice).
SMatrix rkSliceSMatrix(const MatC& TF, const ModeBasis& ref);

}  // namespace fg
