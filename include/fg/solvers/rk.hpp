// FourierGrating - implicit Runge-Kutta integration of dF/dy = M(y) F (Plan section "Loeser A").
#pragma once

#include <functional>
#include <vector>

#include "fg/physics/modes.hpp"
#include "fg/smatrix/smatrix.hpp"
#include "fg/solvers/butcher.hpp"

namespace fg {

struct RKStats {
  long steps = 0;
  long rejected = 0;
};

/// Accepted steps of a transfer in the order applied (start y, step h); the transfer is the product of
/// rkStepPropagator over them (adjoint, M19: the step sequence counts as fixed).
struct RKSteps {
  std::vector<double> y, h;
};

/// Field transfer matrix over a slice of thickness H with constant M: P_{H/k}^k, k doubled
/// until the step-doubling error is below tol (Plan eq:stepcrit). rho: spectral radius estimate.
MatC rkTransferConstant(const MatC& M, double H, const ButcherTableau& tab, double tol,
                        double rho, RKStats& stats, RKSteps* record = nullptr);

/// Field transfer matrix over [y0, y1] for y-dependent M with adaptive steps (step doubling).
MatC rkTransfer(const std::function<MatC(double)>& M, double y0, double y1,
                const ButcherTableau& tab, double tol, double rho, RKStats& stats, RKSteps* record = nullptr);

/// As rkTransfer, but detects the exact TE/TM decoupling of 1D in-plane problems
/// (F = (Ex, Ez, Hx, Hz): blocks (Ex, Hz) and (Ez, Hx)) and integrates both halves separately.
MatC rkTransferAuto(const std::function<MatC(double)>& M, double y0, double y1,
                    const ButcherTableau& tab, double tol, double rho, RKStats& stats);

/// Slice S-matrix from the field transfer matrix: T = Psi_ref^-1 T_F Psi_ref (Plan eq:rkslice).
SMatrix rkSliceSMatrix(const MatC& TF, const ModeBasis& ref);

}  // namespace fg
