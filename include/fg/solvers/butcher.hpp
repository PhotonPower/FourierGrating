// FourierGrating - Butcher tableaus of implicit Runge-Kutta methods (Plan section "Loeser A").
#pragma once

#include <string>

#include "fg/core/types.hpp"

namespace fg {

struct ButcherTableau {
  std::string name;
  MatD A;
  VecD b, c;
  int order = 0;

  int stages() const { return static_cast<int>(b.size()); }
  /// Gauss-Legendre, 2 stages, order 4 (default).
  static ButcherTableau gauss2();
  /// Gauss-Legendre, 3 stages, order 6.
  static ButcherTableau gauss3();
  /// Radau IIA, 3 stages, order 5 (L-stable).
  static ButcherTableau radau3();
  /// Lookup by name: "gauss2", "gauss3", "radau3".
  static ButcherTableau byName(const std::string& name);
};

/// Step propagator P_h of an s-stage RK method for the linear ODE dF/dy = M(y) F
/// (Plan eq:rkprop). M_nodes[i] = M(y + c_i h).
MatC rkStepPropagator(const std::vector<MatC>& M_nodes, double h, const ButcherTableau& tab);

}  // namespace fg
