// FourierGrating - aperiodic (isolated) structures: complex coordinate stretching (PML) in lateral
// zones of the computational cell (Plan section "M18: PML und isolierte Strukturen", eq:pml).
#pragma once

#include <vector>

#include "fg/core/types.hpp"

namespace fg {

struct Stack;

/// PML zones of one axis (JSON "aperiodic": {"pml": {"x": {"width_um", "beta"}}}).
struct PmlAxisSpec {
  double width_um = 0;  ///< zone width w at both cell ends (0: no PML on this axis)
  double beta = 1.0;    ///< absorption strength (Plan eq:pml)
  bool active() const { return width_um > 0; }
};

struct AperiodicSpec {
  bool enabled = false;
  PmlAxisSpec x, z;
};

/// 1/f' of the PML on one axis (Plan eq:pml): the physical region [w, d - w] has f(X) = X; in the
/// zones, with the depth xi in [0, 1] and u = sin^2(pi xi / 2), 1/f' = (1 - u) / (1 + i beta u).
/// 1/f' is continuously differentiable and tends to 0 at the cell boundary.
class PmlMap {
 public:
  PmlMap(double period, double width, double beta);

  double period() const { return d_; }
  double width() const { return w_; }
  double beta() const { return beta_; }
  /// Physical region [w, d - w] (f(X) = X).
  double physicalBegin() const { return w_; }
  double physicalEnd() const { return d_ - w_; }

  cplx inverseDerivative(double X) const;
  /// [[1/f']]_p, p = -pmax..pmax (index p + pmax), composite Gauss-Legendre between 0, w, d - w, d.
  std::vector<cplx> inverseDerivativeCoefficients(int pmax) const;

 private:
  double d_, w_, beta_;
};

/// Checks the condition of Plan section "Voraussetzungen und Einfall": inside the PML zones every layer
/// is homogeneous along the axis (no material jump, polygon vertex, contour control point or profile
/// variation there). Throws fg::Error naming the layer otherwise.
void checkPmlZones(const Stack& stack);

/// Background stack (Plan section "Ausgaben"): every structured layer replaced by the material it has
/// in the PML zones (lamellae: the material at x = 0; polygons, contours, CAD: the background; a 1D
/// profile: two homogeneous layers split at its height in the zones). Same transformation settings.
Stack backgroundStack(const Stack& stack);

}  // namespace fg
