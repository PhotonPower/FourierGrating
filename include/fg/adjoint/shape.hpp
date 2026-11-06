// FourierGrating - shape derivatives of indicator coefficients (Plan section "Formableitungen der
// Koeffizienten", eq:shape2d).
#pragma once

#include "fg/fourier/coefficients.hpp"
#include "fg/geometry/contour.hpp"

namespace fg::adjoint {

/// d chi_pq / d theta of a region whose boundary moves with V(t) = (r_plus(t) - r_minus(t)) / (2 delta)
/// (Plan eq:shape2d, valid for all k including k = 0): plus and minus are the region at theta +- delta with
/// the same segment structure and orientation as base. For control points (polygon vertices, NURBS
/// points, circle centre and radius) r is linear in theta and V exact; for weights the error is O(delta^2).
fourier::CoefficientGrid shapeDerivative(const geometry::Region& base, const geometry::Region& plus,
                                         const geometry::Region& minus, double delta, double dx, double dz, int Px,
                                         int Pz);

}  // namespace fg::adjoint
