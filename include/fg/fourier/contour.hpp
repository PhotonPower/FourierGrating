// FourierGrating - Fourier coefficients of regions with curved boundaries (Plan section
// "M13: Fourierkoeffizienten aus Konturintegralen").
#pragma once

#include "fg/fourier/coefficients.hpp"
#include "fg/geometry/contour.hpp"
#include "fg/geometry/coordinate_transform.hpp"

namespace fg::fourier {

/// Coefficients |p| <= Px, |q| <= Pz of the (periodized) indicator of a normalized region
/// (Plan eq:contour, eq:contour0). Adaptive Gauss-Legendre per Bezier piece: n and 2n points,
/// bisection while the piece contributions differ by more than 1e-14 * area. n0 > 0 overrides
/// the initial node count n0 = max(16, ceil(Phi / 2) + 8), Phi = |k|_max * (piece length).
CoefficientGrid regionCoefficients(const geometry::Region& region, double dx, double dz, int Px, int Pz,
                                   int n0 = 0);

/// Coefficients in transformed coordinates (Plan section "M17: Permittivitaet und Normalenfeld"): the
/// boundary is mapped, X(t) = f^{-1}(x(t)), Z(t) = g^{-1}(z(t)), X' = x'/f'(X), and integrated with
/// eq:contour as above (adaptive bisection of the parameter interval). Null maps: identity on that axis.
CoefficientGrid regionCoefficients(const geometry::Region& region, double dx, double dz, int Px, int Pz,
                                   const CoordinateMap* xmap, const CoordinateMap* zmap);

}  // namespace fg::fourier
