// FourierGrating - angular coefficients of cross sections about the axis (M21, Plan section "Winkelkoeffizienten
// aus Querschnitten"): arcs of the circle |x| = r inside a region, chi_q(r) (Plan eq:cylP_chi), critical radii and
// the check for P-fold rotational symmetry.
#pragma once

#include <vector>

#include "fg/core/types.hpp"
#include "fg/geometry/contour.hpp"

namespace fg::cyl {

/// Arc [phiA, phiB] of the circle |x| = r, phi = atan2(z, x) (from x towards z, Plan sec:cylext), phiA < phiB,
/// phiB - phiA <= 2 pi; the full circle is [0, 2 pi].
struct Arc {
  double phiA = 0, phiB = 0;
};

/// Arcs of |x| = r inside a normalized region (outer counter-clockwise, holes clockwise; Region::normalize), r > 0.
/// Crossings with each rational Bezier segment C = P(t) / w(t) are the roots of |P(t)|^2 - r^2 w(t)^2 (Bernstein
/// form, degree 2p); whether the circle enters or leaves the material follows exactly from the boundary tangent
/// (material on the left). Without crossings the circle is inside iff the region contains (r, 0). Radii on which
/// the circle touches the boundary or passes a segment end are critical radii (criticalRadii) and not meant to be
/// sampled.
std::vector<Arc> arcsInside(const geometry::Region& region, double r);

/// chi_q(r) of Plan eq:cylP_chi for q = 0..qmax (harmonic qP): chi_0 = arc length / (2 pi r),
/// chi_q = 1/(2 pi) sum (e^{-iqP phiB} - e^{-iqP phiA}) / (-iqP). chi_{-q} = conj(chi_q) (real indicator).
VecC angularCoefficients(const geometry::Region& region, double r, int P, int qmax);

/// Critical radii of a region, ascending, starting with 0: radii of the segment ends and of the interior extrema of
/// |C(t)| (tangent points, roots of the Bernstein polynomial X X' W - X^2 W' + Z Z' W - Z^2 W' of degree 3p - 1 in
/// the homogeneous coordinates). chi_q is smooth between them (Plan sec:cylP, quadrature panels in r).
std::vector<double> criticalRadii(const geometry::Region& region);

/// Is the set of regions invariant under the rotation by 2 pi / P about the axis (Plan sec:cylP)? Boundary samples
/// of every region (17 per segment), rotated, must lie on the boundary of a region of the same material within tol
/// (micrometres; Plan: 1e-9 R_s).
bool hasRotationSymmetry(const std::vector<geometry::Region>& regions, int P, double tol);

/// Largest P in [1, maxP] with hasRotationSymmetry (the order the CAD check reports, Plan "CAD und Konfiguration").
int rotationOrder(const std::vector<geometry::Region>& regions, int maxP, double tol);

}  // namespace fg::cyl
