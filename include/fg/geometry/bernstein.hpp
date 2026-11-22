// FourierGrating - polynomials in Bernstein form on [0, 1]: evaluation, products, derivatives and simple roots
// (subdivision with the variation-diminishing property, then bisection). Used for the critical heights of NURBS
// profiles (M24) and the circle-boundary crossings of the angular coefficients (M21).
#pragma once

#include <vector>

namespace fg::geometry {

/// Binomial coefficient C(n, k) as a double.
double binom(int n, int k);

/// Value at t (de Casteljau).
double bernsteinValue(std::vector<double> c, double t);

/// Coefficients of the product (degree n + m) of two Bernstein polynomials of degrees n and m.
std::vector<double> bernsteinProduct(const std::vector<double>& a, const std::vector<double>& b);

/// Coefficients of the derivative (degree n - 1; a constant gives {0}).
std::vector<double> bernsteinDerivative(const std::vector<double>& c);

/// Simple roots in (0, 1), ascending: subdivision until one sign change of the coefficients (exactly one root),
/// then bisection on the full polynomial; a root exactly at a subdivision point is kept. Multiple roots without a
/// sign change are not reported.
std::vector<double> bernsteinRoots(const std::vector<double>& c);

}  // namespace fg::geometry
