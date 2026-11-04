// FourierGrating - Fourier coefficients of continuous functions from samples (Plan section
// "Fourierkoeffizienten", "Allgemeine Funktionen"). Only for continuous functions (normal
// vector components etc.); discontinuous permittivities use the exact formulas.
#pragma once

#include <functional>

#include "fg/fourier/coefficients.hpp"

namespace fg::fourier {

/// Coefficients p = -pmax..pmax of a d-periodic function sampled at M equidistant points.
/// Direct DFT, O(M * pmax); requires M > 2 * pmax.
std::vector<cplx> sampledCoefficients1D(const std::function<cplx(double)>& f, double period,
                                        int pmax, int samples);

/// Coefficients of K piecewise-smooth functions on the intervals [breaks_k, breaks_k+1] of one
/// period (breaks from 0 to period). f(x, k, out) writes the K values at x in interval k.
/// Composite 8-point Gauss-Legendre quadrature with |p K| h <= 1 per sub-interval (exact for
/// constants up to rounding; error ~1e-15 for smooth integrands). Result[i][p + pmax].
std::vector<std::vector<cplx>> piecewiseSmoothCoefficients(
    const std::vector<double>& breaks, int K,
    const std::function<void(double, int, cplx*)>& f, double period, int pmax);

/// 2D coefficient grid |p| <= Px, |q| <= Pz from Mx x Mz samples of f(x, z).
CoefficientGrid sampledCoefficients2D(const std::function<cplx(double, double)>& f, double dx,
                                      double dz, int Px, int Pz, int Mx, int Mz);

/// 2D coefficient grid |p| <= Px, |q| <= Pz from samples f[i Mz + k] = f(dx i/Mx, dz k/Mz)
/// (separable DFT, Plan alg:nfield step 5). Requires Mx > 2 Px, Mz > 2 Pz.
CoefficientGrid gridCoefficients2D(const std::vector<double>& f, int Mx, int Mz, int Px, int Pz);

}  // namespace fg::fourier
