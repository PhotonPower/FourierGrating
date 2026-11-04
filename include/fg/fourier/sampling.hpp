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

/// 2D coefficient grid |p| <= Px, |q| <= Pz from Mx x Mz samples of f(x, z).
CoefficientGrid sampledCoefficients2D(const std::function<cplx(double, double)>& f, double dx,
                                      double dz, int Px, int Pz, int Mx, int Mz);

}  // namespace fg::fourier
