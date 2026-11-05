// FourierGrating - Fourier coefficients of continuous functions from samples (Plan section
// "Fourierkoeffizienten", "Allgemeine Funktionen"). Only for continuous functions (normal
// vector components etc.); discontinuous permittivities use the exact formulas.
#pragma once

#include <functional>

#include "fg/fourier/coefficients.hpp"

namespace fg::fourier {

/// Forward DFT X_p = sum_k a_k exp(-2 pi i p k / n) of any length n in O(n log n): radix 2 for
/// powers of two, Bluestein otherwise (Plan "Allgemeine Funktionen", v2.2).
std::vector<cplx> fft(std::vector<cplx> a);

/// Thread-local switch (SolverOptions::fftSampling, JSON "fft_sampling"): true = the sampling
/// functions below use fft() for every sample count; false (default) = v0.10 paths (direct DFT,
/// radix 2 only in gridCoefficients2D for powers of two), bit-identical to v0.10.
void setFFTSampling(bool on);
bool fftSampling();

/// Coefficients p = -pmax..pmax of a d-periodic function sampled at M equidistant points;
/// requires M > 2 * pmax.
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
/// Complex samples; separable (radix-2 FFT for powers of two or with fftSampling(), else DFT).
CoefficientGrid gridCoefficients2D(const std::vector<cplx>& f, int Mx, int Mz, int Px, int Pz);

/// Coefficients of a product of a smooth function (coefficients a, |p'| <= a.Px()) and a
/// discontinuous function with exact coefficients b: c_pq = sum a_p'q' b_(p-p')(q-q'), |p| <= Px,
/// |q| <= Pz. b must cover |p| <= Px + a.Px(), |q| <= Pz + a.Pz() (Plan section FFF_B, 2D).
CoefficientGrid convolveCoefficients(const CoefficientGrid& a, const CoefficientGrid& b, int Px, int Pz);

}  // namespace fg::fourier
