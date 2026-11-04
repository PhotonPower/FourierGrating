// FourierGrating - exact Fourier coefficients of piecewise constant functions.
// Convention (Plan eq:fcoef): f(x,z) = sum_pq f_pq exp(i(p Kx x + q Kz z)),
// f_pq = 1/(dx dz) * integral over the unit cell of f exp(-i(p Kx x + q Kz z)).
#pragma once

#include <array>
#include <vector>

#include "fg/core/types.hpp"

namespace fg::fourier {

/// Coefficient p of the indicator of [x1, x2] with period d (Plan eq:interval). Requires x1 <= x2.
cplx intervalCoefficient(double x1, double x2, double period, int p);

/// Coefficients p = -pmax..pmax (vector index p + pmax) of a piecewise constant 1D function.
/// The cell [0, period) is split at the sorted breakpoints x_0 = 0 < x_1 < ... < x_K = period;
/// values[k] holds on [x_k, x_{k+1}]. breakpoints.size() == values.size() + 1.
std::vector<cplx> piecewiseConstant1D(const std::vector<double>& breakpoints,
                                      const std::vector<cplx>& values, double period, int pmax);

/// Coefficient (p, q) of the indicator of the rectangle [x1,x2] x [z1,z2].
cplx rectangleCoefficient(double x1, double x2, double z1, double z2, double dx, double dz, int p,
                          int q);

/// Coefficient (p, q) of the indicator of a simple polygon (vertices counter-clockwise, (x, z)),
/// Plan eq:polygon.
cplx polygonCoefficient(const std::vector<std::array<double, 2>>& vertices, double dx, double dz,
                        int p, int q);

/// Dense grid of 2D coefficients for |p| <= Px, |q| <= Pz (Px = 2 Nx, Pz = 2 Nz for Toeplitz use).
class CoefficientGrid {
 public:
  CoefficientGrid() : CoefficientGrid(0, 0) {}
  CoefficientGrid(int Px, int Pz) : Px_(Px), Pz_(Pz), data_((2 * Px + 1) * (2 * Pz + 1), 0.0) {}
  int Px() const { return Px_; }
  int Pz() const { return Pz_; }
  cplx& operator()(int p, int q) { return data_[idx(p, q)]; }
  cplx operator()(int p, int q) const { return data_[idx(p, q)]; }

 private:
  std::size_t idx(int p, int q) const {
    return static_cast<std::size_t>((q + Pz_) + (2 * Pz_ + 1) * (p + Px_));
  }
  int Px_, Pz_;
  std::vector<cplx> data_;
};

}  // namespace fg::fourier
