// FourierGrating - Hankel-Galerkin matrices and factorization in r (M20, Plan sections
// "Galerkin-Matrizen" and "Faktorisierung").
#pragma once

#include <functional>
#include <vector>

#include "fg/core/types.hpp"
#include "fg/cyl/spectrum.hpp"

namespace fg::cyl {

/// Compactly supported piecewise-constant function of r: value[i] on [radii[i-1], radii[i]]
/// (radii[-1] = 0), zero for r > radii.back().
struct RadialSteps {
  std::vector<double> radii;
  std::vector<cplx> value;
};

/// [[g]]_{nu<-nu'} of a compact piecewise-constant g (Plan eq:cyl_gal, without background):
/// Lommel integrals (Plan eq:lommel) for |nu| = |nu'|, Gauss quadrature otherwise.
MatC galerkin(const SpectralGrid& grid, int nu, int nup, const RadialSteps& g);

/// [[g]]_{nu<-nu'} of a compact g on [0, breaks.back()] by composite Gauss-Legendre quadrature
/// (panels <= pi / (2 max|k|) between the breaks, 16 nodes; Plan sec:cylgal).
MatC galerkinQuadrature(const SpectralGrid& grid, int nu, int nup,
                        const std::vector<double>& breaks, const std::function<cplx(double)>& g);

/// Rings eps_i on [r_{i-1}, r_i] (i = 1..n, r_0 = 0) and background epsInf for r > r_n.
struct Rings {
  std::vector<double> radii;
  std::vector<cplx> eps;  ///< n values, inner to outer
  cplx epsInf = 1.0;
};

/// Factorization matrices of one harmonic m: d_pm = Q_pm+ e_+ + Q_pm- e_-, d_y = Q_yy e_y.
struct CylQ {
  MatC Qpp, Qpm, Qmp, Qmm, Qyy;
};

enum class CylFactorization { Laurent, LiR };

/// Laurent (Plan sec:cylfact) or Li in r with constant eps_inf (Plan eq:cyl_li, LU of frak A).
CylQ cylQ(const SpectralGrid& grid, int m, const Rings& rings, CylFactorization f);

}  // namespace fg::cyl
