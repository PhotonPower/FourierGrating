// FourierGrating - Bessel functions of integer order and complex argument (M20, Plan section
// "Galerkin-Matrizen": complex arguments on the deformed spectral path).
#pragma once

#include <vector>

#include "fg/core/types.hpp"

namespace fg::cyl {

/// J_0(z) .. J_nmax(z) for complex z (Miller backward recurrence, normalized with
/// exp(-+iz) = J_0 + 2 sum_n (-+i)^n J_n). Requires |Im z| <= 700.
std::vector<cplx> besselJAll(int nmax, cplx z);

/// J_n(z) for any integer n (J_{-n} = (-1)^n J_n).
cplx besselJ(int n, cplx z);

/// J_n'(z) = (J_{n-1}(z) - J_{n+1}(z)) / 2.
cplx besselJp(int n, cplx z);

}  // namespace fg::cyl
