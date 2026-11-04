// FourierGrating - (block-)Toeplitz matrices [[f]]_{r r'} = f_{n-n', m-m'} (Plan eq:toeplitz).
#pragma once

#include <functional>
#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/core/types.hpp"
#include "fg/fourier/coefficients.hpp"

namespace fg::fourier {

/// 1D Toeplitz matrix of size (2N+1); coeffs holds p = -2N..2N at index p + 2N.
MatC toeplitz1D(const std::vector<cplx>& coeffs, int N);

/// General (block-)Toeplitz matrix for an index map; coef(dp, dq) returns f_{dp,dq}.
MatC toeplitz(const IndexMap& map, const std::function<cplx(int, int)>& coef);

/// Toeplitz matrix from a coefficient grid (requires Px >= 2 Nx, Pz >= 2 Nz).
MatC toeplitz(const IndexMap& map, const CoefficientGrid& grid);

}  // namespace fg::fourier
