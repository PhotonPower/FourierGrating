// FourierGrating - Gauss-Legendre quadrature (contour integrals, M13).
#pragma once

#include <vector>

namespace fg::quadrature {

struct Rule {
  std::vector<double> x, w;  ///< nodes and weights on [-1, 1]
};

/// n-point Gauss-Legendre rule (Newton on the Legendre recursion), cached per thread.
const Rule& gaussLegendre(int n);

}  // namespace fg::quadrature
