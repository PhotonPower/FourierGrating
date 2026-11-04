#include "fg/fourier/toeplitz.hpp"

#include "fg/core/error.hpp"

namespace fg::fourier {

MatC toeplitz1D(const std::vector<cplx>& coeffs, int N) {
  FG_REQUIRE(static_cast<int>(coeffs.size()) == 4 * N + 1, "need coefficients p = -2N..2N");
  const int n = 2 * N + 1;
  MatC T(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) T(i, j) = coeffs[static_cast<std::size_t>(i - j + 2 * N)];
  return T;
}

MatC toeplitz(const IndexMap& map, const std::function<cplx(int, int)>& coef) {
  const int nf = map.size();
  MatC T(nf, nf);
  for (int r = 0; r < nf; ++r) {
    const auto [n1, m1] = map.orders(r);
    for (int s = 0; s < nf; ++s) {
      const auto [n2, m2] = map.orders(s);
      T(r, s) = coef(n1 - n2, m1 - m2);
    }
  }
  return T;
}

MatC toeplitz(const IndexMap& map, const CoefficientGrid& grid) {
  FG_REQUIRE(grid.Px() >= 2 * map.Nx() && grid.Pz() >= 2 * map.Nz(),
             "coefficient grid too small for Toeplitz matrix");
  return toeplitz(map, [&grid](int p, int q) { return grid(p, q); });
}

}  // namespace fg::fourier
