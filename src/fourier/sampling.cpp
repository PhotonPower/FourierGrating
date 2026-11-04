#include "fg/fourier/sampling.hpp"

#include <cmath>

#include "fg/core/error.hpp"

namespace fg::fourier {

std::vector<cplx> sampledCoefficients1D(const std::function<cplx(double)>& f, double period,
                                        int pmax, int M) {
  FG_REQUIRE(M > 2 * pmax, "too few samples for the requested Fourier orders");
  std::vector<cplx> vals(static_cast<std::size_t>(M));
  for (int k = 0; k < M; ++k) vals[static_cast<std::size_t>(k)] = f(period * k / M);
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1));
  for (int p = -pmax; p <= pmax; ++p) {
    cplx s = 0.0;
    const cplx w = std::exp(-I * (2.0 * PI * p / M));
    cplx wk = 1.0;
    for (int k = 0; k < M; ++k) {
      s += vals[static_cast<std::size_t>(k)] * wk;
      wk *= w;
    }
    c[static_cast<std::size_t>(p + pmax)] = s / static_cast<double>(M);
  }
  return c;
}

CoefficientGrid sampledCoefficients2D(const std::function<cplx(double, double)>& f, double dx,
                                      double dz, int Px, int Pz, int Mx, int Mz) {
  FG_REQUIRE(Mx > 2 * Px && Mz > 2 * Pz, "too few samples for the requested Fourier orders");
  // separable two-pass DFT: first along z for every x sample, then along x
  std::vector<std::vector<cplx>> rows(static_cast<std::size_t>(Mx));
  for (int i = 0; i < Mx; ++i) {
    const double x = dx * i / Mx;
    rows[static_cast<std::size_t>(i)] =
        sampledCoefficients1D([&](double z) { return f(x, z); }, dz, Pz, Mz);
  }
  CoefficientGrid g(Px, Pz);
  for (int q = -Pz; q <= Pz; ++q) {
    const auto col = sampledCoefficients1D(
        [&](double x) {
          const int i = static_cast<int>(std::lround(x / dx * Mx)) % Mx;
          return rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(q + Pz)];
        },
        dx, Px, Mx);
    for (int p = -Px; p <= Px; ++p) g(p, q) = col[static_cast<std::size_t>(p + Px)];
  }
  return g;
}

}  // namespace fg::fourier
