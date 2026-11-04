#include "fg/fourier/coefficients.hpp"

#include <cmath>

#include "fg/core/error.hpp"

namespace fg::fourier {

cplx intervalCoefficient(double x1, double x2, double period, int p) {
  FG_REQUIRE(x2 >= x1 && period > 0.0, "invalid interval");
  if (p == 0) return {(x2 - x1) / period, 0.0};
  const double K = 2.0 * PI / period;
  const double k = p * K;
  // (exp(-i k x2) - exp(-i k x1)) / (-i k d)
  return (std::exp(-I * k * x2) - std::exp(-I * k * x1)) / (-I * k * period);
}

std::vector<cplx> piecewiseConstant1D(const std::vector<double>& breakpoints,
                                      const std::vector<cplx>& values, double period, int pmax) {
  FG_REQUIRE(breakpoints.size() == values.size() + 1, "breakpoints/values size mismatch");
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1), 0.0);
  for (std::size_t k = 0; k < values.size(); ++k) {
    for (int p = -pmax; p <= pmax; ++p) {
      c[static_cast<std::size_t>(p + pmax)] +=
          values[k] * intervalCoefficient(breakpoints[k], breakpoints[k + 1], period, p);
    }
  }
  return c;
}

cplx rectangleCoefficient(double x1, double x2, double z1, double z2, double dx, double dz, int p,
                          int q) {
  return intervalCoefficient(x1, x2, dx, p) * intervalCoefficient(z1, z2, dz, q);
}

cplx polygonCoefficient(const std::vector<std::array<double, 2>>& v, double dx, double dz, int p,
                        int q) {
  FG_REQUIRE(v.size() >= 3, "polygon needs at least 3 vertices");
  const double kx = p * 2.0 * PI / dx;
  const double kz = q * 2.0 * PI / dz;
  const std::size_t L = v.size();
  if (p == 0 && q == 0) {  // area (shoelace), counter-clockwise -> positive
    double a = 0.0;
    for (std::size_t j = 0; j < L; ++j) {
      const auto& P0 = v[j];
      const auto& P1 = v[(j + 1) % L];
      a += P0[0] * P1[1] - P1[0] * P0[1];
    }
    return {0.5 * a / (dx * dz), 0.0};
  }
  const double k2 = kx * kx + kz * kz;
  cplx sum = 0.0;
  for (std::size_t j = 0; j < L; ++j) {
    const auto& P0 = v[j];
    const auto& P1 = v[(j + 1) % L];
    const double tx = P1[0] - P0[0], tz = P1[1] - P0[1];
    const double len = std::hypot(tx, tz);
    if (len == 0.0) continue;
    // outward normal of a counter-clockwise polygon in the (x, z) plane
    const double nx = tz / len, nz = -tx / len;
    const double cx = 0.5 * (P0[0] + P1[0]), cz = 0.5 * (P0[1] + P1[1]);
    const double u = 0.5 * (kx * tx + kz * tz);
    const double sinc = (std::abs(u) < 1e-12) ? 1.0 : std::sin(u) / u;
    sum += (kx * nx + kz * nz) * len * std::exp(-I * (kx * cx + kz * cz)) * sinc;
  }
  return I / k2 * sum / (dx * dz);
}

}  // namespace fg::fourier
