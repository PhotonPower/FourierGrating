#include "fg/geometry/normal_field.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/geometry/profile.hpp"

namespace fg {

NormalFieldCoefficients normalField(const Profile1DLayer& layer, double dx_um, const IndexMap& map) {
  // N = (-g', 1, 0) / sqrt(1 + g'^2)  (Plan eq:N1d), only N_i N_j is needed
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  NormalFieldCoefficients nf;
  for (auto& row : nf.NN)
    for (auto& g : row) g = fourier::CoefficientGrid(Px, Pz);
  auto nn = [](double s, int i, int j) {
    const double w = 1.0 / (1.0 + s * s);
    const double N[3] = {-s, 1.0, 0.0};
    return N[i] * N[j] * w;
  };
  for (auto& g : nf.N) g = fourier::CoefficientGrid(Px, Pz);
  nf.hasSigned = true;
  for (int i = 0; i < 2; ++i) {
    auto ni = [i](double s) { const double w = 1.0 / std::sqrt(1.0 + s * s); return i == 0 ? -s * w : w; };
    std::vector<cplx> c;
    if (layer.shape == Profile1DLayer::Shape::PiecewiseLinear) {
      std::vector<double> br{0.0};
      std::vector<cplx> vals;
      std::vector<double> xs;
      for (double x : layer.xs_um) xs.push_back(std::fmod(std::fmod(x, dx_um) + dx_um, dx_um));
      std::sort(xs.begin(), xs.end());
      for (double x : xs)
        if (x > 0 && x < dx_um) br.push_back(x);
      br.push_back(dx_um);
      for (std::size_t k = 0; k + 1 < br.size(); ++k)
        vals.push_back(ni(profileSlope(layer, 0.5 * (br[k] + br[k + 1]), dx_um)));
      c = fourier::piecewiseConstant1D(br, vals, dx_um, Px);
    } else {
      const int M = std::max(1024, 16 * (Px + 1));
      c = fourier::sampledCoefficients1D([&](double x) { return cplx(ni(profileSlope(layer, x, dx_um)), 0.0); },
                                         dx_um, Px, M);
    }
    for (int p = -Px; p <= Px; ++p) nf.N[static_cast<std::size_t>(i)](p, 0) = c[static_cast<std::size_t>(p + Px)];
  }
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      std::vector<cplx> c;
      if (layer.shape == Profile1DLayer::Shape::PiecewiseLinear) {
        // slope is piecewise constant -> exact interval coefficients
        std::vector<double> br{0.0};
        std::vector<cplx> vals;
        std::vector<double> xs;
        for (double x : layer.xs_um) xs.push_back(std::fmod(std::fmod(x, dx_um) + dx_um, dx_um));
        std::sort(xs.begin(), xs.end());
        for (double x : xs)
          if (x > 0 && x < dx_um) br.push_back(x);
        br.push_back(dx_um);
        for (std::size_t k = 0; k + 1 < br.size(); ++k)
          vals.push_back(nn(profileSlope(layer, 0.5 * (br[k] + br[k + 1]), dx_um), i, j));
        c = fourier::piecewiseConstant1D(br, vals, dx_um, Px);
      } else {
        const int M = std::max(1024, 16 * (Px + 1));
        c = fourier::sampledCoefficients1D(
            [&](double x) { return cplx(nn(profileSlope(layer, x, dx_um), i, j), 0.0); }, dx_um, Px, M);
      }
      for (int p = -Px; p <= Px; ++p)
        nf.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)](p, 0) = c[static_cast<std::size_t>(p + Px)];
    }
  return nf;
}

NormalFieldCoefficients normalField(const Polygon2DLayer&, double, double, const IndexMap&) {
  throw NotImplemented("normalField(Polygon2DLayer) - milestone M7, Plan algorithm NF");
}

}  // namespace fg
