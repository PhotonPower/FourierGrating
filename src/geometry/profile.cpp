#include "fg/geometry/profile.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"

namespace fg {

namespace {
double wrap(double x, double d) {
  double r = std::fmod(x, d);
  if (r < 0) r += d;
  return r;
}

// periodic piecewise-linear profile: vertex list closed by (xs[0] + d, ys[0])
void segment(const Profile1DLayer& p, double x, double d, double& x1, double& y1, double& x2,
             double& y2) {
  const auto& xs = p.xs_um;
  const auto& ys = p.ys_um;
  const std::size_t n = xs.size();
  double xr = wrap(x - xs[0], d) + xs[0];
  std::size_t i = 0;
  while (i + 1 < n && xs[i + 1] <= xr) ++i;
  x1 = xs[i];
  y1 = ys[i];
  x2 = (i + 1 < n) ? xs[i + 1] : xs[0] + d;
  y2 = (i + 1 < n) ? ys[i + 1] : ys[0];
}
}  // namespace

double profileValue(const Profile1DLayer& p, double x, double d) {
  switch (p.shape) {
    case Profile1DLayer::Shape::Sinus:
      return 0.5 * p.depth_um * (1.0 + std::cos(2 * PI / d * (x - p.x0_um)));
    case Profile1DLayer::Shape::PiecewiseLinear: {
      double x1, y1, x2, y2;
      segment(p, x, d, x1, y1, x2, y2);
      const double xr = wrap(x - p.xs_um[0], d) + p.xs_um[0];
      return y1 + (y2 - y1) * (xr - x1) / (x2 - x1);
    }
    case Profile1DLayer::Shape::Function:
      FG_REQUIRE(static_cast<bool>(p.g), "profile function missing");
      return p.g(wrap(x, d));
  }
  return 0.0;
}

double profileSlope(const Profile1DLayer& p, double x, double d) {
  switch (p.shape) {
    case Profile1DLayer::Shape::Sinus: {
      const double K = 2 * PI / d;
      return -0.5 * p.depth_um * K * std::sin(K * (x - p.x0_um));
    }
    case Profile1DLayer::Shape::PiecewiseLinear: {
      double x1, y1, x2, y2;
      segment(p, x, d, x1, y1, x2, y2);
      return (y2 - y1) / (x2 - x1);
    }
    case Profile1DLayer::Shape::Function: {
      if (p.dg) return p.dg(wrap(x, d));
      const double h = 1e-6 * d;
      return (p.g(wrap(x + h, d)) - p.g(wrap(x - h, d))) / (2 * h);
    }
  }
  return 0.0;
}

std::vector<double> profileCrossings(const Profile1DLayer& p, double y, double d) {
  std::vector<double> xs;
  switch (p.shape) {
    case Profile1DLayer::Shape::Sinus: {
      const double c = 2.0 * y / p.depth_um - 1.0;
      if (c > -1.0 && c < 1.0) {
        const double u = std::acos(c) / (2 * PI / d);
        xs = {wrap(p.x0_um + u, d), wrap(p.x0_um - u, d)};
      }
      break;
    }
    case Profile1DLayer::Shape::PiecewiseLinear: {
      const std::size_t n = p.xs_um.size();
      for (std::size_t i = 0; i < n; ++i) {
        const double x1 = p.xs_um[i], y1 = p.ys_um[i];
        const double x2 = (i + 1 < n) ? p.xs_um[i + 1] : p.xs_um[0] + d;
        const double y2 = (i + 1 < n) ? p.ys_um[i + 1] : p.ys_um[0];
        if ((y1 - y) * (y2 - y) < 0) xs.push_back(wrap(x1 + (y - y1) * (x2 - x1) / (y2 - y1), d));
      }
      break;
    }
    case Profile1DLayer::Shape::Function: {
      const int M = 4096;
      double xa = 0, fa = profileValue(p, 0, d) - y;
      for (int k = 1; k <= M; ++k) {
        const double xb = d * k / M, fb = profileValue(p, xb, d) - y;
        if (fa * fb < 0) {
          double lo = xa, hi = xb, flo = fa;
          for (int it = 0; it < 80; ++it) {
            const double mid = 0.5 * (lo + hi), fm = profileValue(p, mid, d) - y;
            if ((fm < 0) == (flo < 0)) { lo = mid; flo = fm; } else { hi = mid; }
          }
          xs.push_back(wrap(0.5 * (lo + hi), d));
        }
        xa = xb;
        fa = fb;
      }
      break;
    }
  }
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end(), [d](double a, double b) { return std::abs(a - b) < 1e-14 * d; }),
           xs.end());
  return xs;
}

void profileIntervals(const Profile1DLayer& p, double y, double d, std::vector<double>& breaks,
                      std::vector<bool>& below) {
  const auto xs = profileCrossings(p, y, d);
  breaks.clear();
  below.clear();
  breaks.push_back(0.0);
  for (double x : xs)
    if (x > 1e-14 * d && x < d * (1 - 1e-14)) breaks.push_back(x);
  breaks.push_back(d);
  for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
    const double xm = 0.5 * (breaks[k] + breaks[k + 1]);
    below.push_back(profileValue(p, xm, d) > y);
  }
}

}  // namespace fg
