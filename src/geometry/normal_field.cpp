#include "fg/geometry/normal_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

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

namespace {

using Poly = std::vector<std::array<double, 2>>;

double signedArea(const Poly& v) {
  double a = 0.0;
  for (std::size_t j = 0; j < v.size(); ++j) {
    const auto& p0 = v[j];
    const auto& p1 = v[(j + 1) % v.size()];
    a += p0[0] * p1[1] - p1[0] * p0[1];
  }
  return 0.5 * a;
}

Poly counterClockwise(Poly v) {
  if (signedArea(v) < 0) std::reverse(v.begin(), v.end());
  return v;
}

/// Crossing-number test (boundary points: deterministic, either side).
bool insidePolygon(const Poly& v, double x, double z) {
  bool in = false;
  for (std::size_t j = 0, k = v.size() - 1; j < v.size(); k = j++) {
    const double xj = v[j][0], zj = v[j][1], xk = v[k][0], zk = v[k][1];
    if ((zj > z) != (zk > z) && x < (xk - xj) * (z - zj) / (zk - zj) + xj) in = !in;
  }
  return in;
}

/// Boundary segment with the outward unit normal of the inclusion it belongs to.
struct Seg {
  double ax, az, bx, bz, nx, nz;
};

const std::string& regionMaterial(const Polygon2DLayer& l, int k) {
  return k < 0 ? l.background : l.inclusions[static_cast<std::size_t>(k)].material;
}

/// Material boundaries: polygon edges (sub-pieces) whose two sides carry different materials
/// (periodic lookup). Edges on the cell boundary shared with a periodic image, edges between
/// touching inclusions of the same material and inclusions of the background material vanish.
std::vector<Seg> boundarySegments(const Polygon2DLayer& l, double dx, double dz) {
  const double off = 1e-9 * std::min(dx, dz);
  constexpr int S = 64;
  std::vector<Seg> out;
  for (const auto& inc : l.inclusions) {
    const Poly v = counterClockwise(inc.vertices);
    for (std::size_t j = 0; j < v.size(); ++j) {
      const auto& p0 = v[j];
      const auto& p1 = v[(j + 1) % v.size()];
      const double tx = p1[0] - p0[0], tz = p1[1] - p0[1], len = std::hypot(tx, tz);
      if (len == 0.0) continue;
      const double nx = tz / len, nz = -tx / len;  // outward for counter-clockwise (x, z)
      std::array<bool, S> isB{};
      for (int s = 0; s < S; ++s) {
        const double t = (s + 0.5) / S, mx = p0[0] + t * tx, mz = p0[1] + t * tz;
        const int rin = polygonRegion(l, dx, dz, mx - off * nx, mz - off * nz);
        const int rout = polygonRegion(l, dx, dz, mx + off * nx, mz + off * nz);
        isB[static_cast<std::size_t>(s)] = regionMaterial(l, rin) != regionMaterial(l, rout);
      }
      for (int s = 0; s < S;) {
        if (!isB[static_cast<std::size_t>(s)]) { ++s; continue; }
        int e = s;
        while (e < S && isB[static_cast<std::size_t>(e)]) ++e;
        Seg g{p0[0] + tx * s / S, p0[1] + tz * s / S, p0[0] + tx * e / S, p0[1] + tz * e / S, nx, nz};
        // translate the segment so that its midpoint lies in the cell
        const double sx = std::floor(0.5 * (g.ax + g.bx) / dx) * dx, sz = std::floor(0.5 * (g.az + g.bz) / dz) * dz;
        g.ax -= sx; g.bx -= sx; g.az -= sz; g.bz -= sz;
        out.push_back(g);
        s = e;
      }
    }
  }
  return out;
}

/// Periodic Gaussian smoothing of a Gx x Gz raster (row-major i Gz + k), widths in points.
void gaussSmooth(std::vector<double>& f, int Gx, int Gz, double sx, double sz) {
  auto kernel = [](double s) {
    std::vector<double> w;
    if (s <= 0) return std::vector<double>{1.0};
    const int R = static_cast<int>(std::ceil(4.0 * s));
    double sum = 0;
    for (int k = -R; k <= R; ++k) {
      w.push_back(std::exp(-0.5 * k * k / (s * s)));
      sum += w.back();
    }
    for (double& x : w) x /= sum;
    return w;
  };
  const auto wx = kernel(sx), wz = kernel(sz);
  const int Rx = static_cast<int>(wx.size() / 2), Rz = static_cast<int>(wz.size() / 2);
  std::vector<double> tmp(f.size(), 0.0);
  auto at = [Gz](int i, int k) { return static_cast<std::size_t>(i) * static_cast<std::size_t>(Gz) + static_cast<std::size_t>(k); };
  for (int i = 0; i < Gx; ++i)
    for (int k = 0; k < Gz; ++k) {
      double s = 0;
      for (int r = -Rz; r <= Rz; ++r) s += wz[static_cast<std::size_t>(r + Rz)] * f[at(i, ((k + r) % Gz + Gz) % Gz)];
      tmp[at(i, k)] = s;
    }
  for (int i = 0; i < Gx; ++i)
    for (int k = 0; k < Gz; ++k) {
      double s = 0;
      for (int r = -Rx; r <= Rx; ++r) s += wx[static_cast<std::size_t>(r + Rx)] * tmp[at(((i + r) % Gx + Gx) % Gx, k)];
      f[at(i, k)] = s;
    }
}

}  // namespace

int polygonRegion(const Polygon2DLayer& l, double dx, double dz, double x, double z) {
  for (std::size_t k = 0; k < l.inclusions.size(); ++k) {
    const Poly& v = l.inclusions[k].vertices;
    double x0 = v[0][0], x1 = x0, z0 = v[0][1], z1 = z0;
    for (const auto& p : v) {
      x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]);
      z0 = std::min(z0, p[1]); z1 = std::max(z1, p[1]);
    }
    // periodic images (x - i dx, z - j dz) that can fall into the bounding box
    for (int i = static_cast<int>(std::ceil((x - x1) / dx)); i <= static_cast<int>(std::floor((x - x0) / dx)); ++i)
      for (int j = static_cast<int>(std::ceil((z - z1) / dz)); j <= static_cast<int>(std::floor((z - z0) / dz)); ++j)
        if (insidePolygon(v, x - i * dx, z - j * dz)) return static_cast<int>(k);
  }
  return -1;
}

NormalFieldRaster normalFieldRaster(const Polygon2DLayer& l, double dx, double dz, int Gx, int Gz) {
  // Plan alg:nfield, steps 1-4 (method: see NormalFieldOptions)
  FG_REQUIRE(Gx > 0 && Gz > 0 && dx > 0 && dz > 0, "normal field raster: invalid size");
  const bool vec = l.normalField.method == NormalFieldOptions::Method::Vector;
  const std::vector<Seg> segs = boundarySegments(l, dx, dz);
  const std::size_t n = static_cast<std::size_t>(Gx) * static_cast<std::size_t>(Gz);
  NormalFieldRaster R;
  R.Gx = Gx;
  R.Gz = Gz;
  R.nxx.assign(n, 0.5);
  R.nxz.assign(n, 0.0);
  R.nzz.assign(n, 0.5);
  R.region.assign(n, -1);
  std::vector<double> vx(n, 0.0), vz(n, 0.0);
  const double L = std::max(dx, dz), tie = 1e-12 * L;
  for (int i = 0; i < Gx; ++i)
    for (int k = 0; k < Gz; ++k) {
      const std::size_t idx = static_cast<std::size_t>(i) * static_cast<std::size_t>(Gz) + static_cast<std::size_t>(k);
      const double x = dx * i / Gx, z = dz * k / Gz;
      const int reg = polygonRegion(l, dx, dz, x, z);
      R.region[idx] = reg;
      if (segs.empty()) continue;  // homogeneous cross-section: N irrelevant
      const double sgn = reg < 0 ? 1.0 : -1.0;  // signed distance: > 0 in the background
      double best = std::numeric_limits<double>::infinity();
      double axx = 0, axz = 0, azz = 0, ax = 0, az = 0;
      int cnt = 0;
      for (const Seg& g : segs)
        for (int di = -1; di <= 1; ++di)
          for (int dj = -1; dj <= 1; ++dj) {
            const double ax0 = g.ax + di * dx, az0 = g.az + dj * dz, bx = g.bx + di * dx, bz = g.bz + dj * dz;
            const double tx = bx - ax0, tz = bz - az0;
            const double t = std::clamp(((x - ax0) * tx + (z - az0) * tz) / (tx * tx + tz * tz), 0.0, 1.0);
            const double rx = x - (ax0 + t * tx), rz = z - (az0 + t * tz);
            const double d = std::hypot(rx, rz);
            if (d > best + tie) continue;
            double ux, uz;
            if (d > tie) {
              ux = sgn * rx / d;
              uz = sgn * rz / d;
            } else {
              ux = g.nx;
              uz = g.nz;
            }
            if (d < best - tie) {
              best = d;
              axx = axz = azz = ax = az = 0;
              cnt = 0;
            }
            axx += ux * ux;
            axz += ux * uz;
            azz += uz * uz;
            ax += ux;
            az += uz;
            ++cnt;
          }
      // ties (medial axis): average of the candidates' projectors / vectors
      R.nxx[idx] = axx / cnt;
      R.nxz[idx] = axz / cnt;
      R.nzz[idx] = azz / cnt;
      vx[idx] = ax / cnt;
      vz[idx] = az / cnt;
    }
  if (segs.empty()) return R;
  // step 3: periodic Gaussian with the same physical width along x and z
  const double hphys = l.normalField.sigma * std::min(dx / Gx, dz / Gz);
  const double sx = hphys * Gx / dx, sz = hphys * Gz / dz;
  if (!vec) {
    for (auto* f : {&R.nxx, &R.nxz, &R.nzz}) gaussSmooth(*f, Gx, Gz, sx, sz);
    return R;
  }
  gaussSmooth(vx, Gx, Gz, sx, sz);
  gaussSmooth(vz, Gx, Gz, sx, sz);
  double gmax = 0;
  for (std::size_t j = 0; j < n; ++j) gmax = std::max(gmax, std::hypot(vx[j], vz[j]));
  R.nx.assign(n, 0.0);
  R.nz.assign(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double gn = std::hypot(vx[j], vz[j]);
    if (gn > l.normalField.tau * gmax) {  // step 4
      R.nx[j] = vx[j] / gn;
      R.nz[j] = vz[j] / gn;
      R.nxx[j] = R.nx[j] * R.nx[j];
      R.nxz[j] = R.nx[j] * R.nz[j];
      R.nzz[j] = R.nz[j] * R.nz[j];
    } else {
      R.nxx[j] = 0.5;
      R.nxz[j] = 0.0;
      R.nzz[j] = 0.5;
    }
  }
  return R;
}

NormalFieldCoefficients normalField(const Polygon2DLayer& l, double dx, double dz, const IndexMap& map) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  const int Gx = std::max(l.normalField.grid, 8 * (2 * map.Nx() + 1));
  const int Gz = std::max(l.normalField.grid, 8 * (2 * map.Nz() + 1));
  const NormalFieldRaster R = normalFieldRaster(l, dx, dz, Gx, Gz);
  NormalFieldCoefficients nf;
  for (auto& row : nf.NN)
    for (auto& g : row) g = fourier::CoefficientGrid(Px, Pz);
  for (auto& g : nf.N) g = fourier::CoefficientGrid(Px, Pz);
  // step 5: DFT of the smooth raster; N_y = 0 (y-invariant layer)
  nf.NN[0][0] = fourier::gridCoefficients2D(R.nxx, Gx, Gz, Px, Pz);
  nf.NN[0][2] = fourier::gridCoefficients2D(R.nxz, Gx, Gz, Px, Pz);
  nf.NN[2][0] = nf.NN[0][2];
  nf.NN[2][2] = fourier::gridCoefficients2D(R.nzz, Gx, Gz, Px, Pz);
  nf.hasSigned = false;  // no continuous signed field exists for general cross-sections
  return nf;
}

NormalFieldCoefficients normalField(const Profile2DLayer& l, double dx, double dz, const IndexMap& map) {
  if (l.shape == Profile2DLayer::Shape::SinusX) {
    // z-invariant: identical to the 1D sinus (V9)
    Profile1DLayer p;
    p.shape = Profile1DLayer::Shape::Sinus;
    p.depth_um = l.depth_um;
    p.x0_um = l.x0_um;
    return normalField(p, dx, map);
  }
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  const int Mx = std::max(256, 16 * (Px + 1)), Mz = std::max(256, 16 * (Pz + 1));
  const double Kx = 2 * PI / dx, Kz = 2 * PI / dz, a = l.depth_um;
  auto grad = [&](double x, double z) -> std::array<double, 2> {
    if (l.shape == Profile2DLayer::Shape::EggCrate) {
      const double cx = 1 + std::cos(Kx * (x - l.x0_um)), cz = 1 + std::cos(Kz * (z - l.z0_um));
      return {-0.25 * a * Kx * std::sin(Kx * (x - l.x0_um)) * cz, -0.25 * a * Kz * std::sin(Kz * (z - l.z0_um)) * cx};
    }
    FG_REQUIRE(static_cast<bool>(l.g), "profile2d: function missing");
    if (l.grad) return l.grad(x, z);
    const double hx = 1e-6 * dx, hz = 1e-6 * dz;
    auto gw = [&](double u, double w) { return l.g(std::fmod(std::fmod(u, dx) + dx, dx), std::fmod(std::fmod(w, dz) + dz, dz)); };
    return {(gw(x + hx, z) - gw(x - hx, z)) / (2 * hx), (gw(x, z + hz) - gw(x, z - hz)) / (2 * hz)};
  };
  const std::size_t n = static_cast<std::size_t>(Mx) * static_cast<std::size_t>(Mz);
  std::array<std::vector<double>, 3> N;
  for (auto& v : N) v.resize(n);
  for (int i = 0; i < Mx; ++i)
    for (int k = 0; k < Mz; ++k) {
      const auto gr = grad(dx * i / Mx, dz * k / Mz);
      const double w = 1.0 / std::sqrt(1.0 + gr[0] * gr[0] + gr[1] * gr[1]);
      const std::size_t j = static_cast<std::size_t>(i) * static_cast<std::size_t>(Mz) + static_cast<std::size_t>(k);
      N[0][j] = -gr[0] * w;
      N[1][j] = w;
      N[2][j] = -gr[1] * w;
    }
  NormalFieldCoefficients nf;
  nf.hasSigned = true;
  std::vector<double> prod(n);
  for (std::size_t i = 0; i < 3; ++i) {
    nf.N[i] = fourier::gridCoefficients2D(N[i], Mx, Mz, Px, Pz);
    for (std::size_t j = i; j < 3; ++j) {
      for (std::size_t q = 0; q < n; ++q) prod[q] = N[i][q] * N[j][q];
      nf.NN[i][j] = fourier::gridCoefficients2D(prod, Mx, Mz, Px, Pz);
      nf.NN[j][i] = nf.NN[i][j];
    }
  }
  return nf;
}

}  // namespace fg
