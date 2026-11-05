#include "fg/geometry/normal_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "fg/core/error.hpp"
#include "fg/core/parallel.hpp"
#include "fg/core/prof.hpp"
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
  // periodic neighbour indices (same summation order as the direct modulo form)
  auto wrap = [](int n, int R) {
    std::vector<int> w(static_cast<std::size_t>(n) * static_cast<std::size_t>(2 * R + 1));
    for (int k = 0; k < n; ++k)
      for (int r = -R; r <= R; ++r) w[static_cast<std::size_t>(k) * static_cast<std::size_t>(2 * R + 1) + static_cast<std::size_t>(r + R)] = ((k + r) % n + n) % n;
    return w;
  };
  const std::vector<int> nz = wrap(Gz, Rz), nx = wrap(Gx, Rx);
  // rows are independent (Plan section "M15", step 4); every value keeps its summation order
  const int threads = innerThreads();
  parallelFor(Gx, threads, [&](int i) {
    for (int k = 0; k < Gz; ++k) {
      double s = 0;
      const int* nk = &nz[static_cast<std::size_t>(k) * static_cast<std::size_t>(2 * Rz + 1)];
      for (int r = -Rz; r <= Rz; ++r) s += wz[static_cast<std::size_t>(r + Rz)] * f[at(i, nk[r + Rz])];
      tmp[at(i, k)] = s;
    }
  });
  parallelFor(Gx, threads, [&](int i) {
    const int* ni = &nx[static_cast<std::size_t>(i) * static_cast<std::size_t>(2 * Rx + 1)];
    for (int k = 0; k < Gz; ++k) {
      double s = 0;
      for (int r = -Rx; r <= Rx; ++r) s += wx[static_cast<std::size_t>(r + Rx)] * tmp[at(ni[r + Rx], k)];
      f[at(i, k)] = s;
    }
  });
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
  FG_PROF_SCOPE("nf.raster");
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
  parallelFor(Gx, innerThreads(), [&](int i) {  // rows are independent (M15)
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
  });
  R.nx.assign(n, 1.0);
  R.nz.assign(n, 0.0);
  if (segs.empty()) return R;
  // step 3: periodic Gaussian with the same physical width along x and z
  const double hphys = l.normalField.sigma * std::min(dx / Gx, dz / Gz);
  const double sx = hphys * Gx / dx, sz = hphys * Gz / dz;
  for (auto* f : {&R.nxx, &R.nxz, &R.nzz}) gaussSmooth(*f, Gx, Gz, sx, sz);
  gaussSmooth(vx, Gx, Gz, sx, sz);
  gaussSmooth(vz, Gx, Gz, sx, sz);
  double gmax = 0;
  for (std::size_t j = 0; j < n; ++j) gmax = std::max(gmax, std::hypot(vx[j], vz[j]));
  // rank-1 field for the projector variant B (M8): normalized smoothed signed gradient; where it
  // vanishes, the principal eigenvector of the smoothed tensor
  for (std::size_t j = 0; j < n; ++j) {
    const double gn = std::hypot(vx[j], vz[j]);
    if (gn > l.normalField.tau * gmax) {
      R.nx[j] = vx[j] / gn;
      R.nz[j] = vz[j] / gn;
    } else {
      const double a = R.nxx[j], b = R.nxz[j], c = R.nzz[j];
      const double th = 0.5 * std::atan2(2 * b, a - c);  // principal axis of [[a, b], [b, c]]
      R.nx[j] = std::cos(th);
      R.nz[j] = std::sin(th);
    }
  }
  if (!vec) return R;
  for (std::size_t j = 0; j < n; ++j) {
    const double gn = std::hypot(vx[j], vz[j]);
    if (gn > l.normalField.tau * gmax) {  // step 4, vector method
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

int normalFieldGrid(const NormalFieldOptions& o, int N) {
  int G = std::max(o.grid, 8 * (2 * N + 1));
  int p = 1;
  while (p < G) p <<= 1;  // power of two (FFT)
  return p;
}

NormalFieldCoefficients normalField(const Polygon2DLayer& l, double dx, double dz, const IndexMap& map) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  const int Gx = normalFieldGrid(l.normalField, map.Nx());
  const int Gz = normalFieldGrid(l.normalField, map.Nz());
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

// ------------------------------------------------------------------ Algorithm NF from boundary samples (M13)
namespace {

/// Chord of the boundary (translated so that its midpoint lies in the cell) with the unit normals
/// at both ends.
struct Chord {
  double ax, az, bx, bz;
  Vec3D na, nb;
};

/// Unit normal at parameter t of a chord: the end normal if both are equal (polygon edges: exact),
/// otherwise linear interpolation, normalized.
Vec3D chordNormal(const Chord& c, double t) {
  if (c.na == c.nb) return c.na;
  const Vec3D n = (1 - t) * c.na + t * c.nb;
  return n / n.norm();
}

/// Material boundaries among the chords of the chains (Plan alg:nfield step 1). Exact chords: as
/// Algorithm NF for polygons (64 sub-pieces, probes 1e-9 periods off the chord); approximating
/// chords: one probe pair at the chord midpoint, 1e-2 raster spacings off (sagitta <= 1e-3).
std::vector<Chord> materialChords(const std::vector<BoundaryChain>& chains,
                                  const std::function<int(double, double)>& regionAt,
                                  const std::function<std::string(int)>& materialOf, double dx, double dz,
                                  double h) {
  const double offExact = 1e-9 * std::min(dx, dz), offCurve = 1e-2 * h;
  constexpr int S = 64;
  std::vector<Chord> out;
  auto differ = [&](double mx, double mz, double nx, double nz, double off) {
    return materialOf(regionAt(mx - off * nx, mz - off * nz)) != materialOf(regionAt(mx + off * nx, mz + off * nz));
  };
  auto push = [&](Chord c) {
    const double sx = std::floor(0.5 * (c.ax + c.bx) / dx) * dx, sz = std::floor(0.5 * (c.az + c.bz) / dz) * dz;
    c.ax -= sx; c.bx -= sx; c.az -= sz; c.bz -= sz;
    out.push_back(c);
  };
  for (const BoundaryChain& ch : chains) {
    const std::size_t n = ch.samples.size();
    for (std::size_t j = 0; j + (ch.closed ? 0 : 1) < n; ++j) {
      const BoundarySample& a = ch.samples[j];
      const BoundarySample& b = ch.samples[(j + 1) % n];
      const double tx = b.x - a.x, tz = b.z - a.z, len = std::hypot(tx, tz);
      if (len == 0.0) continue;
      if (a.exactChord) {
        const double nx = tz / len, nz = -tx / len;  // as boundarySegments (polygon edges)
        std::array<bool, S> isB{};
        for (int s = 0; s < S; ++s) {
          const double t = (s + 0.5) / S;
          isB[static_cast<std::size_t>(s)] = differ(a.x + t * tx, a.z + t * tz, nx, nz, offExact);
        }
        for (int s = 0; s < S;) {
          if (!isB[static_cast<std::size_t>(s)]) { ++s; continue; }
          int e = s;
          while (e < S && isB[static_cast<std::size_t>(e)]) ++e;
          // exact chord = flat piece of the boundary: constant normal (that of its start sample;
          // the next sample may belong to the following edge)
          push(Chord{a.x + tx * s / S, a.z + tz * s / S, a.x + tx * e / S, a.z + tz * e / S, a.n, a.n});
          s = e;
        }
      } else {
        const Vec3D nm = chordNormal(Chord{0, 0, 0, 0, a.n, b.n}, 0.5);
        const double nl = std::hypot(nm(0), nm(2));
        if (nl > 0 && differ(a.x + 0.5 * tx, a.z + 0.5 * tz, nm(0) / nl, nm(2) / nl, offCurve))
          push(Chord{a.x, a.z, b.x, b.z, a.n, b.n});
      }
    }
  }
  return out;
}

/// Bounding volume hierarchy over the periodic images (di, dj in -1..1) of all chords (Plan M13
/// step 2; replaces the bucket grid: the query cost does not grow with the distance to the boundary).
class ChordTree {
 public:
  struct Leaf {
    int s, di, dj;
    double x0, z0, x1, z1;  ///< bounding box of the image
  };
  ChordTree(const std::vector<Chord>& c, double dx, double dz) {
    for (std::size_t k = 0; k < c.size(); ++k)
      for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
          leaves_.push_back({static_cast<int>(k), di, dj, std::min(c[k].ax, c[k].bx) + di * dx,
                             std::min(c[k].az, c[k].bz) + dj * dz, std::max(c[k].ax, c[k].bx) + di * dx,
                             std::max(c[k].az, c[k].bz) + dj * dz});
    if (!leaves_.empty()) build(0, leaves_.size());
  }
  /// Calls f(leaf) for every leaf whose box is closer than bound(); bound() may shrink meanwhile.
  template <typename F, typename B>
  void visit(double x, double z, F&& f, B&& bound) const {
    if (nodes_.empty()) return;
    std::size_t stack[128];  // depth <= log2(leaves) + 1
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
      const Node& n = nodes_[stack[--top]];
      const double b2 = bound() * bound();
      if (boxDist2(n.x0, n.z0, n.x1, n.z1, x, z) > b2) continue;
      if (n.left == 0) {
        for (std::size_t k = n.begin; k < n.end; ++k) {
          const Leaf& l = leaves_[k];
          if (boxDist2(l.x0, l.z0, l.x1, l.z1, x, z) <= bound() * bound()) f(l);
        }
        continue;
      }
      // nearer child last (visited first)
      const Node& a = nodes_[n.left];
      const Node& b = nodes_[n.left + 1];
      const bool aFirst = boxDist2(a.x0, a.z0, a.x1, a.z1, x, z) <= boxDist2(b.x0, b.z0, b.x1, b.z1, x, z);
      stack[top++] = aFirst ? n.left + 1 : n.left;
      stack[top++] = aFirst ? n.left : n.left + 1;
    }
  }

 private:
  struct Node {
    double x0, z0, x1, z1;
    std::size_t begin, end, left;  ///< leaves [begin, end); children left, left + 1 (0: leaf node)
  };
  /// squared distance to a box (pruning only; with a relative margin of 1e-12 in bound())
  static double boxDist2(double x0, double z0, double x1, double z1, double x, double z) {
    const double ex = std::max(std::max(x0 - x, x - x1), 0.0), ez = std::max(std::max(z0 - z, z - z1), 0.0);
    return ex * ex + ez * ez;
  }
  std::size_t build(std::size_t b, std::size_t e) {
    Node n{1e300, 1e300, -1e300, -1e300, b, e, 0};
    for (std::size_t k = b; k < e; ++k) {
      n.x0 = std::min(n.x0, leaves_[k].x0);
      n.z0 = std::min(n.z0, leaves_[k].z0);
      n.x1 = std::max(n.x1, leaves_[k].x1);
      n.z1 = std::max(n.z1, leaves_[k].z1);
    }
    const std::size_t id = nodes_.size();
    nodes_.push_back(n);
    if (e - b <= 4) return id;
    const bool alongX = n.x1 - n.x0 >= n.z1 - n.z0;
    const std::size_t mid = b + (e - b) / 2;
    std::nth_element(leaves_.begin() + static_cast<std::ptrdiff_t>(b), leaves_.begin() + static_cast<std::ptrdiff_t>(mid),
                     leaves_.begin() + static_cast<std::ptrdiff_t>(e), [alongX](const Leaf& u, const Leaf& v) {
                       return alongX ? u.x0 + u.x1 < v.x0 + v.x1 : u.z0 + u.z1 < v.z0 + v.z1;
                     });
    const std::size_t l = nodes_.size();
    nodes_.push_back({});
    nodes_.push_back({});
    const std::size_t lc = build(b, mid), rc = build(mid, e);
    nodes_[l] = nodes_[lc];
    nodes_[l + 1] = nodes_[rc];
    nodes_[id].left = l;
    return id;
  }
  std::vector<Leaf> leaves_;
  std::vector<Node> nodes_;
};

}  // namespace

std::vector<BoundaryChain> polygonChains(const std::vector<std::array<double, 2>>& vertices) {
  const Poly v = counterClockwise(vertices);
  std::vector<BoundaryChain> out;
  for (std::size_t j = 0; j < v.size(); ++j) {
    const auto& p0 = v[j];
    const auto& p1 = v[(j + 1) % v.size()];
    const double tx = p1[0] - p0[0], tz = p1[1] - p0[1], len = std::hypot(tx, tz);
    if (len == 0.0) continue;
    const Vec3D n(tz / len, 0.0, -tx / len);
    BoundaryChain c;
    c.closed = false;
    c.samples = {{p0[0], p0[1], n, true}, {p1[0], p1[1], n, false}};
    out.push_back(std::move(c));
  }
  return out;
}

NormalFieldRaster normalFieldRaster(const std::vector<BoundaryChain>& chains,
                                    const std::function<int(double, double)>& regionAt,
                                    const std::function<std::string(int)>& materialOf, const NormalFieldOptions& opt,
                                    double dx, double dz, int Gx, int Gz) {
  FG_PROF_SCOPE("nf.raster");
  FG_REQUIRE(Gx > 0 && Gz > 0 && dx > 0 && dz > 0, "normal field raster: invalid size");
  const bool vec = opt.method == NormalFieldOptions::Method::Vector;
  const double h = std::min(dx / Gx, dz / Gz);
  const std::vector<Chord> segs = materialChords(chains, regionAt, materialOf, dx, dz, h);
  bool threeD = false;
  for (const Chord& c : segs) threeD = threeD || c.na(1) != 0.0 || c.nb(1) != 0.0;
  const std::size_t n = static_cast<std::size_t>(Gx) * static_cast<std::size_t>(Gz);
  NormalFieldRaster R;
  R.Gx = Gx;
  R.Gz = Gz;
  R.nxx.assign(n, 0.5);
  R.nxz.assign(n, 0.0);
  R.nzz.assign(n, 0.5);
  R.region.assign(n, -1);
  if (threeD) {
    R.nxy.assign(n, 0.0);
    R.nyy.assign(n, 0.0);
    R.nyz.assign(n, 0.0);
  }
  std::vector<double> vx(n, 0.0), vy(n, 0.0), vz(n, 0.0);
  const double L = std::max(dx, dz), tie = 1e-12 * L;
  const ChordTree tree(segs, dx, dz);
  struct Cand {
    int s, di, dj;
    double d, t, rx, rz;
  };
  const double stepZ = dz / Gz * (1 + 1e-9);
  parallelFor(Gx, innerThreads(), [&](int i) {  // rows are independent (M15): prev restarts per row
    std::vector<Cand> cand;
    double prev = std::numeric_limits<double>::infinity();  // nearest distance of the previous point
    for (int k = 0; k < Gz; ++k) {
      const std::size_t idx = static_cast<std::size_t>(i) * static_cast<std::size_t>(Gz) + static_cast<std::size_t>(k);
      const double x = dx * i / Gx, z = dz * k / Gz;
      const int reg = regionAt(x, z);
      R.region[idx] = reg;
      if (segs.empty()) continue;
      const double sgn = reg < 0 ? 1.0 : -1.0;
      // candidates: all chord images within the nearest distance + margin (BVH, best first)
      cand.clear();
      double dmin = std::numeric_limits<double>::infinity();
      // the nearest distance is at most that of the previous point plus the step (pruning only)
      const double dup = prev + stepZ;
      auto bound = [&] { return (std::min(dmin, dup) + 9 * tie) * (1 + 1e-12); };
      tree.visit(
          x, z,
          [&](const ChordTree::Leaf& e) {
            const Chord& g = segs[static_cast<std::size_t>(e.s)];
            // identical arithmetic to normalFieldRaster(Polygon2DLayer)
            const double ax0 = g.ax + e.di * dx, az0 = g.az + e.dj * dz, bx = g.bx + e.di * dx, bz = g.bz + e.dj * dz;
            const double tx = bx - ax0, tz = bz - az0;
            const double t = std::clamp(((x - ax0) * tx + (z - az0) * tz) / (tx * tx + tz * tz), 0.0, 1.0);
            const double rx = x - (ax0 + t * tx), rz = z - (az0 + t * tz);
            const double lim = bound();
            if (rx * rx + rz * rz > lim * lim) return;  // cheap rejection before hypot
            const double d = std::hypot(rx, rz);
            if (d > dmin + 8 * tie) return;
            dmin = std::min(dmin, d);
            cand.push_back({e.s, e.di, e.dj, d, t, rx, rz});
          },
          bound);
      prev = dmin;
      // the running tie logic of Algorithm NF in the original order (segment, di, dj) on the near
      // candidates: identical result (farther candidates are always discarded by it)
      std::sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b) {
        return a.s != b.s ? a.s < b.s : (a.di != b.di ? a.di < b.di : a.dj < b.dj);
      });
      double best = std::numeric_limits<double>::infinity();
      double axx = 0, axz = 0, azz = 0, axy = 0, ayy = 0, ayz = 0, ax = 0, ay = 0, az = 0;
      int cnt = 0;
      for (const Cand& c : cand) {
        if (c.d > dmin + 8 * tie) continue;
        if (c.d > best + tie) continue;
        const Chord& g = segs[static_cast<std::size_t>(c.s)];
        double ux, uy, uz;
        if (c.d > tie) {
          ux = sgn * c.rx / c.d;
          uz = sgn * c.rz / c.d;
          uy = 0.0;
          if (threeD) {
            uy = g.na == g.nb ? g.na(1) : (1 - c.t) * g.na(1) + c.t * g.nb(1);
            const double s = std::sqrt(std::max(0.0, 1.0 - uy * uy));
            ux *= s;
            uz *= s;
          }
        } else {
          const Vec3D nn = chordNormal(g, c.t);
          ux = nn(0);
          uy = nn(1);
          uz = nn(2);
        }
        if (c.d < best - tie) {
          best = c.d;
          axx = axz = azz = axy = ayy = ayz = ax = ay = az = 0;
          cnt = 0;
        }
        axx += ux * ux;
        axz += ux * uz;
        azz += uz * uz;
        axy += ux * uy;
        ayy += uy * uy;
        ayz += uy * uz;
        ax += ux;
        ay += uy;
        az += uz;
        ++cnt;
      }
      R.nxx[idx] = axx / cnt;
      R.nxz[idx] = axz / cnt;
      R.nzz[idx] = azz / cnt;
      if (threeD) {
        R.nxy[idx] = axy / cnt;
        R.nyy[idx] = ayy / cnt;
        R.nyz[idx] = ayz / cnt;
      }
      vx[idx] = ax / cnt;
      vy[idx] = ay / cnt;
      vz[idx] = az / cnt;
    }
  });
  R.nx.assign(n, 1.0);
  R.nz.assign(n, 0.0);
  if (threeD) R.ny.assign(n, 0.0);
  if (segs.empty()) return R;
  // step 3: periodic Gaussian with the same physical width along x and z
  const double hphys = opt.sigma * std::min(dx / Gx, dz / Gz);
  const double sx = hphys * Gx / dx, sz = hphys * Gz / dz;
  for (auto* f : {&R.nxx, &R.nxz, &R.nzz}) gaussSmooth(*f, Gx, Gz, sx, sz);
  gaussSmooth(vx, Gx, Gz, sx, sz);
  gaussSmooth(vz, Gx, Gz, sx, sz);
  if (threeD) {
    for (auto* f : {&R.nxy, &R.nyy, &R.nyz}) gaussSmooth(*f, Gx, Gz, sx, sz);
    gaussSmooth(vy, Gx, Gz, sx, sz);
  }
  if (!threeD) {
    // identical to normalFieldRaster(Polygon2DLayer)
    double gmax = 0;
    for (std::size_t j = 0; j < n; ++j) gmax = std::max(gmax, std::hypot(vx[j], vz[j]));
    for (std::size_t j = 0; j < n; ++j) {
      const double gn = std::hypot(vx[j], vz[j]);
      if (gn > opt.tau * gmax) {
        R.nx[j] = vx[j] / gn;
        R.nz[j] = vz[j] / gn;
      } else {
        const double a = R.nxx[j], b = R.nxz[j], c = R.nzz[j];
        const double th = 0.5 * std::atan2(2 * b, a - c);
        R.nx[j] = std::cos(th);
        R.nz[j] = std::sin(th);
      }
    }
    if (!vec) return R;
    for (std::size_t j = 0; j < n; ++j) {
      const double gn = std::hypot(vx[j], vz[j]);
      if (gn > opt.tau * gmax) {
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
  // 3D: rank-1 field from the smoothed vector, else the principal axis of the smoothed tensor
  double gmax = 0;
  for (std::size_t j = 0; j < n; ++j) gmax = std::max(gmax, Vec3D(vx[j], vy[j], vz[j]).norm());
  for (std::size_t j = 0; j < n; ++j) {
    const Vec3D v(vx[j], vy[j], vz[j]);
    const double gn = v.norm();
    Vec3D b;
    if (gn > opt.tau * gmax) {
      b = v / gn;
    } else {
      Eigen::Matrix3d P;
      P << R.nxx[j], R.nxy[j], R.nxz[j], R.nxy[j], R.nyy[j], R.nyz[j], R.nxz[j], R.nyz[j], R.nzz[j];
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(P);
      b = es.eigenvectors().col(2);
    }
    R.nx[j] = b(0);
    R.ny[j] = b(1);
    R.nz[j] = b(2);
    if (vec) {
      const bool ok = gn > opt.tau * gmax;
      R.nxx[j] = ok ? b(0) * b(0) : 0.5;
      R.nxz[j] = ok ? b(0) * b(2) : 0.0;
      R.nzz[j] = ok ? b(2) * b(2) : 0.5;
      R.nxy[j] = ok ? b(0) * b(1) : 0.0;
      R.nyy[j] = ok ? b(1) * b(1) : 0.0;
      R.nyz[j] = ok ? b(1) * b(2) : 0.0;
    }
  }
  return R;
}

NormalFieldCoefficients rasterCoefficients(const NormalFieldRaster& R, const IndexMap& map) {
  FG_PROF_SCOPE("nf.coefficients");
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  NormalFieldCoefficients nf;
  for (auto& row : nf.NN)
    for (auto& g : row) g = fourier::CoefficientGrid(Px, Pz);
  for (auto& g : nf.N) g = fourier::CoefficientGrid(Px, Pz);
  nf.NN[0][0] = fourier::gridCoefficients2D(R.nxx, R.Gx, R.Gz, Px, Pz);
  nf.NN[0][2] = fourier::gridCoefficients2D(R.nxz, R.Gx, R.Gz, Px, Pz);
  nf.NN[2][0] = nf.NN[0][2];
  nf.NN[2][2] = fourier::gridCoefficients2D(R.nzz, R.Gx, R.Gz, Px, Pz);
  if (R.has3D()) {
    nf.NN[0][1] = nf.NN[1][0] = fourier::gridCoefficients2D(R.nxy, R.Gx, R.Gz, Px, Pz);
    nf.NN[1][1] = fourier::gridCoefficients2D(R.nyy, R.Gx, R.Gz, Px, Pz);
    nf.NN[1][2] = nf.NN[2][1] = fourier::gridCoefficients2D(R.nyz, R.Gx, R.Gz, Px, Pz);
  }
  nf.hasSigned = false;
  return nf;
}

std::array<double, 2> profile2DGradient(const Profile2DLayer& l, double dx, double dz, double x, double z) {
  const double Kx = 2 * PI / dx, Kz = 2 * PI / dz, a = l.depth_um;
  switch (l.shape) {
    case Profile2DLayer::Shape::SinusX:
      return {-0.5 * a * Kx * std::sin(Kx * (x - l.x0_um)), 0.0};
    case Profile2DLayer::Shape::EggCrate: {
      const double cx = 1 + std::cos(Kx * (x - l.x0_um)), cz = 1 + std::cos(Kz * (z - l.z0_um));
      return {-0.25 * a * Kx * std::sin(Kx * (x - l.x0_um)) * cz, -0.25 * a * Kz * std::sin(Kz * (z - l.z0_um)) * cx};
    }
    case Profile2DLayer::Shape::Function:
      break;
  }
  FG_REQUIRE(static_cast<bool>(l.g), "profile2d: function missing");
  if (l.grad) return l.grad(x, z);
  const double hx = 1e-6 * dx, hz = 1e-6 * dz;
  auto gw = [&](double u, double w) { return l.g(std::fmod(std::fmod(u, dx) + dx, dx), std::fmod(std::fmod(w, dz) + dz, dz)); };
  return {(gw(x + hx, z) - gw(x - hx, z)) / (2 * hx), (gw(x, z + hz) - gw(x, z - hz)) / (2 * hz)};
}

Vec3D profile2DNormal(const Profile2DLayer& l, double dx, double dz, double x, double z) {
  const auto g = profile2DGradient(l, dx, dz, x, z);
  return Vec3D(-g[0], 1.0, -g[1]) / std::sqrt(1.0 + g[0] * g[0] + g[1] * g[1]);
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
  auto grad = [&](double x, double z) { return profile2DGradient(l, dx, dz, x, z); };
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
