#include "layer_model.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/fourier/contour.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/geometry/profile.hpp"
#include "fg/physics/anisotropic.hpp"

namespace fg::detail {

const Material& LayerModel::material() const { throw Error("layer is not homogeneous"); }

namespace {

using Grid = fourier::CoefficientGrid;
using GridMat = std::array<std::array<Grid, 3>, 3>;

/// N N^T coefficients of a constant normal.
NormalFieldCoefficients constantNN(const Vec3D& N, const IndexMap& map) {
  NormalFieldCoefficients nf;
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  for (int i = 0; i < 3; ++i) {
    nf.N[static_cast<std::size_t>(i)] = Grid(Px, Pz);
    nf.N[static_cast<std::size_t>(i)](0, 0) = N(i);
    for (int j = 0; j < 3; ++j) {
      auto& g = nf.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      g = Grid(Px, Pz);
      g(0, 0) = N(i) * N(j);
    }
  }
  nf.hasSigned = true;
  return nf;
}

/// 1D piecewise-constant isotropic coefficients eps, 1/eps on intervals.
void isoGrids(const std::vector<double>& br, const std::vector<cplx>& eps, double dx, const IndexMap& map,
              SliceCoefficients& sc) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  std::vector<cplx> inv;
  for (cplx e : eps) inv.push_back(1.0 / e);
  const auto ce = fourier::piecewiseConstant1D(br, eps, dx, Px);
  const auto ci = fourier::piecewiseConstant1D(br, inv, dx, Px);
  sc.eps = Grid(Px, Pz);
  sc.invEps = Grid(Px, Pz);
  for (int p = -Px; p <= Px; ++p) {
    sc.eps(p, 0) = ce[static_cast<std::size_t>(p + Px)];
    sc.invEps(p, 0) = ci[static_cast<std::size_t>(p + Px)];
  }
}

/// Number of coefficient functions per anisotropic factorization (Laurent: eps_ij; FFF_A: C, eps C;
/// FFF_B: G1, G2, H1, P_T, 1/xi).
int anisoCount(Factorization f) {
  switch (f) {
    case Factorization::Laurent: return 9;
    case Factorization::FFF_A: return 18;
    case Factorization::FFF_B: return 37;
    default: throw Error("factorization not applicable to an anisotropic layer (use FFF_A, FFF_B or Laurent)");
  }
}

/// Pointwise values of the anisotropic coefficient functions (Plan eq:QaniA / eq:QaniB).
void anisoEval(Factorization f, const Mat3C& e, const Vec3D& N, cplx* out) {
  if (f == Factorization::Laurent) {
    for (int q = 0; q < 9; ++q) out[q] = e(q / 3, q % 3);
  } else if (f == Factorization::FFF_A) {
    Vec3D T1, T2;
    profileTangents(N, T1, T2);
    const Mat3C C = cMatrix(e, N, T1, T2);
    const Mat3C eC = e * C;
    for (int q = 0; q < 9; ++q) {
      out[q] = C(q / 3, q % 3);
      out[9 + q] = eC(q / 3, q % 3);
    }
  } else {
    const ProjectorSet ps = projectorMatrices(e, N);
    for (int q = 0; q < 9; ++q) {
      out[q] = ps.G1(q / 3, q % 3);
      out[9 + q] = ps.G2(q / 3, q % 3);
      out[18 + q] = ps.H1(q / 3, q % 3);
      out[27 + q] = ps.PT(q / 3, q % 3);
    }
    out[36] = ps.invXi;
  }
}

/// Q blocks from the coefficient grids of the anisotropic functions (order as in anisoEval).
QBlocks anisoQFromGrids(Factorization f, std::vector<Grid>& g, const IndexMap& map) {
  SliceCoefficients sc;
  sc.isotropic = false;
  auto take = [&](GridMat& m, int off) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = std::move(g[static_cast<std::size_t>(off + 3 * i + j)]);
  };
  if (f == Factorization::Laurent) {
    take(sc.epsTensor, 0);
    return qLaurent(sc, map);
  }
  if (f == Factorization::FFF_A) {
    take(sc.C, 0);
    take(sc.epsC, 9);
    return qFFFAnisotropicA(sc, map);
  }
  take(sc.G1, 0);
  take(sc.G2, 9);
  take(sc.H1, 18);
  take(sc.PT, 27);
  sc.invXi = std::move(g[36]);
  return qFFFAnisotropicB(sc, map);
}

/// 1D coefficient vectors of the anisotropic functions on intervals with tensors eps (exact for
/// piecewise constant, Gauss quadrature for the smooth N dependence).
std::vector<std::vector<cplx>> anisoLineCoefficients(const std::vector<double>& br, const std::vector<Mat3C>& eps,
                                                     const std::function<Vec3D(double)>& Nfun, Factorization f,
                                                     double dx, int Px) {
  return fourier::piecewiseSmoothCoefficients(
      br, anisoCount(f),
      [&](double x, int k, cplx* out) {
        anisoEval(f, eps[static_cast<std::size_t>(k)], f == Factorization::Laurent ? Vec3D(1, 0, 0) : Nfun(x), out);
      },
      dx, Px);
}

/// Q blocks of an anisotropic 1D slice (Plan eq:QaniA / eq:QaniB; Laurent for comparison).
QBlocks anisotropicQ(const std::vector<double>& br, const std::vector<Mat3C>& eps,
                     const std::function<Vec3D(double)>& Nfun, Factorization fac, double dx,
                     const IndexMap& map) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  const auto c = anisoLineCoefficients(br, eps, Nfun, fac, dx, Px);
  std::vector<Grid> g(c.size(), Grid(Px, Pz));
  for (std::size_t m = 0; m < c.size(); ++m)
    for (int p = -Px; p <= Px; ++p) g[m](p, 0) = c[m][static_cast<std::size_t>(p + Px)];
  return anisoQFromGrids(fac, g, map);
}

// ------------------------------------------------------------------ homogeneous
class HomogeneousModel final : public LayerModel {
 public:
  HomogeneousModel(const HomogeneousLayer& l, const Material& m, const IndexMap& map)
      : h_(l.thickness_um), mat_(m), map_(map), nn_(constantNN(Vec3D(0, 1, 0), map)) {
    epsList = {m.isotropic ? m.eps : m.epsTensor.diagonal().mean()};
  }
  double thickness_um() const override { return h_; }
  bool homogeneous() const override { return true; }
  const Material& material() const override { return mat_; }
  QBlocks Q(double) const override {
    return mat_.isotropic ? homogeneousIsotropicQ(mat_.eps, map_.size()) : homogeneousTensorQ(mat_.epsTensor, map_.size());
  }
  Mat3C epsAt(double, double, double) const override { return mat_.tensor(); }
  Vec3D normalAt(double, double, double) const override { return Vec3D(0, 1, 0); }
  const NormalFieldCoefficients& NN(double) const override { return nn_; }

 private:
  double h_;
  Material mat_;
  IndexMap map_;
  NormalFieldCoefficients nn_;
};

// ------------------------------------------------------------------ lamellar (1D)
class LamellarModel final : public LayerModel {
 public:
  LamellarModel(const Stack& st, const Lamellar1DLayer& l, Factorization fac, const IndexMap& map)
      : l_(l), dx_(st.dx_um), nn_(constantNN(Vec3D(1, 0, 0), map)) {
    FG_REQUIRE(l.breakpoints_um.size() == l.materials.size() + 1, "lamellar: breakpoints/materials");
    bool iso = true;
    for (const auto& name : l.materials) {
      const Material& m = st.material(name);
      iso = iso && m.isotropic;
      tensors_.push_back(m.tensor());
      epsList.push_back(m.isotropic ? m.eps : m.epsTensor.diagonal().mean());
    }
    if (!iso) {
      const Factorization f = fac == Factorization::Auto ? Factorization::FFF_A : fac;
      Q_ = anisotropicQ(l.breakpoints_um, tensors_, [](double) { return Vec3D(1, 0, 0); }, f, dx_, map);
      return;
    }
    SliceCoefficients sc;
    isoGrids(l.breakpoints_um, epsList, dx_, map, sc);
    switch (fac) {
      case Factorization::Laurent: Q_ = qLaurent(sc, map); break;
      case Factorization::Auto:
      case Factorization::LiLamellar:
      case Factorization::FFF_iso:
      case Factorization::FFF_sym:
      case Factorization::FFF_rig: Q_ = qLiLamellar(sc, map); break;  // FFF with N = x (exact)
      default: throw Error("factorization not applicable to an isotropic lamellar layer");
    }
  }
  double thickness_um() const override { return l_.thickness_um; }
  QBlocks Q(double) const override { return Q_; }
  Mat3C epsAt(double x, double, double) const override {
    const double xr = std::fmod(std::fmod(x, dx_) + dx_, dx_);
    std::size_t k = 0;
    while (k + 1 < tensors_.size() && xr >= l_.breakpoints_um[k + 1]) ++k;
    return tensors_[k];
  }
  Vec3D normalAt(double, double, double) const override { return Vec3D(1, 0, 0); }
  const NormalFieldCoefficients& NN(double) const override { return nn_; }

 private:
  Lamellar1DLayer l_;
  double dx_;
  std::vector<Mat3C> tensors_;
  QBlocks Q_;
  NormalFieldCoefficients nn_;
};

// ------------------------------------------------------------------ profile y = g(x) (1D)
class ProfileModel final : public LayerModel {
 public:
  ProfileModel(const Stack& st, const Profile1DLayer& p, Factorization fac, const IndexMap& map, int defaultSlices)
      : p_(p), dx_(st.dx_um), map_(map), nf_(normalField(p, st.dx_um, map)) {
    const Material& mb = st.material(p.below);
    const Material& ma = st.material(p.above);
    aniso_ = !mb.isotropic || !ma.isotropic;
    // z-invariant profile in a 2D stack: exact (q = 0 coefficients only, M7/M8)
    epsB_ = mb.eps;
    epsA_ = ma.eps;
    tb_ = mb.tensor();
    ta_ = ma.tensor();
    epsList = {aniso_ ? tb_.diagonal().mean() : epsB_, aniso_ ? ta_.diagonal().mean() : epsA_};
    rcwaSlices = p.slices > 0 ? p.slices : defaultSlices;
    // default factorization (Plan section "Varianten der isotropen FFF")
    fac_ = fac;
    if (fac_ == Factorization::Auto)
      fac_ = aniso_ ? Factorization::FFF_A
                    : (p.shape == Profile1DLayer::Shape::Sinus ? Factorization::FFF_rig : Factorization::FFF_sym);
    if (fac_ == Factorization::FFF_iso)
      for (const Material* m : {&mb, &ma})
        if (m->eps.real() < 0 && std::abs(m->eps.imag()) < 0.5 * std::abs(m->eps.real()))
          warnings.push_back("FFF_iso (book IV.32) with a weakly absorbing metal ('" + m->name +
                             "'): known to produce unphysical results, use FFF_rig or FFF_sym");
    FG_REQUIRE(aniso_ ? (fac_ == Factorization::FFF_A || fac_ == Factorization::FFF_B || fac_ == Factorization::Laurent)
                      : (fac_ == Factorization::FFF_iso || fac_ == Factorization::FFF_sym ||
                         fac_ == Factorization::FFF_rig || fac_ == Factorization::Laurent),
               "factorization not applicable to this profile layer");
  }
  double thickness_um() const override { return p_.depth_um; }
  bool yInvariant() const override { return false; }
  QBlocks Q(double y) const override {
    const double a = p_.depth_um;
    std::vector<double> br;
    std::vector<bool> below;
    profileIntervals(p_, std::clamp(y, 1e-14 * a, a * (1 - 1e-14)), dx_, br, below);
    if (aniso_) {
      std::vector<Mat3C> tens;
      for (bool b : below) tens.push_back(b ? tb_ : ta_);
      return anisotropicQ(br, tens, [this](double x) { return normalAt(x, 0.0, 0.0); }, fac_, dx_, map_);
    }
    std::vector<cplx> ve;
    for (bool b : below) ve.push_back(b ? epsB_ : epsA_);
    SliceCoefficients sc;
    isoGrids(br, ve, dx_, map_, sc);
    sc.normal = &nf_;
    return fac_ == Factorization::Laurent   ? qLaurent(sc, map_)
           : fac_ == Factorization::FFF_rig ? qFFFIsotropicRigorous(sc, map_)
                                            : qFFFIsotropic(sc, map_, fac_ == Factorization::FFF_sym);
  }
  Mat3C epsAt(double x, double y, double) const override { return profileValue(p_, x, dx_) > y ? tb_ : ta_; }
  Vec3D normalAt(double x, double, double) const override {
    const double sl = profileSlope(p_, x, dx_);
    return Vec3D(-sl, 1.0, 0.0) / std::sqrt(1.0 + sl * sl);
  }
  const NormalFieldCoefficients& NN(double) const override { return nf_; }

 private:
  Profile1DLayer p_;
  double dx_;
  IndexMap map_;
  NormalFieldCoefficients nf_;
  bool aniso_ = false;
  cplx epsB_, epsA_;
  Mat3C tb_, ta_;
  Factorization fac_ = Factorization::Auto;
};

// ------------------------------------------------------------------ helpers (2D)
/// Weak-metal warning of the book variant (Plan section "Varianten der isotropen FFF").
void warnFFFiso(Factorization f, const std::vector<const Material*>& mats, std::vector<std::string>& warnings) {
  if (f != Factorization::FFF_iso) return;
  for (const Material* m : mats)
    if (m->isotropic && m->eps.real() < 0 && std::abs(m->eps.imag()) < 0.5 * std::abs(m->eps.real()))
      warnings.push_back("FFF_iso (book IV.32) with a weakly absorbing metal ('" + m->name +
                         "'): known to produce unphysical results, use FFF_rig or FFF_sym");
}

QBlocks isotropicQ(Factorization f, const SliceCoefficients& sc, const IndexMap& map) {
  switch (f) {
    case Factorization::Laurent: return qLaurent(sc, map);
    case Factorization::FFF_rig: return qFFFIsotropicRigorous(sc, map);
    case Factorization::FFF_iso: return qFFFIsotropic(sc, map, false);
    case Factorization::FFF_sym: return qFFFIsotropic(sc, map, true);
    default: throw Error("factorization not applicable to this isotropic layer");
  }
}

// ------------------------------------------------------------------ cross-section geometries (2D)
/// Geometry of a y-invariant cross-section for CrossSectionModel (Polygon2D, Contour2D, sections of
/// SectionedLayer): background + K inclusions/regions with exact indicator coefficients.
class CrossSectionGeometry {
 public:
  virtual ~CrossSectionGeometry() = default;
  virtual double thickness() const = 0;
  virtual const std::string& background() const = 0;
  virtual std::size_t count() const = 0;
  virtual const std::string& material(std::size_t k) const = 0;
  /// Exact coefficients of the indicator of inclusion k, |p| <= Px, |q| <= Pz.
  virtual Grid indicator(std::size_t k, int Px, int Pz) const = 0;
  /// Index of the inclusion containing (x, z) (periodic), -1 for the background.
  virtual int regionAt(double x, double z) const = 0;
  /// Algorithm NF steps 1-4 on a Gx x Gz raster.
  virtual NormalFieldRaster raster(int Gx, int Gz) const = 0;
  virtual const NormalFieldOptions& nfOptions() const = 0;
};

/// Polygon2D: the original functions (eq:polygon, normalFieldRaster(Polygon2DLayer)).
class PolygonGeometry final : public CrossSectionGeometry {
 public:
  PolygonGeometry(Polygon2DLayer l, double dx, double dz) : l_(std::move(l)), dx_(dx), dz_(dz) {
    for (auto& inc : l_.inclusions) {
      FG_REQUIRE(inc.vertices.size() >= 3, "polygon2d: inclusion needs at least 3 vertices");
      // counter-clockwise orientation (Plan eq:polygon)
      double a = 0;
      for (std::size_t j = 0; j < inc.vertices.size(); ++j) {
        const auto& p0 = inc.vertices[j];
        const auto& p1 = inc.vertices[(j + 1) % inc.vertices.size()];
        a += p0[0] * p1[1] - p1[0] * p0[1];
      }
      if (a < 0) std::reverse(inc.vertices.begin(), inc.vertices.end());
    }
  }
  double thickness() const override { return l_.thickness_um; }
  const std::string& background() const override { return l_.background; }
  std::size_t count() const override { return l_.inclusions.size(); }
  const std::string& material(std::size_t k) const override { return l_.inclusions[k].material; }
  Grid indicator(std::size_t k, int Px, int Pz) const override {
    Grid chi(Px, Pz);
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q) chi(p, q) = fourier::polygonCoefficient(l_.inclusions[k].vertices, dx_, dz_, p, q);
    return chi;
  }
  int regionAt(double x, double z) const override { return polygonRegion(l_, dx_, dz_, x, z); }
  NormalFieldRaster raster(int Gx, int Gz) const override { return normalFieldRaster(l_, dx_, dz_, Gx, Gz); }
  const NormalFieldOptions& nfOptions() const override { return l_.normalField; }

 private:
  Polygon2DLayer l_;
  double dx_, dz_;
};

/// Regions with curved boundaries (Contour2D, sections of SectionedLayer): eq:contour, Algorithm NF
/// from boundary samples (given, or sampled from the contours with n_y = 0).
class ContourGeometry final : public CrossSectionGeometry {
 public:
  ContourGeometry(double thickness, std::string background, std::vector<geometry::Region> regions,
                  std::vector<BoundaryChain> chains, NormalFieldOptions nf, double dx, double dz, const char* name)
      : h_(thickness), bg_(std::move(background)), regions_(std::move(regions)), chains_(std::move(chains)),
        nf_(nf), dx_(dx), dz_(dz) {
    const double tol = 1e-6 * std::min(dx, dz);  // flattening for contains() (NF probes: 1e-2 raster spacings)
    for (auto& r : regions_) r.normalize(tol);
    checkDisjoint(name);
  }
  double thickness() const override { return h_; }
  const std::string& background() const override { return bg_; }
  std::size_t count() const override { return regions_.size(); }
  const std::string& material(std::size_t k) const override { return regions_[k].material; }
  Grid indicator(std::size_t k, int Px, int Pz) const override {
    return fourier::regionCoefficients(regions_[k], dx_, dz_, Px, Pz);  // Plan eq:contour
  }
  int regionAt(double x, double z) const override {
    for (std::size_t k = 0; k < regions_.size(); ++k) {
      const auto b = regions_[k].bounds();
      for (int i = static_cast<int>(std::ceil((x - b[2]) / dx_)); i <= static_cast<int>(std::floor((x - b[0]) / dx_)); ++i)
        for (int j = static_cast<int>(std::ceil((z - b[3]) / dz_)); j <= static_cast<int>(std::floor((z - b[1]) / dz_)); ++j)
          if (regions_[k].contains(x - i * dx_, z - j * dz_)) return static_cast<int>(k);
    }
    return -1;
  }
  NormalFieldRaster raster(int Gx, int Gz) const override {
    const double h = std::min(dx_ / Gx, dz_ / Gz);
    const std::vector<BoundaryChain> chains = chains_.empty() ? contourChains(1e-3 * h) : chains_;
    return normalFieldRaster(
        chains, [this](double x, double z) { return regionAt(x, z); },
        [this](int k) { return k < 0 ? bg_ : regions_[static_cast<std::size_t>(k)].material; }, nf_, dx_, dz_, Gx, Gz);
  }
  const NormalFieldOptions& nfOptions() const override { return nf_; }

 private:
  /// Boundary samples of the contours (Plan M13, NF step 1): pieces until the control polygon
  /// deviates <= tol from the chord; normal (z', 0, -x')/|.| (outward for the normalized orientation).
  std::vector<BoundaryChain> contourChains(double tol) const {
    std::vector<BoundaryChain> out;
    std::function<void(const geometry::BezierSegment&, BoundaryChain&, int)> sample =
        [&](const geometry::BezierSegment& s, BoundaryChain& c, int depth) {
          double dev = 0;
          const auto& a = s.P.front();
          const auto& b = s.P.back();
          const double tx = b[0] - a[0], tz = b[1] - a[1], len = std::hypot(tx, tz);
          for (std::size_t j = 1; j + 1 < s.P.size(); ++j)
            dev = std::max(dev, len > 0 ? std::abs((s.P[j][0] - a[0]) * tz - (s.P[j][1] - a[1]) * tx) / len
                                        : std::hypot(s.P[j][0] - a[0], s.P[j][1] - a[1]));
          if (s.degree() > 1 && dev > tol && depth < 40) {
            geometry::BezierSegment l, r;
            s.split(0.5, l, r);
            sample(l, c, depth + 1);
            sample(r, c, depth + 1);
            return;
          }
          geometry::Point2 r0, d0;
          s.eval(0.0, r0, d0);
          if (s.degree() == 1) d0 = {tx, tz};
          const double dl = std::hypot(d0[0], d0[1]);
          c.samples.push_back({a[0], a[1], Vec3D(d0[1] / dl, 0.0, -d0[0] / dl), s.degree() == 1});
        };
    auto addContour = [&](const geometry::Contour2D& con) {
      BoundaryChain c;
      for (const auto& s : con.segments) sample(s, c, 0);
      out.push_back(std::move(c));
    };
    for (const auto& r : regions_) {
      addContour(r.outer);
      for (const auto& h : r.holes) addContour(h);
    }
    return out;
  }

  /// Regions must be pairwise disjoint (incl. periodic images): no boundary crossings, no
  /// containment of one region's point in another (Plan M13 "Neuer Schichttyp Contour2D").
  void checkDisjoint(const char* name) const {
    const double tol = 1e-6 * std::min(dx_, dz_);
    std::vector<std::vector<geometry::Point2>> poly;
    for (const auto& r : regions_) poly.push_back(r.outer.flatten(tol));
    auto cross = [](const geometry::Point2& a, const geometry::Point2& b, const geometry::Point2& c,
                    const geometry::Point2& d) {
      auto o = [](const geometry::Point2& p, const geometry::Point2& q, const geometry::Point2& r) {
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
      };
      return o(a, b, c) * o(a, b, d) < 0 && o(c, d, a) * o(c, d, b) < 0;
    };
    for (std::size_t a = 0; a < regions_.size(); ++a)
      for (std::size_t b = a; b < regions_.size(); ++b)
        for (int di = -1; di <= 1; ++di)
          for (int dj = -1; dj <= 1; ++dj) {
            if (a == b && di == 0 && dj == 0) continue;
            const auto& P = poly[a];
            std::vector<geometry::Point2> Q = poly[b];
            for (auto& q : Q) q = {q[0] + di * dx_, q[1] + dj * dz_};
            const auto ba = regions_[a].bounds();
            double qb[4] = {1e300, 1e300, -1e300, -1e300};
            for (const auto& q : Q) {
              qb[0] = std::min(qb[0], q[0]); qb[1] = std::min(qb[1], q[1]);
              qb[2] = std::max(qb[2], q[0]); qb[3] = std::max(qb[3], q[1]);
            }
            if (qb[0] > ba[2] || qb[2] < ba[0] || qb[1] > ba[3] || qb[3] < ba[1]) continue;
            FG_REQUIRE(!regions_[a].contains(Q[0][0], Q[0][1]) &&
                           !regions_[b].contains(P[0][0] - di * dx_, P[0][1] - dj * dz_),
                       std::string(name) + ": regions overlap (or a region overlaps its periodic image)");
            for (std::size_t i = 0; i < P.size(); ++i)
              for (std::size_t j = 0; j < Q.size(); ++j)
                FG_REQUIRE(!cross(P[i], P[(i + 1) % P.size()], Q[j], Q[(j + 1) % Q.size()]),
                           std::string(name) + ": region boundaries intersect");
          }
  }

  double h_;
  std::string bg_;
  std::vector<geometry::Region> regions_;
  std::vector<BoundaryChain> chains_;
  NormalFieldOptions nf_;
  double dx_, dz_;
};

// ------------------------------------------------------------------ cross-section model (2D)
class CrossSectionModel final : public LayerModel {
 public:
  CrossSectionModel(const Stack& st, std::shared_ptr<const CrossSectionGeometry> geo, Factorization fac,
                    const IndexMap& map, const char* name)
      : geo_(std::move(geo)), dx_(st.dx_um), dz_(st.dz_um), map_(map) {
    FG_REQUIRE(st.is2D(), std::string(name) + " layers need a 2D stack (period_um.z > 0)");
    std::vector<const Material*> mats{&st.material(geo_->background())};
    for (std::size_t k = 0; k < geo_->count(); ++k) mats.push_back(&st.material(geo_->material(k)));
    bool aniso = false;
    for (const Material* m : mats) {
      aniso = aniso || !m->isotropic;
      tensors_.push_back(m->tensor());
      epsList.push_back(m->isotropic ? m->eps : m->epsTensor.diagonal().mean());
    }
    // defaults: FFF_sym (isotropic, M7), FFF_B (anisotropic, M8; FFF_A needs a continuous signed
    // normal field, which does not exist for general cross-sections)
    fac_ = fac == Factorization::Auto ? (aniso ? Factorization::FFF_B : Factorization::FFF_sym) : fac;
    FG_REQUIRE(aniso ? (fac_ == Factorization::Laurent || fac_ == Factorization::FFF_B)
                     : (fac_ == Factorization::Laurent || fac_ == Factorization::FFF_iso ||
                        fac_ == Factorization::FFF_sym || fac_ == Factorization::FFF_B),
               aniso ? std::string(name) + " (anisotropic): factorization must be laurent or FFF_B"
                     : std::string(name) + ": factorization must be laurent, FFF_iso, FFF_sym or FFF_B (FFF_rig/FFF_A "
                       "need a continuous signed normal field, which does not exist for general cross-sections)");
    warnFFFiso(fac_, mats, warnings);
    rank1_ = aniso || fac_ == Factorization::FFF_B;
    const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
    if (fac_ == Factorization::FFF_B) {
      Q_ = projectorQ(map);
    } else if (aniso) {
      // Laurent with tensors: exact polygon coefficients entry by entry
      std::vector<Grid> g(9, Grid(Px, Pz));
      for (int e = 0; e < 9; ++e) g[static_cast<std::size_t>(e)](0, 0) = tensors_[0](e / 3, e % 3);
      for (std::size_t k = 0; k < geo_->count(); ++k) {
        const Mat3C d = tensors_[k + 1] - tensors_[0];
        if (d.cwiseAbs().maxCoeff() == 0.0) continue;
        const Grid chi = geo_->indicator(k, Px, Pz);
        for (int e = 0; e < 9; ++e)
          for (int p = -Px; p <= Px; ++p)
            for (int q = -Pz; q <= Pz; ++q) g[static_cast<std::size_t>(e)](p, q) += d(e / 3, e % 3) * chi(p, q);
      }
      Q_ = anisoQFromGrids(Factorization::Laurent, g, map);
    } else {
      // exact coefficients: eps = eps_bg + sum_k (eps_k - eps_bg) chi_k (Plan eq:polygon)
      SliceCoefficients sc;
      sc.eps = Grid(Px, Pz);
      sc.invEps = Grid(Px, Pz);
      const cplx eb = mats[0]->eps;
      sc.eps(0, 0) = eb;
      sc.invEps(0, 0) = 1.0 / eb;
      for (std::size_t k = 0; k < geo_->count(); ++k) {
        const cplx ek = mats[k + 1]->eps;
        if (ek == eb) continue;
        const Grid chi = geo_->indicator(k, Px, Pz);
        for (int p = -Px; p <= Px; ++p)
          for (int q = -Pz; q <= Pz; ++q) {
            sc.eps(p, q) += (ek - eb) * chi(p, q);
            sc.invEps(p, q) += (1.0 / ek - 1.0 / eb) * chi(p, q);
          }
      }
      if (fac_ != Factorization::Laurent) {
        nf_ = rasterCoefficients(raster(), map);  // = normalField(Polygon2DLayer) for polygons
        nfReady_ = true;
      }
      sc.normal = &nf_;
      Q_ = isotropicQ(fac_, sc, map);
    }
  }
  double thickness_um() const override { return geo_->thickness(); }
  QBlocks Q(double) const override { return Q_; }
  Mat3C epsAt(double x, double, double z) const override {
    return tensors_[static_cast<std::size_t>(geo_->regionAt(x, z) + 1)];
  }
  /// Rank-1 field of Algorithm NF (variant B, anisotropic reconstruction): principal axis of the
  /// bilinearly interpolated projector b b^T (sign-free interpolation).
  Vec3D normalAt(double x, double, double z) const override {
    if (raster().has3D()) {  // tilted walls (SectionedLayer): principal axis of the 3 x 3 projector
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(interpolate3(x, z, true));
      return es.eigenvectors().col(2);
    }
    const auto P = interpolate(x, z, true);
    const double th = 0.5 * std::atan2(2 * P[1], P[0] - P[2]);
    return Vec3D(std::cos(th), 0.0, std::sin(th));
  }
  Eigen::Matrix3d projectorAt(double x, double y, double z) const override {
    if (rank1_) return LayerModel::projectorAt(x, y, z);
    if (raster().has3D()) return interpolate3(x, z, false);
    const auto P = interpolate(x, z, false);
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    M(0, 0) = P[0];
    M(0, 2) = M(2, 0) = P[1];
    M(2, 2) = P[2];
    return M;
  }
  bool rank1Projector() const override { return rank1_; }
  /// Coefficients of the reconstruction projector: tensor field (isotropic) or b b^T (variant B).
  const NormalFieldCoefficients& NN(double) const override {
    if (nfReady_) return nf_;
    if (!rank1_) {
      nf_ = rasterCoefficients(raster(), map_);
    } else {
      const NormalFieldRaster& R = raster();
      const std::size_t n = R.nx.size();
      std::vector<double> xx(n), xz(n), zz(n);
      for (std::size_t j = 0; j < n; ++j) {
        xx[j] = R.nx[j] * R.nx[j];
        xz[j] = R.nx[j] * R.nz[j];
        zz[j] = R.nz[j] * R.nz[j];
      }
      const int Px = 2 * map_.Nx(), Pz = 2 * map_.Nz();
      for (auto& row : nf_.NN)
        for (auto& g : row) g = Grid(Px, Pz);
      for (auto& g : nf_.N) g = Grid(Px, Pz);
      nf_.NN[0][0] = fourier::gridCoefficients2D(xx, R.Gx, R.Gz, Px, Pz);
      nf_.NN[0][2] = fourier::gridCoefficients2D(xz, R.Gx, R.Gz, Px, Pz);
      nf_.NN[2][0] = nf_.NN[0][2];
      nf_.NN[2][2] = fourier::gridCoefficients2D(zz, R.Gx, R.Gz, Px, Pz);
      if (R.has3D()) {
        std::vector<double> xy(n), yy(n), yz(n);
        for (std::size_t j = 0; j < n; ++j) {
          xy[j] = R.nx[j] * R.ny[j];
          yy[j] = R.ny[j] * R.ny[j];
          yz[j] = R.ny[j] * R.nz[j];
        }
        nf_.NN[0][1] = nf_.NN[1][0] = fourier::gridCoefficients2D(xy, R.Gx, R.Gz, Px, Pz);
        nf_.NN[1][1] = fourier::gridCoefficients2D(yy, R.Gx, R.Gz, Px, Pz);
        nf_.NN[1][2] = nf_.NN[2][1] = fourier::gridCoefficients2D(yz, R.Gx, R.Gz, Px, Pz);
      }
    }
    nfReady_ = true;
    return nf_;
  }

 private:

  /// Variant B (Plan eq:QaniB) for a cross-section, M8. The coefficient functions F(eps, N) jump
  /// with the material; with the rank-1 field N of Algorithm NF (smooth, sign-free use):
  ///   F = F(eps_bg, N) + sum_k chi_k [F(eps_k, N) - F(eps_bg, N)],
  /// smooth parts by DFT of the raster (|p| <= 2 Pf), the products with chi_k by convolution with
  /// the exact polygon coefficients (Laurent rule for smooth x discontinuous, Pf = 2P orders of the
  /// smooth factor). Exact when N is constant (z-invariant strips: Li formula).
  QBlocks projectorQ(const IndexMap& map) const {
    const int Px = 2 * map.Nx(), Pz = 2 * map.Nz(), Pfx = 2 * Px, Pfz = 2 * Pz;
    const NormalFieldRaster& R = raster();
    const int Gx = R.Gx, Gz = R.Gz;
    FG_REQUIRE(Gx > 2 * Pfx && Gz > 2 * Pfz, "polygon2d FFF_B: normal field raster too coarse");
    const std::size_t n = R.nx.size();
    const std::size_t K = geo_->count();
    std::vector<Grid> chi;
    std::vector<bool> differs;
    for (std::size_t k = 0; k < K; ++k) {
      differs.push_back((tensors_[k + 1] - tensors_[0]).cwiseAbs().maxCoeff() != 0.0);
      chi.push_back(differs.back() ? geo_->indicator(k, Px + Pfx, Pz + Pfz) : Grid());
    }
    std::vector<Grid> g(37);
    // blocks of functions (G1, G2, H1, P_T, 1/xi) to bound the raster memory
    const std::array<std::array<int, 2>, 5> blocks{{{0, 9}, {9, 18}, {18, 27}, {27, 36}, {36, 37}}};
    std::array<cplx, 37> vb{}, vk{};
    for (const auto& b : blocks) {
      const int m0 = b[0], m1 = b[1], nb = m1 - m0;
      std::vector<std::vector<cplx>> bg(static_cast<std::size_t>(nb), std::vector<cplx>(n));
      std::vector<std::vector<std::vector<cplx>>> diff(K);
      for (std::size_t k = 0; k < K; ++k)
        if (differs[k] && m0 != 27)  // P_T does not depend on the material
          diff[k].assign(static_cast<std::size_t>(nb), std::vector<cplx>(n));
      for (std::size_t j = 0; j < n; ++j) {
        const Vec3D N(R.nx[j], R.has3D() ? R.ny[j] : 0.0, R.nz[j]);
        anisoEval(Factorization::FFF_B, tensors_[0], N, vb.data());
        for (int m = m0; m < m1; ++m) bg[static_cast<std::size_t>(m - m0)][j] = vb[static_cast<std::size_t>(m)];
        for (std::size_t k = 0; k < K; ++k) {
          if (diff[k].empty()) continue;
          anisoEval(Factorization::FFF_B, tensors_[k + 1], N, vk.data());
          for (int m = m0; m < m1; ++m)
            diff[k][static_cast<std::size_t>(m - m0)][j] = vk[static_cast<std::size_t>(m)] - vb[static_cast<std::size_t>(m)];
        }
      }
      for (int m = m0; m < m1; ++m) {
        Grid c = fourier::gridCoefficients2D(bg[static_cast<std::size_t>(m - m0)], Gx, Gz, Px, Pz);
        for (std::size_t k = 0; k < K; ++k) {
          if (diff[k].empty()) continue;
          const Grid d = fourier::gridCoefficients2D(diff[k][static_cast<std::size_t>(m - m0)], Gx, Gz, Pfx, Pfz);
          const Grid cd = fourier::convolveCoefficients(d, chi[k], Px, Pz);
          for (int p = -Px; p <= Px; ++p)
            for (int q = -Pz; q <= Pz; ++q) c(p, q) += cd(p, q);
        }
        g[static_cast<std::size_t>(m)] = std::move(c);
      }
    }
    return anisoQFromGrids(Factorization::FFF_B, g, map);
  }

  /// Algorithm NF raster (computed once, also for the real-space reconstruction).
  const NormalFieldRaster& raster() const {
    if (R_.nx.empty())
      R_ = geo_->raster(normalFieldGrid(geo_->nfOptions(), map_.Nx()), normalFieldGrid(geo_->nfOptions(), map_.Nz()));
    return R_;
  }
  /// Bilinear periodic interpolation of (P_xx, P_xz, P_zz): rank-1 b b^T or the tensor field.
  /// 3D raster: bilinear periodic interpolation of the full 3 x 3 projector (rank-1 b b^T or tensor).
  Eigen::Matrix3d interpolate3(double x, double z, bool rank1) const {
    const NormalFieldRaster& R = raster();
    const double u = x / dx_ * R.Gx, w = z / dz_ * R.Gz;
    const double fu = std::floor(u), fw = std::floor(w);
    const double tu = u - fu, tw = w - fw;
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) {
        const long i = ((static_cast<long>(fu) + a) % R.Gx + R.Gx) % R.Gx;
        const long k = ((static_cast<long>(fw) + b) % R.Gz + R.Gz) % R.Gz;
        const std::size_t j = static_cast<std::size_t>(i * R.Gz + k);
        const double wt = (a ? tu : 1 - tu) * (b ? tw : 1 - tw);
        if (rank1) {
          const Vec3D n(R.nx[j], R.ny[j], R.nz[j]);
          M += wt * n * n.transpose();
        } else {
          Eigen::Matrix3d P;
          P << R.nxx[j], R.nxy[j], R.nxz[j], R.nxy[j], R.nyy[j], R.nyz[j], R.nxz[j], R.nyz[j], R.nzz[j];
          M += wt * P;
        }
      }
    return M;
  }
  std::array<double, 3> interpolate(double x, double z, bool rank1) const {
    const NormalFieldRaster& R = raster();
    const double u = x / dx_ * R.Gx, w = z / dz_ * R.Gz;
    const double fu = std::floor(u), fw = std::floor(w);
    const double tu = u - fu, tw = w - fw;
    std::array<double, 3> out{0, 0, 0};
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) {
        const long i = ((static_cast<long>(fu) + a) % R.Gx + R.Gx) % R.Gx;
        const long k = ((static_cast<long>(fw) + b) % R.Gz + R.Gz) % R.Gz;
        const std::size_t j = static_cast<std::size_t>(i * R.Gz + k);
        const double wt = (a ? tu : 1 - tu) * (b ? tw : 1 - tw);
        if (rank1) {
          out[0] += wt * R.nx[j] * R.nx[j];
          out[1] += wt * R.nx[j] * R.nz[j];
          out[2] += wt * R.nz[j] * R.nz[j];
        } else {
          out[0] += wt * R.nxx[j];
          out[1] += wt * R.nxz[j];
          out[2] += wt * R.nzz[j];
        }
      }
    return out;
  }

  std::shared_ptr<const CrossSectionGeometry> geo_;
  double dx_, dz_;
  IndexMap map_;
  std::vector<Mat3C> tensors_;  ///< background, inclusions
  Factorization fac_ = Factorization::Auto;
  bool rank1_ = false;
  mutable bool nfReady_ = false;
  mutable NormalFieldCoefficients nf_;
  mutable NormalFieldRaster R_;
  QBlocks Q_;
};

// ------------------------------------------------------------------ height-dependent cross-sections (M13)
/// SectionedLayer (Plan section "M13: Hoehenabhaengige Schichtmodelle"): one CrossSectionModel per
/// height, built on demand from the provider's cross-section and cached (key: y bit-exact; the RK
/// stage nodes and RCWA slice centres are reused by the near-field reconstruction).
class SectionedModel final : public LayerModel {
 public:
  SectionedModel(const Stack& st, const SectionedLayer& l, Factorization fac, const IndexMap& map, int defaultSlices)
      : l_(l), fac_(fac), map_(map) {
    FG_REQUIRE(st.is2D(), "sectioned layers need a 2D stack (period_um.z > 0)");
    FG_REQUIRE(static_cast<bool>(l.provider), "sectioned layer: no section provider");
    h_ = l.provider->thickness_um();
    FG_REQUIRE(h_ > 0, "sectioned layer: thickness must be positive");
    stack_.dx_um = st.dx_um;
    stack_.dz_um = st.dz_um;
    stack_.materials = st.materials;
    rcwaSlices = l.slices > 0 ? l.slices : defaultSlices;
    // permittivities for the spectral radius estimate: background and the regions at a few heights
    std::vector<std::string> names{l.background};
    for (double t : {0.02, 0.5, 0.98})
      for (const auto& r : l.provider->section(t * h_).regions)
        if (std::find(names.begin(), names.end(), r.material) == names.end()) names.push_back(r.material);
    for (const auto& n : names) {
      const Material& m = st.material(n);
      epsList.push_back(m.isotropic ? m.eps : m.epsTensor.diagonal().mean());
    }
    warnings = at(0.5 * h_).warnings;
  }
  double thickness_um() const override { return h_; }
  bool yInvariant() const override { return false; }
  QBlocks Q(double y) const override {
    if (l_.chebyshevNodes > 0 && interpolant()) return interpolate(*cheb_, y);
    return at(y).Q(0.0);
  }
  Mat3C epsAt(double x, double y, double z) const override { return at(y).epsAt(x, 0.0, z); }
  Vec3D normalAt(double x, double y, double z) const override { return at(y).normalAt(x, 0.0, z); }
  Eigen::Matrix3d projectorAt(double x, double y, double z) const override { return at(y).projectorAt(x, 0.0, z); }
  bool rank1Projector() const override { return at(0.5 * h_).rank1Projector(); }
  const NormalFieldCoefficients& NN(double y) const override { return at(y).NN(0.0); }

 private:
  const LayerModel& at(double y) const {
    const double yc = std::clamp(y, 1e-14 * h_, h_ * (1 - 1e-14));
    const auto it = cache_.find(yc);
    if (it != cache_.end()) return *it->second;
    if (cache_.size() >= 4096) cache_.clear();  // bounds the memory; references are used immediately
    CrossSection cs = l_.provider->section(yc);
    auto geo = std::make_shared<ContourGeometry>(h_, l_.background, std::move(cs.regions), std::move(cs.chains),
                                                 l_.normalField, stack_.dx_um, stack_.dz_um, "sectioned layer");
    auto m = std::make_unique<CrossSectionModel>(stack_, std::move(geo), fac_, map_, "sectioned layer");
    return *cache_.emplace(yc, std::move(m)).first->second;
  }

  /// Q blocks at K Chebyshev nodes (first kind) of the band and the barycentric weights.
  struct Chebyshev {
    std::vector<double> x, w;  ///< nodes in [-1, 1], weights
    std::vector<QBlocks> Q;
  };
  Chebyshev chebyshev(int K) const {
    Chebyshev c;
    for (int k = 0; k < K; ++k) {
      const double th = PI * (k + 0.5) / K;
      c.x.push_back(-std::cos(th));
      c.w.push_back((k % 2 ? -1.0 : 1.0) * std::sin(th));
      c.Q.push_back(at(0.5 * h_ * (1.0 + c.x.back())).Q(0.0));
    }
    return c;
  }
  static QBlocks interpolate(const Chebyshev& c, double y, double h) {
    const double x = std::clamp(2.0 * y / h - 1.0, -1.0, 1.0);
    std::vector<double> a(c.x.size());
    double sum = 0;
    for (std::size_t k = 0; k < c.x.size(); ++k) {
      if (x == c.x[k]) return c.Q[k];
      a[k] = c.w[k] / (x - c.x[k]);
      sum += a[k];
    }
    QBlocks Q = c.Q[0];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        MatC m = MatC::Zero(c.Q[0](i, j).rows(), c.Q[0](i, j).cols());
        for (std::size_t k = 0; k < c.x.size(); ++k) m += (a[k] / sum) * c.Q[k](i, j);
        Q(i, j) = std::move(m);
      }
    return Q;
  }
  QBlocks interpolate(const Chebyshev& c, double y) const { return interpolate(c, y, h_); }
  static double maxAbs(const QBlocks& Q) {
    double m = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) m = std::max(m, Q(i, j).size() ? Q(i, j).cwiseAbs().maxCoeff() : 0.0);
    return m;
  }
  /// Builds the interpolant on first use: K and 2K nodes compared at the K midpoints (1e-10 max|Q|).
  bool interpolant() const {
    if (chebReady_) return static_cast<bool>(cheb_);
    chebReady_ = true;
    for (int K = l_.chebyshevNodes; K <= 256; K *= 2) {
      const Chebyshev a = chebyshev(K), b = chebyshev(2 * K);
      double err = 0, ref = 0;
      for (int k = 0; k < K; ++k) {
        const double y = 0.5 * h_ * (1.0 - std::cos(PI * (k + 1.0) / (K + 1.0)));
        const QBlocks qa = interpolate(a, y), qb = interpolate(b, y);
        ref = std::max(ref, maxAbs(qb));
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            if (qa(i, j).size()) err = std::max(err, (qa(i, j) - qb(i, j)).cwiseAbs().maxCoeff());
      }
      if (err <= 1e-10 * ref) {
        cheb_ = std::make_unique<Chebyshev>(b);
        return true;
      }
    }
    warnings.push_back("sectioned layer: Chebyshev interpolation in y did not reach 1e-10 with 256 nodes; exact sections used");
    return false;
  }

  SectionedLayer l_;
  Factorization fac_;
  IndexMap map_;
  Stack stack_;  ///< periods and materials
  double h_ = 0;
  mutable std::map<double, std::unique_ptr<LayerModel>> cache_;
  mutable bool chebReady_ = false;
  mutable std::unique_ptr<Chebyshev> cheb_;
};

// ------------------------------------------------------------------ profile y = g(x, z) (2D)
class Profile2DModel final : public LayerModel {
 public:
  Profile2DModel(const Stack& st, const Profile2DLayer& p, Factorization fac, const IndexMap& map, int defaultSlices)
      : p_(p), dx_(st.dx_um), dz_(st.dz_um), map_(map) {
    FG_REQUIRE(st.is2D(), "profile2d layers need a 2D stack (period_um.z > 0)");
    FG_REQUIRE(p.depth_um > 0, "profile2d: depth must be positive");
    FG_REQUIRE(p.shape != Profile2DLayer::Shape::Function || static_cast<bool>(p.g), "profile2d: function missing");
    const Material& mb = st.material(p.below);
    const Material& ma = st.material(p.above);
    aniso_ = !mb.isotropic || !ma.isotropic;
    epsB_ = mb.eps;
    epsA_ = ma.eps;
    tb_ = mb.tensor();
    ta_ = ma.tensor();
    epsList = {aniso_ ? tb_.diagonal().mean() : epsB_, aniso_ ? ta_.diagonal().mean() : epsA_};
    rcwaSlices = p.slices > 0 ? p.slices : defaultSlices;
    // defaults as in 1D: FFF_rig for the smooth analytic shapes, FFF_sym for user functions, FFF_A
    // for anisotropic media (N from eq:N2d is continuous and signed)
    fac_ = fac;
    if (fac_ == Factorization::Auto)
      fac_ = aniso_ ? Factorization::FFF_A
                    : (p.shape == Profile2DLayer::Shape::Function ? Factorization::FFF_sym : Factorization::FFF_rig);
    FG_REQUIRE(aniso_ ? (fac_ == Factorization::Laurent || fac_ == Factorization::FFF_A || fac_ == Factorization::FFF_B)
                      : (fac_ == Factorization::Laurent || fac_ == Factorization::FFF_iso ||
                         fac_ == Factorization::FFF_sym || fac_ == Factorization::FFF_rig ||
                         fac_ == Factorization::FFF_A || fac_ == Factorization::FFF_B),
               "profile2d: factorization not applicable");
    tensorPath_ = aniso_ || fac_ == Factorization::FFF_A || fac_ == Factorization::FFF_B;
    warnFFFiso(fac_, {&mb, &ma}, warnings);
    nf_ = normalField(p_, dx_, dz_, map);  // isotropic FFF and the FFF-consistent reconstruction
    const int Pz = 2 * map.Nz();
    Mz_ = p.zSamples > 0 ? p.zSamples : std::max(128, 16 * (Pz + 1));
    FG_REQUIRE(Mz_ > 2 * Pz, "profile2d: z_samples must exceed 4 Nz");
  }
  double thickness_um() const override { return p_.depth_um; }
  bool yInvariant() const override { return false; }
  QBlocks Q(double y) const override {
    const double a = p_.depth_um;
    const double yc = std::clamp(y, 1e-14 * a, a * (1 - 1e-14));
    const int Px = 2 * map_.Nx(), Pz = 2 * map_.Nz();
    if (tensorPath_) {
      if (p_.shape == Profile2DLayer::Shape::SinusX) {
        // z-invariant: identical to the 1D anisotropic profile (V8/V9)
        std::vector<double> br;
        std::vector<bool> below;
        profileIntervals(line(0.0), yc, dx_, br, below);
        std::vector<Mat3C> tens;
        for (bool b : below) tens.push_back(b ? tb_ : ta_);
        return anisotropicQ(br, tens, [this](double x) { return profile2DNormal(p_, dx_, dz_, x, 0.0); }, fac_, dx_,
                            map_);
      }
      std::vector<Grid> g(static_cast<std::size_t>(anisoCount(fac_)), Grid(Px, Pz));
      std::vector<double> br;
      std::vector<bool> below;
      zIntegrate(yc, [&](double z, double w) {
        lineIntervals(z, yc, br, below);
        std::vector<Mat3C> tens;
        for (bool b : below) tens.push_back(b ? tb_ : ta_);
        const auto c = anisoLineCoefficients(br, tens, [&](double x) { return profile2DNormal(p_, dx_, dz_, x, z); },
                                             fac_, dx_, Px);
        accumulate(g, c, z, w);
      });
      return anisoQFromGrids(fac_, g, map_);
    }
    SliceCoefficients sc;
    if (p_.shape == Profile2DLayer::Shape::SinusX) {
      // z-invariant: exact intervals, q = 0 only (identical to the 1D sinus, V9)
      std::vector<double> br;
      std::vector<bool> below;
      profileIntervals(line(0.0), yc, dx_, br, below);
      std::vector<cplx> ve;
      for (bool b : below) ve.push_back(b ? epsB_ : epsA_);
      isoGrids(br, ve, dx_, map_, sc);
    } else {
      std::vector<Grid> g(2, Grid(Px, Pz));
      std::vector<double> br;
      std::vector<bool> below;
      zIntegrate(yc, [&](double z, double w) {
        lineIntervals(z, yc, br, below);
        std::vector<cplx> ve, vi;
        for (bool b : below) {
          ve.push_back(b ? epsB_ : epsA_);
          vi.push_back(b ? 1.0 / epsB_ : 1.0 / epsA_);
        }
        accumulate(g, {fourier::piecewiseConstant1D(br, ve, dx_, Px), fourier::piecewiseConstant1D(br, vi, dx_, Px)},
                   z, w);
      });
      sc.eps = std::move(g[0]);
      sc.invEps = std::move(g[1]);
    }
    sc.normal = &nf_;
    return isotropicQ(fac_, sc, map_);
  }
  Mat3C epsAt(double x, double y, double z) const override {
    return profileValue(line(z), x, dx_) > y ? tb_ : ta_;
  }
  Vec3D normalAt(double x, double, double z) const override { return profile2DNormal(p_, dx_, dz_, x, z); }
  const NormalFieldCoefficients& NN(double) const override { return nf_; }

 private:
  /// The 1D profile x -> g(x, z) of the line z = const.
  Profile1DLayer line(double z) const {
    Profile1DLayer l;
    l.x0_um = p_.x0_um;
    switch (p_.shape) {
      case Profile2DLayer::Shape::SinusX:
        l.shape = Profile1DLayer::Shape::Sinus;
        l.depth_um = p_.depth_um;
        break;
      case Profile2DLayer::Shape::EggCrate:
        l.shape = Profile1DLayer::Shape::Sinus;
        l.depth_um = 0.5 * p_.depth_um * (1.0 + std::cos(2 * PI / dz_ * (z - p_.z0_um)));
        break;
      case Profile2DLayer::Shape::Function: {
        l.shape = Profile1DLayer::Shape::Function;
        l.depth_um = p_.depth_um;
        const auto g = p_.g;
        l.g = [g, z](double x) { return g(x, z); };
        if (p_.grad) {
          const auto gr = p_.grad;
          l.dg = [gr, z](double x) { return gr(x, z)[0]; };
        }
        break;
      }
    }
    return l;
  }
  /// Intervals (breaks, below flags) of the line z at height y.
  void lineIntervals(double z, double y, std::vector<double>& br, std::vector<bool>& below) const {
    const Profile1DLayer l = line(z);
    if (l.shape == Profile1DLayer::Shape::Sinus && l.depth_um <= y) {
      br = {0.0, dx_};
      below = {false};
    } else {
      profileIntervals(l, y, dx_, br, below);
    }
  }
  /// grids(p, q) += w e^{-i q Kz z} c[m](p)
  void accumulate(std::vector<Grid>& g, const std::vector<std::vector<cplx>>& c, double z, double w) const {
    const int Px = 2 * map_.Nx(), Pz = 2 * map_.Nz();
    const double Kz = 2 * PI / dz_;
    for (int q = -Pz; q <= Pz; ++q) {
      const cplx f = w * std::exp(-I * (Kz * q * z));
      for (std::size_t m = 0; m < c.size(); ++m)
        for (int p = -Px; p <= Px; ++p) g[m](p, q) += f * c[m][static_cast<std::size_t>(p + Px)];
    }
  }
  /// Quadrature along z for the line coefficients at height y (Plan section "Permittivitaet bei
  /// beliebiger Hoehe", Profile2D): the line coefficients are continuous in z with sqrt-type
  /// behaviour where the number of crossings changes (contour tangent to the line). Without such
  /// critical lines: trapezoid rule (periodic, smooth; exact for z-invariant profiles). Otherwise:
  /// Gauss-Legendre between the critical lines with z = z_a + L (1 - cos(pi t))/2, which removes the
  /// endpoint singularities. add(z, w) receives nodes and weights (sum of weights = 1).
  void zIntegrate(double y, const std::function<void(double, double)>& add) const {
    std::vector<double> br;
    std::vector<bool> below;
    const int Ms = 2 * Mz_;
    auto count = [&](double z) {
      lineIntervals(z, y, br, below);
      return br.size();
    };
    std::vector<double> zc;
    const std::size_t c0 = count(0.0);
    std::size_t cprev = c0;
    for (int k = 1; k <= Ms; ++k) {
      const double zb = dz_ * k / Ms;
      const std::size_t cb = k == Ms ? c0 : count(zb);
      if (cb != cprev) {
        double lo = dz_ * (k - 1) / Ms, hi = zb;
        for (int it = 0; it < 60 && hi - lo > 1e-15 * dz_; ++it) {
          const double mid = 0.5 * (lo + hi);
          (count(mid) == cprev ? lo : hi) = mid;
        }
        zc.push_back(0.5 * (lo + hi));
      }
      cprev = cb;
    }
    if (zc.empty()) {
      for (int k = 0; k < Mz_; ++k) add(dz_ * k / Mz_, 1.0 / Mz_);
      return;
    }
    std::sort(zc.begin(), zc.end());
    const int ng = 16;
    const auto [gx, gw] = gaussLegendre(ng);
    for (std::size_t j = 0; j < zc.size(); ++j) {
      const double za = zc[j], zb = j + 1 < zc.size() ? zc[j + 1] : zc[0] + dz_;
      const double L = zb - za;
      const int panels = std::max(1, static_cast<int>(std::ceil(static_cast<double>(Mz_) * L / dz_ / ng)));
      for (int m = 0; m < panels; ++m)
        for (int i = 0; i < ng; ++i) {
          const double t = (m + 0.5 * (1.0 + gx[static_cast<std::size_t>(i)])) / panels;
          const double z = za + 0.5 * L * (1.0 - std::cos(PI * t));
          const double w = 0.5 * gw[static_cast<std::size_t>(i)] / panels * 0.5 * L * PI * std::sin(PI * t) / dz_;
          add(std::fmod(z, dz_), w);
        }
    }
  }

  /// Gauss-Legendre nodes and weights on [-1, 1].
  static std::pair<std::vector<double>, std::vector<double>> gaussLegendre(int n) {
    std::vector<double> x(static_cast<std::size_t>(n)), w(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      double z = std::cos(PI * (i + 0.75) / (n + 0.5)), dp = 1;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1, p1 = z;
        for (int k = 2; k <= n; ++k) {
          const double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1);
        const double dz = p1 / dp;
        z -= dz;
        if (std::abs(dz) < 1e-16) break;
      }
      x[static_cast<std::size_t>(i)] = z;
      w[static_cast<std::size_t>(i)] = 2.0 / ((1 - z * z) * dp * dp);
    }
    return {x, w};
  }

  Profile2DLayer p_;
  double dx_, dz_;
  IndexMap map_;
  NormalFieldCoefficients nf_;
  bool aniso_ = false, tensorPath_ = false;
  cplx epsB_, epsA_;
  Mat3C tb_, ta_;
  Factorization fac_ = Factorization::Auto;
  int Mz_ = 128;
};

}  // namespace

std::unique_ptr<LayerModel> buildLayerModel(const Stack& stack, std::size_t k, const IndexMap& map, int defaultSlices) {
  FG_REQUIRE(k < stack.layers.size(), "layer index out of range");
  const LayerEntry& le = stack.layers[k];
  if (const auto* hl = std::get_if<HomogeneousLayer>(&le.spec))
    return std::make_unique<HomogeneousModel>(*hl, stack.material(hl->material), map);
  if (const auto* ll = std::get_if<Lamellar1DLayer>(&le.spec))
    return std::make_unique<LamellarModel>(stack, *ll, le.factorization, map);
  if (const auto* pl = std::get_if<Profile1DLayer>(&le.spec))
    return std::make_unique<ProfileModel>(stack, *pl, le.factorization, map, defaultSlices);
  if (const auto* gl = std::get_if<Polygon2DLayer>(&le.spec))
    return std::make_unique<CrossSectionModel>(stack, std::make_shared<PolygonGeometry>(*gl, stack.dx_um, stack.dz_um),
                                               le.factorization, map, "polygon2d");
  if (const auto* sl = std::get_if<SectionedLayer>(&le.spec))
    return std::make_unique<SectionedModel>(stack, *sl, le.factorization, map, defaultSlices);
  if (const auto* cl = std::get_if<Contour2DLayer>(&le.spec))
    return std::make_unique<CrossSectionModel>(
        stack,
        std::make_shared<ContourGeometry>(cl->thickness_um, cl->background, cl->regions, std::vector<BoundaryChain>{},
                                          cl->normalField, stack.dx_um, stack.dz_um, "contour2d"),
        le.factorization, map, "contour2d");
  if (const auto* ql = std::get_if<Profile2DLayer>(&le.spec))
    return std::make_unique<Profile2DModel>(stack, *ql, le.factorization, map, defaultSlices);
  throw NotImplemented("layer type");
}

LayerModels buildLayerModels(const Stack& stack, const IndexMap& map, int defaultSlices) {
  LayerModels out;
  for (std::size_t k = 0; k < stack.layers.size(); ++k) out.push_back(buildLayerModel(stack, k, map, defaultSlices));
  return out;
}

}  // namespace fg::detail
