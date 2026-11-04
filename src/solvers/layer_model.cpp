#include "layer_model.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/fourier/coefficients.hpp"
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

/// Q blocks of an anisotropic 1D slice (Plan eq:QaniA / eq:QaniB; Laurent for comparison).
QBlocks anisotropicQ(const std::vector<double>& br, const std::vector<Mat3C>& eps,
                     const std::function<Vec3D(double)>& Nfun, Factorization fac, double dx,
                     const IndexMap& map) {
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  SliceCoefficients sc;
  sc.isotropic = false;
  auto fill = [&](GridMat& g, const std::vector<std::vector<cplx>>& c, int offset) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        auto& grid = g[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        grid = Grid(Px, Pz);
        for (int p = -Px; p <= Px; ++p)
          grid(p, 0) = c[static_cast<std::size_t>(offset + 3 * i + j)][static_cast<std::size_t>(p + Px)];
      }
  };
  if (fac == Factorization::Laurent) {
    const auto c = fourier::piecewiseSmoothCoefficients(
        br, 9,
        [&](double, int k, cplx* out) {
          for (int q = 0; q < 9; ++q) out[q] = eps[static_cast<std::size_t>(k)](q / 3, q % 3);
        },
        dx, Px);
    fill(sc.epsTensor, c, 0);
    return qLaurent(sc, map);
  }
  if (fac == Factorization::FFF_A) {
    const auto c = fourier::piecewiseSmoothCoefficients(
        br, 18,
        [&](double x, int k, cplx* out) {
          const Vec3D N = Nfun(x);
          Vec3D T1, T2;
          profileTangents(N, T1, T2);
          const Mat3C& e = eps[static_cast<std::size_t>(k)];
          const Mat3C C = cMatrix(e, N, T1, T2);
          const Mat3C eC = e * C;
          for (int q = 0; q < 9; ++q) {
            out[q] = C(q / 3, q % 3);
            out[9 + q] = eC(q / 3, q % 3);
          }
        },
        dx, Px);
    fill(sc.C, c, 0);
    fill(sc.epsC, c, 9);
    return qFFFAnisotropicA(sc, map);
  }
  if (fac == Factorization::FFF_B) {
    const auto c = fourier::piecewiseSmoothCoefficients(
        br, 37,
        [&](double x, int k, cplx* out) {
          const ProjectorSet ps = projectorMatrices(eps[static_cast<std::size_t>(k)], Nfun(x));
          for (int q = 0; q < 9; ++q) {
            out[q] = ps.G1(q / 3, q % 3);
            out[9 + q] = ps.G2(q / 3, q % 3);
            out[18 + q] = ps.H1(q / 3, q % 3);
            out[27 + q] = ps.PT(q / 3, q % 3);
          }
          out[36] = ps.invXi;
        },
        dx, Px);
    fill(sc.G1, c, 0);
    fill(sc.G2, c, 9);
    fill(sc.H1, c, 18);
    fill(sc.PT, c, 27);
    sc.invXi = Grid(Px, Pz);
    for (int p = -Px; p <= Px; ++p) sc.invXi(p, 0) = c[36][static_cast<std::size_t>(p + Px)];
    return qFFFAnisotropicB(sc, map);
  }
  throw Error("factorization not applicable to an anisotropic layer (use FFF_A, FFF_B or Laurent)");
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
  Mat3C epsAt(double, double) const override { return mat_.tensor(); }
  Vec3D normalAt(double) const override { return Vec3D(0, 1, 0); }
  const NormalFieldCoefficients& NN() const override { return nn_; }

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
  Mat3C epsAt(double x, double) const override {
    const double xr = std::fmod(std::fmod(x, dx_) + dx_, dx_);
    std::size_t k = 0;
    while (k + 1 < tensors_.size() && xr >= l_.breakpoints_um[k + 1]) ++k;
    return tensors_[k];
  }
  Vec3D normalAt(double) const override { return Vec3D(1, 0, 0); }
  const NormalFieldCoefficients& NN() const override { return nn_; }

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
    // z-invariant profile in a 2D stack: exact (q = 0 coefficients only); anisotropic 2D: M8
    FG_REQUIRE(!st.is2D() || !aniso_, "anisotropic layers in 2D stacks: milestone M8");
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
      return anisotropicQ(br, tens, [this](double x) { return normalAt(x); }, fac_, dx_, map_);
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
  Mat3C epsAt(double x, double y) const override { return profileValue(p_, x, dx_) > y ? tb_ : ta_; }
  Vec3D normalAt(double x) const override {
    const double sl = profileSlope(p_, x, dx_);
    return Vec3D(-sl, 1.0, 0.0) / std::sqrt(1.0 + sl * sl);
  }
  const NormalFieldCoefficients& NN() const override { return nf_; }

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
const Material& isotropic2D(const Stack& st, const std::string& name) {
  const Material& m = st.material(name);
  if (!m.isotropic) throw NotImplemented("anisotropic material '" + name + "' in a 2D layer (milestone M8)");
  return m;
}

/// Weak-metal warning of the book variant (Plan section "Varianten der isotropen FFF").
void warnFFFiso(Factorization f, const std::vector<const Material*>& mats, std::vector<std::string>& warnings) {
  if (f != Factorization::FFF_iso) return;
  for (const Material* m : mats)
    if (m->eps.real() < 0 && std::abs(m->eps.imag()) < 0.5 * std::abs(m->eps.real()))
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

// ------------------------------------------------------------------ polygon cross-section (2D)
class Polygon2DModel final : public LayerModel {
 public:
  Polygon2DModel(const Stack& st, const Polygon2DLayer& l, Factorization fac, const IndexMap& map)
      : l_(l), dx_(st.dx_um), dz_(st.dz_um) {
    FG_REQUIRE(st.is2D(), "polygon2d layers need a 2D stack (period_um.z > 0)");
    std::vector<const Material*> mats{&isotropic2D(st, l.background)};
    for (auto& inc : l_.inclusions) {
      FG_REQUIRE(inc.vertices.size() >= 3, "polygon2d: inclusion needs at least 3 vertices");
      mats.push_back(&isotropic2D(st, inc.material));
      // counter-clockwise orientation (Plan eq:polygon)
      double a = 0;
      for (std::size_t j = 0; j < inc.vertices.size(); ++j) {
        const auto& p0 = inc.vertices[j];
        const auto& p1 = inc.vertices[(j + 1) % inc.vertices.size()];
        a += p0[0] * p1[1] - p1[0] * p0[1];
      }
      if (a < 0) std::reverse(inc.vertices.begin(), inc.vertices.end());
    }
    for (const Material* m : mats) {
      epsList.push_back(m->eps);
      tensors_.push_back(m->tensor());
    }
    // exact coefficients: eps = eps_bg + sum_k (eps_k - eps_bg) chi_k (Plan eq:polygon)
    const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
    SliceCoefficients sc;
    sc.eps = Grid(Px, Pz);
    sc.invEps = Grid(Px, Pz);
    const cplx eb = mats[0]->eps;
    sc.eps(0, 0) = eb;
    sc.invEps(0, 0) = 1.0 / eb;
    for (std::size_t k = 0; k < l_.inclusions.size(); ++k) {
      const cplx ek = mats[k + 1]->eps;
      if (ek == eb) continue;
      for (int p = -Px; p <= Px; ++p)
        for (int q = -Pz; q <= Pz; ++q) {
          const cplx c = fourier::polygonCoefficient(l_.inclusions[k].vertices, dx_, dz_, p, q);
          sc.eps(p, q) += (ek - eb) * c;
          sc.invEps(p, q) += (1.0 / ek - 1.0 / eb) * c;
        }
    }
    // default: FFF_sym (sign-free normal field; Plan section "Varianten der isotropen FFF", M7)
    fac_ = fac == Factorization::Auto ? Factorization::FFF_sym : fac;
    FG_REQUIRE(fac_ == Factorization::Laurent || fac_ == Factorization::FFF_iso || fac_ == Factorization::FFF_sym,
               "polygon2d: factorization must be laurent, FFF_iso or FFF_sym (FFF_rig needs a continuous "
               "signed normal field, which does not exist for general cross-sections)");
    warnFFFiso(fac_, mats, warnings);
    if (fac_ != Factorization::Laurent) nf_ = normalField(l_, dx_, dz_, map);
    sc.normal = &nf_;
    Q_ = isotropicQ(fac_, sc, map);
  }
  double thickness_um() const override { return l_.thickness_um; }
  QBlocks Q(double) const override { return Q_; }
  /// Permittivity in the plane z = 0 (vertical cut; full (x, z) dependence: M9).
  Mat3C epsAt(double x, double) const override {
    return tensors_[static_cast<std::size_t>(polygonRegion(l_, dx_, dz_, x, 0.0) + 1)];
  }
  Vec3D normalAt(double) const override {
    throw NotImplemented("real-space normal field of 2D cross-sections (near field 2D, milestone M9)");
  }
  const NormalFieldCoefficients& NN() const override { return nf_; }

 private:
  Polygon2DLayer l_;
  double dx_, dz_;
  std::vector<Mat3C> tensors_;  ///< background, inclusions
  Factorization fac_ = Factorization::Auto;
  NormalFieldCoefficients nf_;
  QBlocks Q_;
};

// ------------------------------------------------------------------ profile y = g(x, z) (2D)
class Profile2DModel final : public LayerModel {
 public:
  Profile2DModel(const Stack& st, const Profile2DLayer& p, Factorization fac, const IndexMap& map, int defaultSlices)
      : p_(p), dx_(st.dx_um), dz_(st.dz_um), map_(map) {
    FG_REQUIRE(st.is2D(), "profile2d layers need a 2D stack (period_um.z > 0)");
    FG_REQUIRE(p.depth_um > 0, "profile2d: depth must be positive");
    FG_REQUIRE(p.shape != Profile2DLayer::Shape::Function || static_cast<bool>(p.g), "profile2d: function missing");
    const Material& mb = isotropic2D(st, p.below);
    const Material& ma = isotropic2D(st, p.above);
    epsB_ = mb.eps;
    epsA_ = ma.eps;
    epsList = {epsB_, epsA_};
    rcwaSlices = p.slices > 0 ? p.slices : defaultSlices;
    // defaults as in 1D: FFF_rig for the smooth analytic shapes, FFF_sym for user functions
    fac_ = fac == Factorization::Auto
               ? (p.shape == Profile2DLayer::Shape::Function ? Factorization::FFF_sym : Factorization::FFF_rig)
               : fac;
    FG_REQUIRE(fac_ == Factorization::Laurent || fac_ == Factorization::FFF_iso || fac_ == Factorization::FFF_sym ||
                   fac_ == Factorization::FFF_rig,
               "profile2d: factorization must be laurent, FFF_iso, FFF_sym or FFF_rig");
    warnFFFiso(fac_, {&mb, &ma}, warnings);
    nf_ = normalField(p_, dx_, dz_, map);
    const int Pz = 2 * map.Nz();
    Mz_ = p.zSamples > 0 ? p.zSamples : std::max(128, 16 * (Pz + 1));
    FG_REQUIRE(Mz_ > 2 * Pz, "profile2d: z_samples must exceed 4 Nz");
  }
  double thickness_um() const override { return p_.depth_um; }
  bool yInvariant() const override { return false; }
  QBlocks Q(double y) const override {
    const double a = p_.depth_um;
    const double yc = std::clamp(y, 1e-14 * a, a * (1 - 1e-14));
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
      sc.eps = Grid(2 * map_.Nx(), 2 * map_.Nz());
      sc.invEps = sc.eps;
      zLineCoefficients(yc, sc.eps, sc.invEps);
    }
    sc.normal = &nf_;
    return isotropicQ(fac_, sc, map_);
  }
  /// Permittivity in the plane z = 0 (vertical cut; full (x, z) dependence: M9).
  Mat3C epsAt(double x, double y) const override {
    return (profileValue(line(0.0), x, dx_) > y ? epsB_ : epsA_) * Mat3C::Identity();
  }
  Vec3D normalAt(double) const override {
    throw NotImplemented("real-space normal field of 2D profiles (near field 2D, milestone M9)");
  }
  const NormalFieldCoefficients& NN() const override { return nf_; }

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
  /// eps and 1/eps coefficients at height y (Plan section "Permittivitaet bei beliebiger Hoehe",
  /// Profile2D): exact interval coefficients along x on lines z = const; along z the line
  /// coefficients are continuous with sqrt-type behaviour where the number of crossings changes
  /// (contour tangent to the line). Without such critical lines: trapezoid rule (periodic, smooth;
  /// exact for z-invariant profiles). Otherwise: Gauss-Legendre on the sub-intervals between the
  /// critical lines with z = z_a + L (1 - cos(pi t))/2, which removes the endpoint singularities.
  void zLineCoefficients(double y, Grid& eps, Grid& inv) const {
    const int Px = 2 * map_.Nx(), Pz = 2 * map_.Nz();
    const double Kz = 2 * PI / dz_;
    std::vector<double> br;
    std::vector<bool> below;
    auto add = [&](double z, double w) {
      lineIntervals(z, y, br, below);
      std::vector<cplx> ve, vi;
      for (bool b : below) {
        ve.push_back(b ? epsB_ : epsA_);
        vi.push_back(b ? 1.0 / epsB_ : 1.0 / epsA_);
      }
      const auto ce = fourier::piecewiseConstant1D(br, ve, dx_, Px);
      const auto ci = fourier::piecewiseConstant1D(br, vi, dx_, Px);
      for (int q = -Pz; q <= Pz; ++q) {
        const cplx f = w * std::exp(-I * (Kz * q * z));
        for (int p = -Px; p <= Px; ++p) {
          eps(p, q) += f * ce[static_cast<std::size_t>(p + Px)];
          inv(p, q) += f * ci[static_cast<std::size_t>(p + Px)];
        }
      }
    };
    // critical lines: change of the number of intervals between neighbouring sample lines
    const int Ms = 2 * Mz_;
    auto count = [&](double z) {
      lineIntervals(z, y, br, below);
      return br.size();
    };
    std::vector<double> zc;
    std::size_t c0 = count(0.0), cprev = c0;
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
  cplx epsB_, epsA_;
  Factorization fac_ = Factorization::Auto;
  int Mz_ = 128;
};

}  // namespace

std::vector<std::unique_ptr<LayerModel>> buildLayerModels(const Stack& stack, const IndexMap& map,
                                                          int defaultSlices) {
  std::vector<std::unique_ptr<LayerModel>> out;
  for (const LayerEntry& le : stack.layers) {
    if (const auto* hl = std::get_if<HomogeneousLayer>(&le.spec))
      out.push_back(std::make_unique<HomogeneousModel>(*hl, stack.material(hl->material), map));
    else if (const auto* ll = std::get_if<Lamellar1DLayer>(&le.spec))
      out.push_back(std::make_unique<LamellarModel>(stack, *ll, le.factorization, map));
    else if (const auto* pl = std::get_if<Profile1DLayer>(&le.spec))
      out.push_back(std::make_unique<ProfileModel>(stack, *pl, le.factorization, map, defaultSlices));
    else if (const auto* gl = std::get_if<Polygon2DLayer>(&le.spec))
      out.push_back(std::make_unique<Polygon2DModel>(stack, *gl, le.factorization, map));
    else if (const auto* ql = std::get_if<Profile2DLayer>(&le.spec))
      out.push_back(std::make_unique<Profile2DModel>(stack, *ql, le.factorization, map, defaultSlices));
    else
      throw NotImplemented("layer type");
  }
  return out;
}

}  // namespace fg::detail
