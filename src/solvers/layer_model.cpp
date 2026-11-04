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
    FG_REQUIRE(!st.is2D(), "Profile1D layers are only supported in 1D stacks (2D: Profile2D, M7)");
    const Material& mb = st.material(p.below);
    const Material& ma = st.material(p.above);
    aniso_ = !mb.isotropic || !ma.isotropic;
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
    else
      throw NotImplemented("layer type (2D layers: milestone M7)");
  }
  return out;
}

}  // namespace fg::detail
