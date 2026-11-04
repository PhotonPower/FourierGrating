// FourierGrating - stack assembly shared by the RK and RCWA solvers (Plan sections 7, 9, 10).
#include "fg/solvers/solver.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/geometry/profile.hpp"
#include "fg/physics/efficiency.hpp"
#include "fg/physics/anisotropic.hpp"
#include "fg/physics/modes.hpp"
#include "fg/physics/qmatrix.hpp"
#include "fg/solvers/rcwa.hpp"
#include "fg/solvers/rk.hpp"

namespace fg {

namespace {

/// Normalized problem data for one wavelength.
struct Setup {
  IndexMap map{0, 0};
  double k0 = 0;
  VecC alpha, gamma;
  IncidentWave inc;
  cplx epsSup = 1.0, epsSub = 1.0;
  ModeBasis ref;  // reference medium = superstrate
};

/// Spectral radius estimate rho_max = max_r,j |sqrt(alpha_r^2 + gamma_r^2 - eps_j)| (Plan eq:slicecrit).
double rhoEstimate(const Setup& s, const std::vector<cplx>& eps) {
  double rho = 0;
  for (Eigen::Index r = 0; r < s.alpha.size(); ++r)
    for (const cplx e : eps)
      rho = std::max(rho, std::abs(std::sqrt(s.alpha(r) * s.alpha(r) + s.gamma(r) * s.gamma(r) - e)));
  return rho;
}

double spectralRadius(const MatC& M) {
  Eigen::ComplexEigenSolver<MatC> es(M, false);
  return es.eigenvalues().cwiseAbs().maxCoeff();
}

class StackSolver final : public Solver {
 public:
  explicit StackSolver(SolverOptions o) : opt_(std::move(o)) {}

  Result solve(const Stack& stack, const Wavelength& wl) override {
    const Setup s = setup(stack, wl);
    const int nf = s.map.size();
    const int nh = 2 * nf;
    Result res;

    // substrate -> reference interface (Plan eq:Stot)
    const Material& sub = stack.material(stack.substrate);
    const ModeBasis subModes = sub.isotropic ? isotropicModes(s.epsSub, s.alpha, s.gamma)
                                             : anisotropicModes(sub.epsTensor, s.alpha, s.gamma);
    SMatrix S = fromTransfer(s.ref.PsiInv * subModes.Psi);

    for (const LayerEntry& le : stack.layers) S = star(S, layerSMatrix(stack, le, s, res.diag));

    // incident wave: tangential E of the (0,0) down-going order at the top interface
    VecC a_in = VecC::Zero(nh);
    const int c = s.map.center();
    a_in(c) = s.inc.A(0);
    a_in(nf + c) = s.inc.A(2);
    const VecC b_top = S.S12 * a_in;
    const VecC a_sub = S.S22 * a_in;

    const auto incOrder = planeWaveOrder(s.epsSup, s.alpha(c), s.gamma(c), s.ref.beta(c), -1,
                                         a_in(c), a_in(nf + c));
    const double Sinc = std::abs(incOrder.Sy);
    double sumR = 0, sumT = 0;
    const bool subLossless =
        sub.isotropic ? (std::abs(s.epsSub.imag()) < 1e-14 && s.epsSub.real() > 0)
                      : (sub.epsTensor - sub.epsTensor.adjoint()).cwiseAbs().maxCoeff() < 1e-14;
    const VecC Fsub = subModes.Psi.leftCols(nh) * a_sub;  // field just below y = 0
    auto split = [](OrderEfficiency& o, double e) {
      o.efficiency = e;
      const double ws = std::norm(o.amp_s), wp = std::norm(o.amp_p);
      if (ws + wp > 0) {
        o.eff_s = e * ws / (ws + wp);
        o.eff_p = e * wp / (ws + wp);
      }
    };
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      const auto ro = planeWaveOrder(s.epsSup, s.alpha(r), s.gamma(r), s.ref.beta(r), +1, b_top(r),
                                     b_top(nf + r));
      if (ro.propagating) {
        OrderEfficiency o{n, m, 0, ro.amp_s, ro.amp_p};
        split(o, ro.Sy / Sinc);
        res.reflected.push_back(o);
        sumR += ro.Sy / Sinc;
      }
      if (!subLossless) continue;
      if (sub.isotropic) {
        const auto to = planeWaveOrder(s.epsSub, s.alpha(r), s.gamma(r), subModes.beta(r), -1,
                                       a_sub(r), a_sub(nf + r));
        if (to.propagating) {
          OrderEfficiency o{n, m, 0, to.amp_s, to.amp_p};
          split(o, -to.Sy / Sinc);
          res.transmitted.push_back(o);
          sumT += -to.Sy / Sinc;
        }
      } else {
        // anisotropic substrate: flux of the order from the modal field (Plan eq:Sy)
        const double sy = 0.5 * std::real(Fsub(nf + r) * std::conj(Fsub(2 * nf + r)) -
                                          Fsub(r) * std::conj(Fsub(3 * nf + r)));
        if (std::abs(sy) > 1e-14 * Sinc) {
          OrderEfficiency o{n, m, -sy / Sinc, 0, 0};
          res.transmitted.push_back(o);
          sumT += -sy / Sinc;
        }
      }
    }
    res.absorption = 1.0 - sumR - sumT;
    res.diag.energyBalance = sumR + sumT;
    return res;
  }

 private:
  Setup setup(const Stack& stack, const Wavelength& wl) const {
    Setup s;
    s.map = IndexMap(opt_.Nx, stack.is2D() ? opt_.Nz : 0);
    s.k0 = 2.0 * PI / wl.lambda_um;
    const Material& sup = stack.material(stack.superstrate);
    FG_REQUIRE(sup.isotropic && std::abs(sup.eps.imag()) < 1e-14 && sup.eps.real() > 0,
               "superstrate must be isotropic and lossless");
    s.epsSup = sup.eps;
    s.epsSub = stack.material(stack.substrate).eps;
    s.inc = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
    const double Kx = 2.0 * PI / (s.k0 * stack.dx_um);
    const double Kz = stack.is2D() ? 2.0 * PI / (s.k0 * stack.dz_um) : 0.0;
    const int nf = s.map.size();
    s.alpha.resize(nf);
    s.gamma.resize(nf);
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      s.alpha(r) = s.inc.alpha0 + n * Kx;
      s.gamma(r) = s.inc.gamma0 + m * Kz;
    }
    s.ref = isotropicModes(s.epsSup, s.alpha, s.gamma);
    return s;
  }

  SMatrix layerSMatrix(const Stack& stack, const LayerEntry& le, const Setup& s, Diagnostics& diag) {
    if (const auto* hl = std::get_if<HomogeneousLayer>(&le.spec)) {
      const Material& mat = stack.material(hl->material);
      const ModeBasis mb = mat.isotropic ? isotropicModes(mat.eps, s.alpha, s.gamma)
                                         : anisotropicModes(mat.epsTensor, s.alpha, s.gamma);
      return stableLayerSMatrix(layerModesHomogeneous(mb), hl->thickness_um * s.k0, s.ref);
    }
    if (const auto* ll = std::get_if<Lamellar1DLayer>(&le.spec)) return lamellar(stack, *ll, le, s, diag);
    if (const auto* pl = std::get_if<Profile1DLayer>(&le.spec)) return profile(stack, *pl, le, s, diag);
    throw NotImplemented("layer type (2D layers: milestone M7)");
  }

  SMatrix lamellar(const Stack& stack, const Lamellar1DLayer& ll, const LayerEntry& le,
                   const Setup& s, Diagnostics& diag) {
    FG_REQUIRE(ll.breakpoints_um.size() == ll.materials.size() + 1, "lamellar: breakpoints/materials");
    std::vector<cplx> eps, inv;
    bool iso = true;
    for (const auto& name : ll.materials) {
      const Material& m = stack.material(name);
      iso = iso && m.isotropic;
      eps.push_back(m.eps);
      inv.push_back(1.0 / m.eps);
    }
    const int Px = 2 * s.map.Nx(), Pz = 2 * s.map.Nz();
    if (!iso) {
      std::vector<Mat3C> tens;
      for (const auto& name : ll.materials) tens.push_back(stack.material(name).tensor());
      const Factorization fac = le.factorization == Factorization::Auto ? Factorization::FFF_A : le.factorization;
      const QBlocks Q = anisotropicQ(ll.breakpoints_um, tens, [](double) { return Vec3D(1, 0, 0); }, fac,
                                     stack.dx_um, s.map);
      return constantSlice(buildPropagationMatrix(Q, s.alpha, s.gamma), ll.thickness_um * s.k0, eps, s, diag);
    }
    const auto ce = fourier::piecewiseConstant1D(ll.breakpoints_um, eps, stack.dx_um, Px);
    const auto ci = fourier::piecewiseConstant1D(ll.breakpoints_um, inv, stack.dx_um, Px);
    SliceCoefficients sc;
    sc.eps = fourier::CoefficientGrid(Px, Pz);
    sc.invEps = fourier::CoefficientGrid(Px, Pz);
    for (int p = -Px; p <= Px; ++p) {
      sc.eps(p, 0) = ce[static_cast<std::size_t>(p + Px)];
      sc.invEps(p, 0) = ci[static_cast<std::size_t>(p + Px)];
    }
    QBlocks Q;
    switch (le.factorization) {
      case Factorization::Laurent: Q = qLaurent(sc, s.map); break;
      case Factorization::Auto:
      case Factorization::LiLamellar:
      case Factorization::FFF_iso:
      case Factorization::FFF_sym:
      case Factorization::FFF_rig: Q = qLiLamellar(sc, s.map); break;  // FFF with N = x (exact: NN = diag)
      default: throw Error("factorization not applicable to an isotropic lamellar layer");
    }
    return constantSlice(buildPropagationMatrix(Q, s.alpha, s.gamma), ll.thickness_um * s.k0, eps, s, diag);
  }

  /// y-invariant slice of thickness h with constant M: RCWA modes, or RK sub-slices
  /// (Plan eq:slicecrit with the spectral radius of M).
  SMatrix constantSlice(const MatC& M, double h, const std::vector<cplx>& eps, const Setup& s, Diagnostics& diag) {
    if (opt_.kind == SolverKind::RCWA) return stableLayerSMatrix(layerModesAuto(M), h, s.ref);
    const double rho = std::max(rhoEstimate(s, eps), spectralRadius(M));
    const int nsub = std::max(1, static_cast<int>(std::ceil(rho * h / std::log(opt_.growthLimit))));
    RKStats st;
    const MatC TF = rkTransferConstant(M, h / nsub, ButcherTableau::byName(opt_.rkTableau), opt_.rkTol, rho, st);
    diag.rkSteps += st.steps;
    const SMatrix Ssub = rkSliceSMatrix(TF, s.ref);
    SMatrix S = Ssub;
    for (int k = 1; k < nsub; ++k) S = star(S, Ssub);
    return S;
  }

  /// Q blocks of an anisotropic slice: piecewise constant tensors on [br_k, br_k+1] and a normal
  /// field N(x) (Plan eq:QaniA / eq:QaniB; Laurent for comparison).
  static QBlocks anisotropicQ(const std::vector<double>& br, const std::vector<Mat3C>& eps,
                              const std::function<Vec3D(double)>& Nfun, Factorization fac, double dx,
                              const IndexMap& map) {
    const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
    SliceCoefficients sc;
    sc.isotropic = false;
    auto fill = [&](std::array<std::array<fourier::CoefficientGrid, 3>, 3>& g,
                    const std::vector<std::vector<cplx>>& c, int offset) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          auto& grid = g[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
          grid = fourier::CoefficientGrid(Px, Pz);
          for (int p = -Px; p <= Px; ++p)
            grid(p, 0) = c[static_cast<std::size_t>(offset + 3 * i + j)][static_cast<std::size_t>(p + Px)];
        }
    };
    if (fac == Factorization::Laurent) {
      const auto c = fourier::piecewiseSmoothCoefficients(
          br, 9, [&](double, int k, cplx* out) { for (int q = 0; q < 9; ++q) out[q] = eps[static_cast<std::size_t>(k)](q / 3, q % 3); },
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
      sc.invXi = fourier::CoefficientGrid(Px, Pz);
      for (int p = -Px; p <= Px; ++p) sc.invXi(p, 0) = c[36][static_cast<std::size_t>(p + Px)];
      return qFFFAnisotropicB(sc, map);
    }
    throw Error("factorization not applicable to an anisotropic layer (use FFF_A, FFF_B or Laurent)");
  }

  /// RK slice S-matrix over [y0, y1]; bisects if the transfer matrix grows beyond growthLimit.
  SMatrix rkSlice(const std::function<MatC(double)>& M, double y0, double y1, double rho,
                  const ButcherTableau& tab, const ModeBasis& ref, Diagnostics& diag, int depth) {
    RKStats st;
    const MatC TF = rkTransferAuto(M, y0, y1, tab, opt_.rkTol, rho, st);
    diag.rkSteps += st.steps;
    const double growth = TF.cwiseAbs().maxCoeff();
    if ((!TF.allFinite() || growth > opt_.growthLimit) && depth < 12) {
      const double ym = 0.5 * (y0 + y1);
      return star(rkSlice(M, y0, ym, rho, tab, ref, diag, depth + 1),
                  rkSlice(M, ym, y1, rho, tab, ref, diag, depth + 1));
    }
    FG_REQUIRE(TF.allFinite(), "RK transfer matrix not finite");
    return rkSliceSMatrix(TF, ref);
  }

  /// Profile layer y = g(x) (Plan sections "Permittivitaet bei beliebiger Hoehe", Loeser A/B).
  SMatrix profile(const Stack& stack, const Profile1DLayer& pl, const LayerEntry& le, const Setup& s,
                  Diagnostics& diag) {
    FG_REQUIRE(!stack.is2D(), "Profile1D layers are only supported in 1D stacks (2D: Profile2D, M7)");
    const Material& mb = stack.material(pl.below);
    const Material& ma = stack.material(pl.above);
    const bool aniso = !mb.isotropic || !ma.isotropic;
    const int Px = 2 * s.map.Nx(), Pz = 2 * s.map.Nz();
    const NormalFieldCoefficients nf = normalField(pl, stack.dx_um, s.map);
    // Default (Plan section "Varianten der isotropen FFF"): FFF_rig for smooth profiles (signed
    // normal field continuous), FFF_sym for profiles with corners (piecewise linear, user functions).
    Factorization fac = le.factorization;
    if (fac == Factorization::Auto)
      fac = aniso ? Factorization::FFF_A
                  : (pl.shape == Profile1DLayer::Shape::Sinus ? Factorization::FFF_rig : Factorization::FFF_sym);
    if (fac == Factorization::FFF_iso)
      for (const Material* m : {&mb, &ma})
        if (m->eps.real() < 0 && std::abs(m->eps.imag()) < 0.5 * std::abs(m->eps.real()))
          diag.warnings.push_back("FFF_iso (book IV.32) with a weakly absorbing metal ('" + m->name +
                                  "'): known to produce unphysical results, use FFF_rig or FFF_sym");
    FG_REQUIRE(aniso ? (fac == Factorization::FFF_A || fac == Factorization::FFF_B || fac == Factorization::Laurent)
                     : (fac == Factorization::FFF_iso || fac == Factorization::FFF_sym ||
                        fac == Factorization::FFF_rig || fac == Factorization::Laurent),
               "factorization not applicable to this profile layer");
    const Mat3C tb = mb.tensor(), ta = ma.tensor();
    auto Nfun = [&](double x) -> Vec3D {
      const double sl = profileSlope(pl, x, stack.dx_um);
      return Vec3D(-sl, 1.0, 0.0) / std::sqrt(1.0 + sl * sl);
    };
    const double a = pl.depth_um;
    // M at height y_um inside the layer
    auto Mat = [&](double y_um) {
      std::vector<double> br;
      std::vector<bool> below;
      profileIntervals(pl, std::clamp(y_um, 1e-14 * a, a * (1 - 1e-14)), stack.dx_um, br, below);
      if (aniso) {
        std::vector<Mat3C> tens;
        for (bool b : below) tens.push_back(b ? tb : ta);
        return buildPropagationMatrix(anisotropicQ(br, tens, Nfun, fac, stack.dx_um, s.map), s.alpha, s.gamma);
      }
      std::vector<cplx> ve, vi;
      for (bool b : below) {
        const cplx e = b ? mb.eps : ma.eps;
        ve.push_back(e);
        vi.push_back(1.0 / e);
      }
      const auto ce = fourier::piecewiseConstant1D(br, ve, stack.dx_um, Px);
      const auto ci = fourier::piecewiseConstant1D(br, vi, stack.dx_um, Px);
      SliceCoefficients sc;
      sc.eps = fourier::CoefficientGrid(Px, Pz);
      sc.invEps = fourier::CoefficientGrid(Px, Pz);
      for (int p = -Px; p <= Px; ++p) {
        sc.eps(p, 0) = ce[static_cast<std::size_t>(p + Px)];
        sc.invEps(p, 0) = ci[static_cast<std::size_t>(p + Px)];
      }
      sc.normal = &nf;
      const QBlocks Q = fac == Factorization::Laurent  ? qLaurent(sc, s.map)
                        : fac == Factorization::FFF_rig ? qFFFIsotropicRigorous(sc, s.map)
                                                        : qFFFIsotropic(sc, s.map, fac == Factorization::FFF_sym);
      return buildPropagationMatrix(Q, s.alpha, s.gamma);
    };
    const double h = a * s.k0;  // normalized depth
    SMatrix S = SMatrix::identity(2 * s.map.size());
    if (opt_.kind == SolverKind::RCWA) {
      const int K = pl.slices > 0 ? pl.slices : opt_.rcwaSlices;
      for (int k = 0; k < K; ++k) {
        const MatC M = Mat((k + 0.5) * a / K);  // staircase: eps at slice centre, smooth N
        const LayerModes lm = layerModesAuto(M);
        S = star(S, stableLayerSMatrix(lm, h / K, s.ref));
      }
      return S;
    }
    // RK: stability sub-slices (Plan eq:slicecrit) with the spectral radius of the truncated M
    // (FFF can create eigenvalues far above the material estimate), plus a posteriori bisection
    // whenever the field transfer matrix grows beyond growthLimit.
    double rho = rhoEstimate(s, {mb.eps, ma.eps});
    for (double t : {0.02, 0.25, 0.5, 0.75, 0.98}) rho = std::max(rho, spectralRadius(Mat(t * a)));
    const int K = std::max(pl.slices > 0 ? pl.slices : 1,
                           static_cast<int>(std::ceil(rho * h / std::log(opt_.growthLimit))));
    const auto tab = ButcherTableau::byName(opt_.rkTableau);
    for (int k = 0; k < K; ++k)
      S = star(S, rkSlice([&](double yn) { return Mat(yn / s.k0); }, h * k / K, h * (k + 1) / K, rho,
                          tab, s.ref, diag, 0));
    return S;
  }

  SolverOptions opt_;
};

}  // namespace

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt) {
  return std::make_unique<StackSolver>(opt);
}

}  // namespace fg
