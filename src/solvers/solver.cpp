// FourierGrating - stack assembly shared by the RK and RCWA solvers (Plan sections 7, 9, 10).
#include "fg/solvers/solver.hpp"

#include <algorithm>
#include <cmath>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/fourier/coefficients.hpp"
#include "fg/physics/efficiency.hpp"
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
    if (!sub.isotropic) throw NotImplemented("anisotropic substrate (milestone M5)");
    const ModeBasis subModes = isotropicModes(s.epsSub, s.alpha, s.gamma);
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
    const bool subLossless = std::abs(s.epsSub.imag()) < 1e-14 && s.epsSub.real() > 0;
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      const auto ro = planeWaveOrder(s.epsSup, s.alpha(r), s.gamma(r), s.ref.beta(r), +1, b_top(r),
                                     b_top(nf + r));
      if (ro.propagating) {
        res.reflected.push_back({n, m, ro.Sy / Sinc, ro.amp_s, ro.amp_p});
        sumR += ro.Sy / Sinc;
      }
      const auto to = planeWaveOrder(s.epsSub, s.alpha(r), s.gamma(r), subModes.beta(r), -1,
                                     a_sub(r), a_sub(nf + r));
      if (subLossless && to.propagating) {
        res.transmitted.push_back({n, m, -to.Sy / Sinc, to.amp_s, to.amp_p});
        sumT += -to.Sy / Sinc;
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
      if (!mat.isotropic) throw NotImplemented("anisotropic homogeneous layer (milestone M5)");
      const ModeBasis mb = isotropicModes(mat.eps, s.alpha, s.gamma);
      return stableLayerSMatrix(layerModesHomogeneous(mb), hl->thickness_um * s.k0, s.ref);
    }
    if (const auto* ll = std::get_if<Lamellar1DLayer>(&le.spec)) return lamellar(stack, *ll, le, s, diag);
    throw NotImplemented("layer type (profile layers: M3, 2D layers: M7)");
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
    if (!iso) throw NotImplemented("anisotropic lamellar layer (milestone M5)");
    const int Px = 2 * s.map.Nx(), Pz = 2 * s.map.Nz();
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
      case Factorization::FFF_iso: Q = qLiLamellar(sc, s.map); break;  // FFF with N = x
      default: throw Error("factorization not applicable to an isotropic lamellar layer");
    }
    const MatC M = buildPropagationMatrix(Q, s.alpha, s.gamma);
    const double h = ll.thickness_um * s.k0;
    if (opt_.kind == SolverKind::RCWA) {
      const LayerModes lm = isBlockAntiDiagonal(M) ? layerModesReduced(M) : layerModesGeneral(M);
      return stableLayerSMatrix(lm, h, s.ref);
    }
    // RK: split into sub-slices with rho * h_sub <= ln(growthLimit) (Plan eq:slicecrit)
    const double rho = rhoEstimate(s, eps);
    const int nsub = std::max(1, static_cast<int>(std::ceil(rho * h / std::log(opt_.growthLimit))));
    RKStats st;
    const MatC TF = rkTransferConstant(M, h / nsub, ButcherTableau::byName(opt_.rkTableau),
                                       opt_.rkTol, rho, st);
    diag.rkSteps += st.steps;
    const SMatrix Ssub = rkSliceSMatrix(TF, s.ref);
    SMatrix S = Ssub;
    for (int k = 1; k < nsub; ++k) S = star(S, Ssub);
    return S;
  }

  SolverOptions opt_;
};

}  // namespace

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt) {
  return std::make_unique<StackSolver>(opt);
}

}  // namespace fg
