// FourierGrating - stack assembly shared by the RK and RCWA solvers, efficiencies and near fields
// (Plan sections "S-Matrix-Formalismus", "Loeser A/B", "Nahfeld").
#include <algorithm>
#include <array>
#include <cstdlib>
#include <exception>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/core/parallel.hpp"
#include "fg/core/prof.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/fourier/toeplitz.hpp"
#include "fg/geometry/profile.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/physics/anisotropic.hpp"
#include "fg/physics/efficiency.hpp"
#include "fg/cyl/propagation.hpp"
#include "cyl_basis.hpp"
#include "fg/physics/modes.hpp"
#include "fg/solvers/rcwa.hpp"
#include "fg/solvers/resident.hpp"
#include "fg/solvers/rk.hpp"
#include "fg/solvers/solver.hpp"

#include "layer_model.hpp"

namespace fg {

namespace detail {

/// Applies an n_f x n_f matrix to each of the four field components of F.
VecC componentwise(const MatC& T, const VecC& F) {
  const Eigen::Index n = T.rows();
  VecC out(F.size());
  for (Eigen::Index k = 0; k < F.size() / n; ++k) out.segment(k * n, n) = T * F.segment(k * n, n);
  return out;
}

/// Modes of a homogeneous isotropic medium in transformed coordinates. There P R = (A^2 + G^2 - eps) (x) 1_2
/// (one n_f block per tangential E component): the TE/TM pairs are exactly degenerate, and the
/// eigenvectors of the full reduced matrix (Schur back substitution of a non-normal matrix) become
/// ill-conditioned for small eta. The n_f block is diagonalized once and both components take its
/// eigenvectors. Any other structure: layerModesAuto.
LayerModes homogeneousModes(const MatC& M) {
  if (!isBlockAntiDiagonal(M)) return layerModesAuto(M);
  const Eigen::Index n = M.rows() / 2, nf = n / 2;
  const MatC P = M.topRightCorner(n, n), R = M.bottomLeftCorner(n, n);
  const MatC PR = P * R;
  const MatC B = PR.topLeftCorner(nf, nf);
  const double tol = 1e-12 * std::max(1.0, PR.cwiseAbs().maxCoeff());
  if ((PR.bottomRightCorner(nf, nf) - B).cwiseAbs().maxCoeff() > tol ||
      PR.topRightCorner(nf, nf).cwiseAbs().maxCoeff() > tol || PR.bottomLeftCorner(nf, nf).cwiseAbs().maxCoeff() > tol)
    return layerModesAuto(M);
  Eigen::ComplexEigenSolver<MatC> es(B, true);
  FG_REQUIRE(es.info() == Eigen::Success, "eigen decomposition failed (transformed homogeneous medium)");
  LayerModes m;
  m.Wu.resize(2 * n, n);
  m.Wd.resize(2 * n, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index c = 0; c < 2; ++c)
    for (Eigen::Index k = 0; k < nf; ++k) {
      cplx kap = std::sqrt(-es.eigenvalues()(k));  // branch as in layerModesReduced
      if (kap.imag() < 0 || (kap.imag() == 0 && kap.real() < 0)) kap = -kap;
      if (std::abs(kap) < 1e-10) kap += cplx(0, 1e-10);
      VecC u = VecC::Zero(n);
      u.segment(c * nf, nf) = es.eigenvectors().col(k);
      const VecC Ru = R * u;
      const Eigen::Index j = c * nf + k;
      m.lu(j) = I * kap;
      m.ld(j) = -I * kap;
      m.Wu.col(j) << u, Ru / m.lu(j);
      m.Wd.col(j) << u, Ru / m.ld(j);
    }
  return m;
}

/// Mode basis of a homogeneous half space from the eigen-decomposition of M (Plan section "M17": in
/// transformed coordinates the orders are no eigenvectors): columns (down | up), classified by Re lambda
/// and, for |Re lambda| <= tau, by the flux (Plan eq:classify). The reduced path takes the sign of the
/// rounding noise Im kappa for propagating modes, which is harmless inside a layer but not here.
ModeBasis numericBasis(const MatC& M) {
  const LayerModes m = homogeneousModes(M);
  const Eigen::Index n = m.Wu.cols();
  MatC W(M.rows(), 2 * n);
  W << m.Wu, m.Wd;
  VecC lam(2 * n);
  lam << m.lu, m.ld;
  // propagating modes (|Re lambda| <= 1e-6 |lambda|, rounding noise of the non-Hermitian M with A, G; with
  // PML up to 1e-7 |lambda|): by the flux, also when the noise exceeds the global tolerance of modeOrder
  VecC key = lam;
  for (Eigen::Index k = 0; k < key.size(); ++k)
    if (std::abs(key(k).real()) <= 1e-6 * std::abs(key(k))) key(k) = cplx(0.0, key(k).imag());
  const std::vector<Eigen::Index> idx = modeOrder(key, [&](Eigen::Index k) { return modeFlux(W.col(k)); });
  ModeBasis b;
  b.Psi.resize(M.rows(), 2 * n);
  b.lambda.resize(2 * n);
  for (Eigen::Index k = 0; k < n; ++k) {
    const Eigen::Index up = idx[static_cast<std::size_t>(k)], dn = idx[static_cast<std::size_t>(n + k)];
    b.Psi.col(k) = W.col(dn);
    b.lambda(k) = lam(dn);
    b.Psi.col(n + k) = W.col(up);
    b.lambda(n + k) = lam(up);
  }
  b.PsiInv = b.Psi.partialPivLu().inverse();
  return b;
}

/// Separable n_f x n_f matrix T(r, r') = Tx(n, n') Tz(m, m') (identity blocks for null factors).
MatC separable(const IndexMap& map, const MatC* Tx, const MatC* Tz) {
  const int nf = map.size();
  MatC T(nf, nf);
  for (int r = 0; r < nf; ++r)
    for (int c = 0; c < nf; ++c) {
      const auto [n, m] = map.orders(r);
      const auto [n2, m2] = map.orders(c);
      const cplx tx = Tx ? (*Tx)(n + map.Nx(), n2 + map.Nx()) : cplx(n == n2 ? 1.0 : 0.0);
      const cplx tz = Tz ? (*Tz)(m + map.Nz(), m2 + map.Nz()) : cplx(m == m2 ? 1.0 : 0.0);
      T(r, c) = tx * tz;
    }
  return T;
}

TransformedCoordinates transformedCoordinates(const IndexMap& map, const ResolvedTransform& tr, const PmlMap* pmlX,
                                              const PmlMap* pmlZ, const VecC& alpha, const VecC& gamma, double kx0,
                                              double kz0) {
  const int nf = map.size(), Nx = map.Nx(), Nz = map.Nz();
  std::vector<cplx> cx, cz;
  if (tr.x) cx = tr.x->inverseDerivativeCoefficients(2 * Nx);
  if (tr.z) cz = tr.z->inverseDerivativeCoefficients(2 * Nz);
  if (pmlX) cx = pmlX->inverseDerivativeCoefficients(2 * Nx);  // Plan eq:pml
  if (pmlZ) cz = pmlZ->inverseDerivativeCoefficients(2 * Nz);
  TransformedCoordinates t;
  t.A = MatC::Zero(nf, nf);
  t.G = MatC::Zero(nf, nf);
  for (int r = 0; r < nf; ++r)
    for (int c = 0; c < nf; ++c) {
      const auto [n, m] = map.orders(r);
      const auto [n2, m2] = map.orders(c);
      // [[1/f']] is Toeplitz in n and the identity in m (and vice versa)
      if (m == m2) t.A(r, c) = (!cx.empty() ? cx[static_cast<std::size_t>(n - n2 + 2 * Nx)] : cplx(n == n2 ? 1.0 : 0.0)) * alpha(c);
      if (n == n2) t.G(r, c) = (!cz.empty() ? cz[static_cast<std::size_t>(m - m2 + 2 * Nz)] : cplx(m == m2 ? 1.0 : 0.0)) * gamma(c);
    }
  MatC Cx, Cz, Dx, Dz;
  if (tr.x) {
    Cx = tr.x->toPhysical(kx0, Nx);
    Dx = tr.x->fromPhysical(kx0, Nx);
  }
  if (tr.z) {
    Cz = tr.z->toPhysical(kz0, Nz);
    Dz = tr.z->fromPhysical(kz0, Nz);
  }
  t.toPhys = separable(map, tr.x ? &Cx : nullptr, tr.z ? &Cz : nullptr);
  t.fromPhys = separable(map, tr.x ? &Dx : nullptr, tr.z ? &Dz : nullptr);
  return t;
}

}  // namespace detail

namespace {

using detail::componentwise;
using detail::homogeneousModes;
using detail::LayerModel;
using detail::numericBasis;
using detail::separable;

/// Normalized problem data for one wavelength.
struct Setup {
  IndexMap map{0, 0};
  double k0 = 0;
  VecC alpha, gamma;
  IncidentWave inc;
  cplx epsSup = 1.0;
  ModeBasis ref;  // reference medium = superstrate
  const Material* sub = nullptr;
  ModeBasis subModes;
  // coordinate transformation (Plan section "M17"): A = [[1/f']] alpha, G = [[1/g']] gamma; ref and
  // subModes are then numerical bases in X, refPhys/subPhys the analytic bases of the physical orders
  ResolvedTransform tr;
  bool transformed = false;
  MatC A, G;
  MatC toPhys, fromPhys;  ///< n_f x n_f per field component: eq:ct_rayleigh, eq:ct_inc (separable)
  ModeBasis refPhys, subPhys;
  // aperiodic structure (Plan section "M18"): PML zones enter only through [[1/f']] (f(X) = X wherever
  // the structure is); the transformed path above with these coefficients
  bool aperiodic = false;
  std::shared_ptr<const PmlMap> pmlX, pmlZ;
  // structured half spaces (Plan "Einfall (b)"): numerical bases from their cross-sections
  std::optional<QBlocks> supQ, subQ;
  // scattered-field formulation (Plan section "Streufeldformulierung"): the background stack at the
  // incidence order only (Airy), per background layer as four exponential modes (normalized y)
  bool scattered = false;
  struct BgLayer {
    double y0 = 0, y1 = 0;
    MatC Psi;          ///< 4 x 4 single-order modes (down x, down z | up x, up z)
    VecC lambda, amp;  ///< amplitudes: down referenced at y1, up at y0
    MatC M, MB;        ///< propagation matrix of the background material (computation / Bloch basis)
    QBlocks Q, Q1;     ///< background permittivity (n_f orders / the incidence order)
  };
  std::vector<BgLayer> bg;
  VecC bgTop, bgSub;        ///< single order: up amplitudes at the top (superstrate), down at y = 0 (substrate)
  double incIntensity = 0;  ///< |S_y| of the incident plane wave
  /// Oblique incidence on a PML axis: the scattered field has no Bloch phase there; the sources are
  /// formed in the Bloch basis (alphaB, gammaB) and re-expanded over the physical region (Tsrc).
  bool blochSource = false;
  VecC alphaB, gammaB;
  MatC Tsrc;
  /// Cylindrical plane wave (M20): sigma of a layer for the background (Ex, Ez, Hx, Hz) of the
  /// incidence order (Plan eq:cyl_scat); replaces (M - M_bg) F_bg.
  std::function<VecC(const detail::StackLayer&, const BgLayer&, const VecC&)> cylSource;
  const ModeBasis& physRef() const { return transformed ? refPhys : ref; }
  const ModeBasis& physSub() const { return transformed ? subPhys : subModes; }
};

/// Cartesian propagation matrix (Plan eq:Mmu1a): full 4 n_f system or one polarization (2 n_f).
MatC propagation(const Setup& s, const std::optional<Polarization>& pol, const QBlocks& Q) {
  if (pol) return buildPropagationMatrix(Q, s.alpha, *pol);
  return s.transformed ? buildPropagationMatrix(Q, s.A, s.G) : buildPropagationMatrix(Q, s.alpha, s.gamma);
}

/// The layers of a Cartesian stack are LayerModels (only the Cartesian paths downcast).
const LayerModel& cartesian(const detail::StackLayer& l) { return static_cast<const LayerModel&>(l); }

/// Field system of a Cartesian assembly (Plan section "M20: Architektur"): the full 4 n_f problem, or one
/// polarization of the decoupled 1D in-plane problem (2 n_f, Plan eq:TE / eq:TM).
class CartesianBasis final : public detail::ModalBasis {
 public:
  explicit CartesianBasis(const Setup& s, std::optional<Polarization> pol = std::nullopt)
      : s_(s), pol_(pol), ref_(basis(s.ref)), sub_(basis(s.subModes)) {}
  const ModeBasis& ref() const override { return ref_; }
  const ModeBasis& sub() const override { return sub_; }
  MatC M(const detail::StackLayer& layer, double y_um) const override {
    return propagation(s_, pol_, cartesian(layer).Q(y_um));
  }
  LayerModes modes(const MatC& M) const override {
    return pol_ ? layerModesPolarization(M, *pol_) : layerModesAuto(M);
  }
  LayerModes homogeneousModes(const detail::StackLayer& layer) const override {
    const LayerModel& m = cartesian(layer);
    const Material& mat = m.material();
    // analytic modes; in transformed coordinates numerical modes of M (Plan section "M17")
    return s_.transformed ? detail::homogeneousModes(propagation(s_, pol_, m.Q(0.0)))
                          : layerModesHomogeneous(basis(mat.isotropic ? isotropicModes(mat.eps, s_.alpha, s_.gamma)
                                                                      : anisotropicModes(mat.epsTensor, s_.alpha, s_.gamma)));
  }
  double rhoEstimate(const detail::StackLayer& layer) const override;
  std::optional<Polarization> polarization() const override { return pol_; }

 private:
  ModeBasis basis(const ModeBasis& full) const { return pol_ ? polarizationBasis(full, *pol_) : full; }
  const Setup& s_;
  std::optional<Polarization> pol_;
  ModeBasis ref_, sub_;
};

/// Layer pointers of the assembly.
std::vector<const detail::StackLayer*> layerPtrs(const detail::LayerModels& models) {
  std::vector<const detail::StackLayer*> v;
  for (const auto& m : models) v.push_back(m.get());
  return v;
}

/// Materials referenced by a 1D layer; nullopt for 2D layer types.
std::optional<std::vector<std::string>> layerMaterials(const LayerSpec& spec) {
  if (const auto* h = std::get_if<HomogeneousLayer>(&spec)) return std::vector<std::string>{h->material};
  if (const auto* l = std::get_if<Lamellar1DLayer>(&spec)) return l->materials;
  if (const auto* p = std::get_if<Profile1DLayer>(&spec)) return std::vector<std::string>{p->below, p->above};
  return std::nullopt;
}

/// Eigen's own threads (blas level) and nested OpenMP levels (Plan section "M15", step 2). Both are
/// process-wide settings: changed only outside parallel regions.
void setEigenThreads(const ThreadPlan& p) {
#ifdef FG_HAVE_OPENMP
  if (omp_in_parallel()) return;
  if (Eigen::nbThreads() != p.blas) Eigen::setNbThreads(p.blas);
#if _OPENMP >= 200805  // OpenMP 3.0
  if (p.scan > 1 && p.slices > 1 && omp_get_max_active_levels() < 2) omp_set_max_active_levels(2);
#else  // MSVC /openmp (2.0)
  if (p.scan > 1 && p.slices > 1) omp_set_nested(1);
#endif
#else
  (void)p;
#endif
}

double rhoEstimate(const Setup& s, const std::vector<cplx>& eps) {
  // transformed coordinates: |A| up to |alpha| / min f' = |alpha| / (1 - eta)
  const double sx = s.tr.x ? 1.0 / (1.0 - s.tr.x->eta()) : 1.0, sz = s.tr.z ? 1.0 / (1.0 - s.tr.z->eta()) : 1.0;
  double rho = 0;
  for (Eigen::Index r = 0; r < s.alpha.size(); ++r)
    for (const cplx e : eps) {
      const cplx a = sx * s.alpha(r), g = sz * s.gamma(r);
      rho = std::max(rho, std::abs(std::sqrt(a * a + g * g - e)));
    }
  return rho;
}

double CartesianBasis::rhoEstimate(const detail::StackLayer& layer) const {
  return ::fg::rhoEstimate(s_, cartesian(layer).epsList);
}

/// Spectral radius for the slice criterion (Plan eq:slicecrit). Exact eigenvalues up to
/// n = 800 (all 1D cases); above (2D) a deterministic power iteration, which is cheap and errs on
/// the safe side (transient growth); the a-posteriori growth check bisects slices if needed.
double spectralRadius(const MatC& M) {
  FG_PROF_SCOPE("rk.spectral_radius");
  if (M.rows() <= 800) {
    Eigen::ComplexEigenSolver<MatC> es(M, false);
    return es.eigenvalues().cwiseAbs().maxCoeff();
  }
  VecC v(M.rows());
  for (Eigen::Index k = 0; k < v.size(); ++k) v(k) = std::exp(I * (0.7 * static_cast<double>(k)));
  v.normalize();
  double rho = 0;
  for (int it = 0; it < 60; ++it) {
    VecC w = M * v;
    const double nw = w.norm();
    if (nw == 0) break;
    if (it >= 40) rho = std::max(rho, nw);
    v = w / nw;
  }
  return rho;
}

/// One element of the assembled stack (interface, modal slice or RK slice), normalized y.
struct Element {
  enum class Kind { Interface, Modes, RK } kind = Kind::Interface;
  double y0 = 0, y1 = 0;
  SMatrix S;
  resident::SMatrixD dS;  ///< resident chain (M15, step 5): S on the device instead of S
  // modal slices (RCWA, homogeneous layers)
  LayerModes modes;
  MatC G;
  // RK slices: transfer matrices from y0 to the recorded sample heights
  std::vector<double> ys;
  std::vector<MatC> T;
  // for E_y / H_y and the FFF-consistent reconstruction
  const LayerModel* model = nullptr;
  double layerY0 = 0;       ///< normalized y of the layer bottom
  double qHeight_um = -1;   ///< RCWA staircase: height at which Q/eps were taken (layer-relative)
  VecC su, sv;              ///< scattered-field formulation: source amplitudes (u at the top, v at the bottom)
  // scattered-field formulation, near fields: particular solution F_p(y) = sum p e^{kappa (y - yref)}
  // (modal slices) and the affine parts F(ys) = T F(y0) + f (RK slices)
  std::vector<VecC> pVec, f;
  std::vector<cplx> pKap;
  std::vector<double> pRef;
  VecC particular(double y) const {
    VecC F = VecC::Zero(pVec.empty() ? 0 : pVec.front().size());
    for (std::size_t j = 0; j < pVec.size(); ++j) F += pVec[j] * std::exp(pKap[j] * (y - pRef[j]));
    return F;
  }
};

struct Assembly {
  SMatrix total;
  SourcedSMatrix totalSrc;     // scattered-field formulation
  std::vector<Element> elems;  // only filled when recording
  double yTop = 0;             // normalized
};

class StackSolver final : public Solver {
 public:
  explicit StackSolver(SolverOptions o) : opt_(std::move(o)) {}

  // ------------------------------------------------------------------ efficiencies
  Result solve(const Stack& stack, const Wavelength& wl) override {
    if (stack.cylindrical.enabled) return cylindricalSolve(stack, wl);
    const Setup s = setup(stack, wl);
    return solveModels(stack, s, buildModels(stack, s));
  }

  /// solve() with prebuilt layer models (wavelength independent; Simulation, M12).
  Result solveWithModels(const Stack& stack, const Wavelength& wl, const detail::LayerModels& models) {
    FG_REQUIRE(!stack.cylindrical.enabled, "cylindrical: prebuilt layer models are not available (M20)");
    return solveModels(stack, setup(stack, wl), models);
  }

  Result solveModels(const Stack& stack, const Setup& s, const detail::LayerModels& models) {
    FG_PROF_SCOPE("solve");  // without building the layer models ("solve.layer_models")
    linalg::RefinementCounter counter;
    linalg::setRefinementCounter(&counter);
    Result res = solveCounted(stack, s, models, counter);
    linalg::setRefinementCounter(nullptr);
    // safeguard (Plan section "M16: Schutzmechanismen"): lossless structure, energy off by > 1e-3
    const bool fallback = opt_.precisionFallback < 0 ? linalg::precision() == linalg::Precision::FP32
                                                     : opt_.precisionFallback > 0;
    if (linalg::precision() != linalg::Precision::FP64 && fallback && lossless(stack) &&
        std::abs(res.diag.energyBalance - 1.0) > 1e-3) {
      SolverOptions o = opt_;
      o.precision = "fp64";
      StackSolver fp64(o);
      detail::applyThreadSettings(o);
      const std::string used = res.diag.precision;
      const double balance = res.diag.energyBalance;
      res = fp64.solveModels(stack, s, models);
      res.diag.precisionFallback = true;
      res.diag.warnings.push_back("precision " + used + ": energy balance " + std::to_string(balance) +
                                  " (lossless structure), repeated in fp64");
      detail::applyThreadSettings(opt_);
    }
    return res;
  }

  /// All materials lossless (isotropic: Im eps = 0; tensors: Hermitian).
  static bool lossless(const Stack& stack) {
    for (const auto& [name, m] : stack.materials) {
      const bool ok = m.isotropic ? std::abs(m.eps.imag()) < 1e-14
                                  : (m.epsTensor - m.epsTensor.adjoint()).cwiseAbs().maxCoeff() < 1e-14;
      if (!ok) return false;
    }
    return true;
  }

  Result solveCounted(const Stack& stack, const Setup& s, const detail::LayerModels& models,
                      const linalg::RefinementCounter& counter) {
    Result res;
    res.diag.precision = linalg::precisionName(linalg::precision());
    if (linalg::precision() == linalg::Precision::FP32 && opt_.rkTol < 1e-5 && opt_.kind == SolverKind::RK)
      res.diag.warnings.push_back("precision fp32: rk.tol raised to 1e-5");
    for (const auto& m : models) res.diag.warnings.insert(res.diag.warnings.end(), m->warnings.begin(), m->warnings.end());
    if (stack.aperiodic.modeExcitation) {
      modalOutputs(stack, s, models, res);
      res.diag.refinementSteps = counter.steps;
      res.diag.refinementFallbacks = counter.fallbacks;
      return res;
    }
    if (s.scattered && s.aperiodic) {
      scatteredOutputs(stack, s, models, res);
      res.diag.refinementSteps = counter.steps;
      res.diag.refinementFallbacks = counter.fallbacks;
      return res;
    }
    if (s.scattered) {
      const Assembly as = assemble(s, CartesianBasis(s), layerPtrs(models), false, {}, res.diag);
      {
        // total field = background + scattered (no incoming scattered amplitudes)
        const int nf = s.map.size(), c = s.map.center();
        VecC b = as.totalSrc.u, a = as.totalSrc.v;
        b(c) += s.bgTop(0);
        b(nf + c) += s.bgTop(1);
        a(c) += s.bgSub(0);
        a(nf + c) += s.bgSub(1);
        efficienciesFrom(s, b, a, res);
      }
      res.diag.refinementSteps = counter.steps;
      res.diag.refinementFallbacks = counter.fallbacks;
      return res;
    }
    if (s.aperiodic) {
      scattering(stack, s, models, res);
      res.diag.refinementSteps = counter.steps;
      res.diag.refinementFallbacks = counter.fallbacks;
      return res;
    }
    if (decoupled(stack, s)) {
      // 1D in-plane isotropic: TE and TM separately with 2 n_f (Plan section "Spezialfaelle"),
      // embedded block-diagonally into the S-matrix of the (a_x, a_z) amplitudes
      const int nf = s.map.size();
      SMatrix S;
      for (MatC* B : {&S.S11, &S.S12, &S.S21, &S.S22}) *B = MatC::Zero(2 * nf, 2 * nf);
      const VecC a_in = incidentAmplitudes(s);
      for (const Polarization p : {Polarization::TE, Polarization::TM}) {
        const int k = polarizationBlocks(p).amp[0];  // 0: x (TM), 1: z (TE)
        // not excited (|A| <= 1e-14 |A_inc|, e.g. cos 90 deg = 6e-17): its block stays zero
        if (a_in.segment(k * nf, nf).norm() <= 1e-14 * a_in.norm()) continue;
        const SMatrix Sp = assemble(s, CartesianBasis(s, p), layerPtrs(models), false, {}, res.diag).total;
        S.S11.block(k * nf, k * nf, nf, nf) = Sp.S11;
        S.S12.block(k * nf, k * nf, nf, nf) = Sp.S12;
        S.S21.block(k * nf, k * nf, nf, nf) = Sp.S21;
        S.S22.block(k * nf, k * nf, nf, nf) = Sp.S22;
      }
      efficiencies(s, S, res);
      res.diag.refinementSteps = counter.steps;
      res.diag.refinementFallbacks = counter.fallbacks;
      return res;
    }
    const Assembly as = assemble(s, CartesianBasis(s), layerPtrs(models), false, {}, res.diag);
    efficiencies(s, as.total, res);
    res.diag.refinementSteps = counter.steps;
    res.diag.refinementFallbacks = counter.fallbacks;
    return res;
  }

  /// Growth limit of the slice criterion: explicit value, else eps^-1/2 of the working precision.
  double growthLimit() const {
    if (opt_.growthLimit > 0) return opt_.growthLimit;
    return linalg::precision() == linalg::Precision::FP64 ? 1e8 : 4e3;
  }
  /// rk.tol, raised to >= 1e-5 in FP32 (Plan section "M16: Angepasste Parameter").
  double rkTol() const {
    return linalg::precision() == linalg::Precision::FP32 ? std::max(opt_.rkTol, 1e-5) : opt_.rkTol;
  }

  // ------------------------------------------------------------------ field coefficients
  FieldSampleSet fieldCoefficients(const Stack& stack, const Wavelength& wl, const std::vector<double>& y_um) override {
    if (stack.cylindrical.enabled) throw NotImplemented("cylindrical near fields (M20)");
    const Setup s = setup(stack, wl);
    const auto models = buildModels(stack, s);
    return samples(s, models, y_um, nullptr);
  }

  // ------------------------------------------------------------------ near field cuts
  NearFieldResult nearField(const Stack& stack, const Wavelength& wl, const NearFieldRequest& req) override {
    FG_PROF_SCOPE("nearfield");
    if (stack.cylindrical.enabled) throw NotImplemented("cylindrical near fields (M20)");
    const Setup s = setup(stack, wl);
    const auto models = buildModels(stack, s);
    std::vector<double> ys;
    for (int k = 0; k < req.ny; ++k)
      ys.push_back(req.ny == 1 ? req.y_min_um : req.y_min_um + (req.y_max_um - req.y_min_um) * k / (req.ny - 1));
    ys.push_back(req.horizontal_y_um);
    std::vector<const Element*> where;
    const FieldSampleSet fs = samples(s, models, ys, &where);
    FG_PROF_SCOPE("nearfield.synthesize");
    std::vector<double> xs, zs;
    for (int k = 0; k < req.nx; ++k) xs.push_back(stack.dx_um * req.periods * k / req.nx);
    if (stack.is2D())
      for (int k = 0; k < req.nz; ++k) zs.push_back(stack.dz_um * req.periods * k / req.nz);
    NearFieldResult out;
    // vertical cut: map (x, y) at z = z0; horizontal cut: line u(x) (1D) or map (x, z) (2D) at y = y0
    out.vertical = synthesize(stack, s, fs, where, xs, {req.vertical_z_um}, 0, req.ny, req, "vertical");
    out.horizontal = synthesize(stack, s, fs, where, xs, stack.is2D() ? zs : std::vector<double>{req.vertical_z_um}, req.ny,
                                req.ny + 1, req, "horizontal");
    if (stack.is2D()) {
      out.horizontal.axis2Name = "z";
      out.horizontal.axis2_um = zs;
    } else {
      out.horizontal.axis2_um.clear();
    }
    return out;
  }

 private:
  // ------------------------------------------------------------------ cylindrical (M20)
  /// Beam on the axis or plane wave through a stack with ring layers (Plan sections "M20: Anregung",
  /// "M20: S-Matrix, Loeser und Nahfeld"): one assembly per azimuthal harmonic with the Hankel basis,
  /// powers from Plan eq:cyl_flux.
  Result cylindricalSolve(const Stack& stack, const Wavelength& wl) {
    FG_PROF_SCOPE("solve.cylindrical");
    detail::applyThreadSettings(opt_);
    setEigenThreads({1, sliceThreads(), std::max(1, opt_.threadsBlas)});
    const CylindricalSpec& c = stack.cylindrical;
    FG_REQUIRE(!stack.aperiodic.enabled && !stack.superstrateStructure && !stack.substrateStructure &&
                   !opt_.scatteredField,
               "cylindrical: aperiodic, structured half spaces and scattered_field are not available (M20)");
    FG_REQUIRE(c.planeWave || c.waist_um > 0, "cylindrical: incidence gaussian_beam with waist_um > 0");
    Result res;
    res.diag.precision = linalg::precisionName(linalg::precision());
    Setup s;
    s.k0 = 2.0 * PI / wl.lambda_um;
    const Material& sup = stack.material(stack.superstrate);
    const Material& sub = stack.material(stack.substrate);
    FG_REQUIRE(sup.isotropic && std::abs(sup.eps.imag()) < 1e-14 && sup.eps.real() > 0,
               "superstrate must be isotropic and lossless");
    FG_REQUIRE(sub.isotropic, "cylindrical: the substrate must be isotropic (M20)");
    const auto layers = detail::buildCylLayers(stack, s.k0);
    std::vector<const detail::StackLayer*> ptrs;
    double yTop = 0, Rs_um = 2.0 * c.waist_um, nLayers = 0;
    std::vector<cplx> metals, dielectrics{sup.eps, sub.eps};
    for (std::size_t k = 0; k < layers.size(); ++k) {
      ptrs.push_back(layers[k].get());
      yTop += layers[k]->thickness_um() * s.k0;
      if (const auto* r = std::get_if<RingsLayer>(&stack.layers[k].spec)) Rs_um = std::max(Rs_um, r->radii_um.back());
      for (const cplx e : layers[k]->epsList) {
        nLayers = std::max(nLayers, std::sqrt(e).real());
        (e.real() < 0 ? metals : dielectrics).push_back(e);
      }
    }
    // spectrum (Plan sec:cylgrid): branch points of the half spaces, complex path over [k_a, k_b] when
    // the background can carry guided modes (a layer of higher index than both half spaces, or a metal)
    const double nSup = std::sqrt(sup.eps).real(), nSub = std::sqrt(sub.eps).real(), ka = std::max(nSup, nSub);
    const double nMax = std::max(ka, nLayers);
    cyl::SpectrumOptions o;
    o.kmax = c.kappaMax * nMax;
    o.Rs = Rs_um * s.k0;
    o.dkFactor = c.dkFactor;
    o.nodesPerPanel = c.nodesPerPanel;
    o.branchPoints = {nSub, nSup};
    double kb = 1.2 * nMax;
    for (const cplx em : metals)
      for (const cplx ed : dielectrics) kb = std::max(kb, 1.2 * std::sqrt(em * ed / (em + ed)).real());
    const bool guided = nLayers > ka || !metals.empty();
    if (c.deform == CylindricalSpec::Deform::On || (c.deform == CylindricalSpec::Deform::Auto && guided)) {
      FG_REQUIRE(kb < o.kmax, "cylindrical: k_b of the deformed path must be below k_max (raise kappa_max)");
      o.deformA = ka;
      o.deformB = kb;
      o.delta = c.deltaRel * (kb - ka);
    }
    const cyl::SpectralGrid grid = cyl::spectralGrid(o);
    const Eigen::Index n = grid.k.size();
    const cyl::CylFactorization fac = c.liR ? cyl::CylFactorization::LiR : cyl::CylFactorization::Laurent;
    if (c.planeWave) {
      FG_REQUIRE(Rs_um > 0, "cylindrical plane wave: the stack needs a rings layer");
      cylindricalPlaneWave(stack, wl, s, layers, ptrs, grid, fac, res);
      return res;
    }

    // incident beam (Plan sec:cylexc): spectra of order 0 (w0^2/2 e^{-k^2 w0^2/4}) and order 1
    // (k w0^4/4 e^{-k^2 w0^2/4}, of r e^{-r^2/w0^2}); phase referenced at the top, aperture <= 1
    const double w0 = c.waist_um * s.k0, yf = c.focus_y_um * s.k0;
    VecC g0 = VecC::Zero(n), g1 = VecC::Zero(n);
    for (Eigen::Index j = 0; j < n; ++j) {
      const cplx k = grid.k(j);
      if (k.imag() != 0.0 || k.real() >= nSup) continue;
      const cplx ph = std::exp(I * std::sqrt(sup.eps - k * k) * (yf - yTop)) * std::exp(-k * k * w0 * w0 / 4.0);
      g0(j) = 0.5 * w0 * w0 * ph;
      g1(j) = 0.25 * k * w0 * w0 * w0 * w0 * ph;
    }
    // harmonic -> (a_+, a_-) by Plan eq:cyl_pm
    auto incident = [&](int m, VecC& ap, VecC& am) {
      ap = VecC::Zero(n);
      am = VecC::Zero(n);
      const std::string& p = c.polarization;
      if (p == "x") {
        if (m == -1) ap = g0;
        if (m == 1) am = g0;
      } else if (p == "z") {
        if (m == -1) ap = I * g0;
        if (m == 1) am = -I * g0;
      } else if (p == "circular_plus") {  // E_z = +i E_x: only E_x - i E_z = 2 f
        if (m == 1) am = 2.0 * g0;
      } else if (p == "circular_minus") {
        if (m == -1) ap = 2.0 * g0;
      } else if (p == "radial") {  // E_r = f: E_+ = E_- = f; spectrum of order -1 = -(order 1)
        if (m == 0) ap = g1, am = -g1;
      } else if (p == "azimuthal") {  // E_phi = f: E_pm = +-i f
        if (m == 0) ap = I * g1, am = I * g1;
      } else {
        throw Error("cylindrical: polarization must be x, z, circular_plus, circular_minus, radial or azimuthal");
      }
    };
    std::vector<int> ms = c.harmonics;
    if (ms.empty()) {
      const std::string& p = c.polarization;
      ms = (p == "radial" || p == "azimuthal") ? std::vector<int>{0}
           : p == "circular_plus"              ? std::vector<int>{1}
           : p == "circular_minus"             ? std::vector<int>{-1}
                                               : std::vector<int>{-1, 1};
    }
    double pInc = 0, pR = 0, pT = 0;
    std::vector<double> rm, tm;
    for (const int m : ms) {
      VecC ap, am;
      incident(m, ap, am);
      VecC aIn(2 * n);
      aIn << ap, am;
      if (aIn.isZero(0.0)) {
        rm.push_back(0);
        tm.push_back(0);
        continue;
      }
      const detail::CylBasis basis(grid, m, sup.eps, sub.eps, layers, fac);
      const Assembly as = assemble(s, basis, ptrs, false, {}, res.diag);
      const VecC b = as.total.S12 * aIn, aSub = as.total.S22 * aIn;
      VecC V = VecC::Zero(4 * n);
      V.head(2 * n) = aIn;
      pInc -= cyl::flux(grid, basis.ref().Psi * V);
      V.setZero();
      V.tail(2 * n) = b;
      rm.push_back(cyl::flux(grid, basis.ref().Psi * V));
      V.setZero();
      V.head(2 * n) = aSub;
      tm.push_back(-cyl::flux(grid, basis.sub().Psi * V));
      pR += rm.back();
      pT += tm.back();
    }
    FG_REQUIRE(pInc > 0, "cylindrical: the incident beam carries no power (harmonics?)");
    CylindricalResult& cr = res.cylindrical;
    cr.present = true;
    cr.harmonics = ms;
    for (double& v : rm) v /= pInc;
    for (double& v : tm) v /= pInc;
    cr.reflectedPerHarmonic = rm;
    cr.transmittedPerHarmonic = tm;
    cr.reflected = pR / pInc;
    cr.transmitted = pT / pInc;
    cr.absorbed = 1.0 - cr.reflected - cr.transmitted;
    cr.nodes = static_cast<int>(n);
    res.absorption = cr.absorbed;
    res.diag.energyBalance = cr.reflected + cr.transmitted;
    return res;
  }

  /// Plane wave on a cylindrical stack (Plan sec:cylexc, "Ebene Welle (Streufeld)"): the background
  /// stack (rings -> outer material) at the incidence order (Airy), the sources of the ring layers from
  /// Jacobi-Anger (Plan eq:cyl_jacobi, eq:cyl_scat), scattered powers up and down per harmonic.
  void cylindricalPlaneWave(const Stack& stack, const Wavelength& wl, Setup& s,
                            const std::vector<std::unique_ptr<detail::CylLayer>>& layers,
                            const std::vector<const detail::StackLayer*>& ptrs, const cyl::SpectralGrid& grid,
                            cyl::CylFactorization fac, Result& res) {
    const CylindricalSpec& c = stack.cylindrical;
    const Material& sup = stack.material(stack.superstrate);
    const Material& sub = stack.material(stack.substrate);
    const Eigen::Index n = grid.k.size();
    Stack cs = stack;  // Cartesian background problem (one order, the period is not used)
    cs.cylindrical = {};
    cs.dx_um = 1.0;
    cs.dz_um = 0.0;
    s.map = IndexMap(0, 0);
    s.epsSup = sup.eps;
    s.inc = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
    s.alpha = VecC::Constant(1, s.inc.alpha0);
    s.gamma = VecC::Constant(1, s.inc.gamma0);
    s.sub = &sub;
    s.scattered = true;
    background(cs, wl, s);
    const double kpar = std::hypot(s.inc.alpha0, s.inc.gamma0);
    const double Psi = kpar > 0 ? std::atan2(s.inc.gamma0, s.inc.alpha0) : 0.0;
    double Rs = 0;
    for (const auto& L : layers)
      if (!L->homogeneous()) Rs = std::max(Rs, L->rings.radii.back());
    std::vector<int> ms = c.harmonics;
    if (ms.empty()) {
      // the source of harmonic m vanishes for |m| > m_cut (Plan eq:cyl_scat); normal incidence: m = +-1
      const int mc = kpar < 1e-12 ? 1 : cyl::harmonicCutoff(kpar * Rs);
      for (int m = -mc; m <= mc; ++m)
        if (kpar >= 1e-12 || m != 0) ms.push_back(m);
    }
    CylindricalResult& cr = res.cylindrical;
    cr.present = true;
    cr.planeWave = true;
    cr.harmonics = ms;
    cr.nodes = static_cast<int>(n);
    const double toArea = 1.0 / (s.incIntensity * s.k0 * s.k0);  // normalized power -> um^2
    for (const int m : ms) {
      detail::CylBasis basis(grid, m, sup.eps, sub.eps, layers, fac);
      basis.setPlaneWave(kpar, Psi);
      s.cylSource = [&](const detail::StackLayer& l, const Setup::BgLayer& L, const VecC& w) {
        const cplx ey = computeEy(L.Q1, s.alpha, s.gamma, w)(0);
        return basis.source(l, w(0) + I * w(1), w(0) - I * w(1), ey);
      };
      const Assembly as = assemble(s, basis, ptrs, false, {}, res.diag);
      VecC V = VecC::Zero(4 * n);
      V.tail(2 * n) = as.totalSrc.u;  // scattered field, nothing comes in
      const double up = cyl::flux(grid, basis.ref().Psi * V);
      V.setZero();
      V.head(2 * n) = as.totalSrc.v;
      const double down = -cyl::flux(grid, basis.sub().Psi * V);
      cr.reflectedPerHarmonic.push_back(up * toArea);
      cr.transmittedPerHarmonic.push_back(down * toArea);
      cr.scatteredUp_um2 += up * toArea;
      cr.scatteredDown_um2 += down * toArea;
    }
    s.cylSource = nullptr;
    cr.scattering_um2 = cr.scatteredUp_um2 + cr.scatteredDown_um2;
  }

  detail::LayerModels buildModels(const Stack& stack, const Setup& s) const {
    FG_PROF_SCOPE("solve.layer_models");
    return detail::buildLayerModels(stack, s.map, opt_.rcwaSlices);
  }

  // ------------------------------------------------------------------ setup
  Setup setup(const Stack& stack, const Wavelength& wl) const {
    detail::applyThreadSettings(opt_);
    setEigenThreads({1, sliceThreads(), std::max(1, opt_.threadsBlas)});
    Setup s;
    s.map = IndexMap(opt_.Nx, stack.is2D() ? opt_.Nz : 0);
    s.k0 = 2.0 * PI / wl.lambda_um;
    const Material& sup = stack.material(stack.superstrate);
    FG_REQUIRE(sup.isotropic && std::abs(sup.eps.imag()) < 1e-14 && sup.eps.real() > 0,
               "superstrate must be isotropic and lossless");
    s.epsSup = sup.eps;
    s.inc = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
    const double Kx = 2.0 * PI / (s.k0 * stack.dx_um);
    const double Kz = stack.is2D() ? 2.0 * PI / (s.k0 * stack.dz_um) : 0.0;
    const int nf = s.map.size();
    // oblique incidence on a PML axis: scattered-field formulation without Bloch phase on that axis
    const AperiodicSpec& ap = stack.aperiodic;
    const bool oblX = ap.enabled && ap.x.active() && std::abs(s.inc.alpha0) >= 1e-12;
    const bool oblZ = ap.enabled && ap.z.active() && std::abs(s.inc.gamma0) >= 1e-12;
    s.scattered = opt_.scatteredField || oblX || oblZ;
    s.blochSource = oblX || oblZ;
    s.alpha.resize(nf);
    s.gamma.resize(nf);
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      s.alpha(r) = (oblX ? 0.0 : s.inc.alpha0) + n * Kx;
      s.gamma(r) = (oblZ ? 0.0 : s.inc.gamma0) + m * Kz;
    }
    if (s.blochSource) {
      s.alphaB.resize(nf);
      s.gammaB.resize(nf);
      for (int r = 0; r < nf; ++r) {
        const auto [n, m] = s.map.orders(r);
        s.alphaB(r) = s.inc.alpha0 + n * Kx;
        s.gammaB(r) = s.inc.gamma0 + m * Kz;
      }
      // T(m, n) = 1/d int_phys exp(i k0 (alpha0 + (n - m) K) x) dx on the oblique PML axes
      auto restrict = [&](int N, double d, double w, double a0, double K) {
        MatC T(2 * N + 1, 2 * N + 1);
        const double xa = w, xb = d - w;
        for (int m = -N; m <= N; ++m)
          for (int n = -N; n <= N; ++n) {
            const double kap = s.k0 * (a0 + (n - m) * K);
            T(m + N, n + N) = std::abs(kap) * (xb - xa) < 1e-12
                                  ? cplx((xb - xa) / d)
                                  : (std::exp(I * (kap * xb)) - std::exp(I * (kap * xa))) / (I * kap * d);
          }
        return T;
      };
      MatC Tx, Tz;
      if (oblX) Tx = restrict(s.map.Nx(), stack.dx_um, ap.x.width_um, s.inc.alpha0, Kx);
      if (oblZ) Tz = restrict(s.map.Nz(), stack.dz_um, ap.z.width_um, s.inc.gamma0, Kz);
      s.Tsrc = separable(s.map, oblX ? &Tx : nullptr, oblZ ? &Tz : nullptr);
    }
    s.ref = isotropicModes(s.epsSup, s.alpha, s.gamma);
    s.sub = &stack.material(stack.substrate);
    s.subModes = s.sub->isotropic ? isotropicModes(s.sub->eps, s.alpha, s.gamma)
                                  : anisotropicModes(s.sub->epsTensor, s.alpha, s.gamma);
    s.tr = resolveTransform(stack);
    if (ap.enabled) {
      FG_REQUIRE(!s.tr.x && !s.tr.z, "aperiodic: PML and coordinate_transform (ASR) cannot be combined");
      FG_REQUIRE(!ap.z.active() || stack.is2D(), "aperiodic: pml.z needs a 2D grating");
      checkPmlZones(stack);
      // oblique incidence: scattered-field formulation (set above); the total-field formulation (Plan
      // "Einfall (a)") would damp the incident wave in the zones
      s.aperiodic = true;
      if (ap.x.active()) s.pmlX = std::make_shared<const PmlMap>(stack.dx_um, ap.x.width_um, ap.x.beta);
      if (ap.z.active()) s.pmlZ = std::make_shared<const PmlMap>(stack.dz_um, ap.z.width_um, ap.z.beta);
    }
    if (s.tr.active() || s.aperiodic) transform(s);
    if (s.scattered) {
      FG_REQUIRE(!s.tr.x && !s.tr.z, "scattered-field formulation: coordinate_transform (ASR) is not supported");
      background(stack, wl, s);
    }
    if (stack.superstrateStructure || stack.substrateStructure) {
      FG_REQUIRE(s.aperiodic, "structured superstrate/substrate need an aperiodic cell (PML)");
      FG_REQUIRE(!s.scattered, "structured superstrate/substrate: the scattered-field formulation is not available");
      auto halfSpace = [&](const LayerEntry& le, std::optional<QBlocks>& Qout) {
        Stack t = stack;
        t.layers = {le};
        const auto m = detail::buildLayerModel(t, 0, s.map, opt_.rcwaSlices);
        FG_REQUIRE(m->yInvariant(), "a structured superstrate/substrate must be y-invariant");
        Qout = m->Q(0.0);
        return numericBasis(propagation(s, std::nullopt, *Qout));
      };
      if (stack.superstrateStructure) s.ref = halfSpace(*stack.superstrateStructure, s.supQ);
      if (stack.substrateStructure) s.subModes = halfSpace(*stack.substrateStructure, s.subQ);
    }
    FG_REQUIRE(!ap.modeExcitation || stack.superstrateStructure,
               "aperiodic.excitation \"mode\" needs a structured superstrate (the incident waveguide)");
    FG_REQUIRE(!ap.modeExcitation || !s.scattered, "modal excitation: phi must be 0 (the angles are ignored)");
    return s;
  }

  /// Background field (Plan section "Streufeldformulierung"): the background stack (no PML, no
  /// transformation) at the incidence order only, as exponential modes per layer, and the propagation
  /// matrices of the background materials for sigma = (M - M_bg) F_bg.
  void background(const Stack& stack, const Wavelength& wl, Setup& s) const {
    FG_PROF_SCOPE("solve.background");
    Stack bg = backgroundStack(stack);
    bg.aperiodic = {};
    bg.transform = {};
    SolverOptions o = opt_;
    o.Nx = 0;
    o.Nz = 0;
    o.kind = SolverKind::RCWA;
    o.scatteredField = false;
    o.splitPolarizations = false;
    StackSolver sub(o);
    std::vector<double> ys{0.0};
    for (const auto& le : bg.layers) ys.push_back(ys.back() + std::get<HomogeneousLayer>(le.spec).thickness_um);
    const FieldSampleSet fs = sub.fieldCoefficients(bg, wl, ys);
    detail::applyThreadSettings(opt_);
    const VecC a0 = VecC::Constant(1, s.inc.alpha0), g0 = VecC::Constant(1, s.inc.gamma0);
    auto modes1 = [&](const Material& m) {
      return m.isotropic ? isotropicModes(m.eps, a0, g0) : anisotropicModes(m.epsTensor, a0, g0);
    };
    const int nf = s.map.size();
    s.bg.clear();
    for (std::size_t k = 0; k < bg.layers.size(); ++k) {
      const Material& mat = bg.material(std::get<HomogeneousLayer>(bg.layers[k].spec).material);
      Setup::BgLayer L;
      L.y0 = ys[k] * s.k0;
      L.y1 = ys[k + 1] * s.k0;
      const ModeBasis mb = modes1(mat);
      L.Psi = mb.Psi;
      L.lambda = mb.lambda;
      const VecC c0 = mb.PsiInv * fs.samples[k].F, c1 = mb.PsiInv * fs.samples[k + 1].F;
      L.amp.resize(4);
      L.amp << c1.head(2), c0.tail(2);
      const QBlocks Q = mat.isotropic ? homogeneousIsotropicQ(mat.eps, nf) : homogeneousTensorQ(mat.epsTensor, nf);
      L.M = propagation(s, std::nullopt, Q);
      if (s.blochSource) L.MB = buildPropagationMatrix(Q, s.alphaB, s.gammaB);
      L.Q = Q;
      L.Q1 = mat.isotropic ? homogeneousIsotropicQ(mat.eps, 1) : homogeneousTensorQ(mat.epsTensor, 1);
      s.bg.push_back(std::move(L));
    }
    const ModeBasis sup = isotropicModes(s.epsSup, a0, g0), sb = modes1(*s.sub);
    s.bgTop = (sup.PsiInv * fs.samples.back().F).tail(2);
    s.bgSub = (sb.PsiInv * fs.samples.front().F).head(2);
    VecC ainc(2);
    ainc << s.inc.A(0), s.inc.A(2);
    const VecC Fi = sup.Psi.leftCols(2) * ainc;
    const double beta0 = std::sqrt(std::max(s.epsSup.real() - s.inc.alpha0 * s.inc.alpha0 - s.inc.gamma0 * s.inc.gamma0, 1e-300));
    s.incIntensity = std::abs(0.5 * std::real(Fi(1) * std::conj(Fi(2)) - Fi(0) * std::conj(Fi(3)))) *
                     std::sqrt(s.epsSup.real()) / beta0;
  }

  // ------------------------------------------------------------------ scattered-field sources
  /// Background layer containing the normalized height y.
  static const Setup::BgLayer& bgLayer(const Setup& s, double y) {
    for (const auto& L : s.bg)
      if (y <= L.y1) return L;
    return s.bg.back();
  }
  /// Background field (4 components, incidence order) at normalized y.
  static VecC bgField(const Setup::BgLayer& L, double y) {
    VecC f = VecC::Zero(4);
    for (int k = 0; k < 4; ++k) f += L.Psi.col(k) * (L.amp(k) * std::exp(L.lambda(k) * (y - (k < 2 ? L.y1 : L.y0))));
    return f;
  }
  /// sigma = (M - M_bg) F_bg for a single-order background field w (Plan eq:scat_ode in the form
  /// F_s' = M F_s + (M - M_bg) F_bg); oblique PML: in the Bloch basis (MB), then re-expanded.
  static VecC sourceOf(const Setup& s, const detail::StackLayer& layer, const MatC& M, const MatC* MB,
                       const Setup::BgLayer& L, const VecC& w) {
    if (s.cylSource) return s.cylSource(layer, L, w);
    const int nf = s.map.size(), c = s.map.center();
    VecC sig = VecC::Zero(4 * nf);
    for (int k = 0; k < 4; ++k)
      if (w(k) != 0.0)
        sig += w(k) * (s.blochSource ? VecC(MB->col(k * nf + c) - L.MB.col(k * nf + c)) : VecC(M.col(k * nf + c) - L.M.col(k * nf + c)));
    return s.blochSource ? componentwise(s.Tsrc, sig) : sig;
  }
  /// Background layers boundaries strictly inside (ya, yb) (normalized): elements must not cross them.
  static std::vector<double> bgCuts(const Setup& s, double ya, double yb) {
    std::vector<double> c{ya};
    if (s.scattered)
      for (const auto& L : s.bg)
        if (L.y1 > ya + 1e-12 * std::max(1.0, yb) && L.y1 < yb - 1e-12 * std::max(1.0, yb)) c.push_back(L.y1);
    c.push_back(yb);
    return c;
  }
  /// Source amplitudes of a modal slice [y0, y1] (Plan "RCWA: partikulaere Loesung je Exponentialterm"):
  /// F_p = sum (kappa - M)^-1 sigma_kappa e^{kappa y}; with g = Psi_ref^-1 F_p at both ends,
  /// u = g1_b - S11 g0_b - S12 g1_a, v = g0_a - S21 g0_b - S22 g1_a. False if kappa is within 1e-8 max|lambda|
  /// of an eigenvalue of M (then the slice is integrated with RK).
  static bool modalSources(const Setup& s, const detail::StackLayer& layer, const MatC& M, const MatC* MB, const LayerModes& modes, const SMatrix& S,
                           double y0, double y1, const ModeBasis& ref, VecC& u, VecC& v, Element* rec = nullptr) {
    const Setup::BgLayer& L = bgLayer(s, 0.5 * (y0 + y1));
    const double lmax = std::max(modes.lu.cwiseAbs().maxCoeff(), modes.ld.cwiseAbs().maxCoeff());
    const Eigen::Index N = M.rows();
    VecC F0 = VecC::Zero(N), F1 = VecC::Zero(N);
    for (int k = 0; k < 4; ++k) {
      if (L.amp(k) == 0.0) continue;
      const cplx kap = L.lambda(k);
      double dist = 1e300;
      for (Eigen::Index j = 0; j < modes.lu.size(); ++j)
        dist = std::min({dist, std::abs(kap - modes.lu(j)), std::abs(kap - modes.ld(j))});
      if (dist < 1e-8 * lmax) return false;
      const VecC p = (kap * MatC::Identity(N, N) - M).partialPivLu().solve(sourceOf(s, layer, M, MB, L, L.Psi.col(k))) * L.amp(k);
      const double yr = k < 2 ? L.y1 : L.y0;
      F0 += p * std::exp(kap * (y0 - yr));
      F1 += p * std::exp(kap * (y1 - yr));
      if (rec) {
        rec->pVec.push_back(p);
        rec->pKap.push_back(kap);
        rec->pRef.push_back(yr);
      }
    }
    const Eigen::Index n = S.nh();
    const VecC g0 = ref.PsiInv * F0, g1 = ref.PsiInv * F1;
    u = g1.tail(n) - S.S11 * g0.tail(n) - S.S12 * g1.head(n);
    v = g0.head(n) - S.S21 * g0.tail(n) - S.S22 * g1.head(n);
    return true;
  }
  /// Source amplitudes from an affine transfer F(y1) = T F(y0) + f (RK): u = g_b - S12 g_a, v = -S22 g_a.
  static void affineSources(const SMatrix& S, const VecC& f, const ModeBasis& ref, VecC& u, VecC& v) {
    const Eigen::Index n = S.nh();
    const VecC g = ref.PsiInv * f;
    u = g.tail(n) - S.S12 * g.head(n);
    v = -S.S22 * g.head(n);
  }
  /// RK transfer of dF/dy = M F + sigma(y) over [ya, yb] as the augmented linear system (F; 1).
  void affineTransfer(const std::function<MatC(double)>& Mfun, const std::function<VecC(double, const MatC&)>& sigma,
                      double ya, double yb, const ButcherTableau& tab, double rho, MatC& T, VecC& f, Diagnostics& diag) const {
    auto aug = [&](double y) {
      const MatC Mm = Mfun(y);
      const Eigen::Index n = Mm.rows();
      MatC A = MatC::Zero(n + 1, n + 1);
      A.topLeftCorner(n, n) = Mm;
      A.topRightCorner(n, 1) = sigma(y, Mm);
      return A;
    };
    RKStats st;
    const MatC Ta = rkTransfer(aug, ya, yb, tab, rkTol(), rho, st);
    diag.rkSteps += st.steps;
    const Eigen::Index n = Ta.rows() - 1;
    T = Ta.topLeftCorner(n, n);
    f = Ta.topRightCorner(n, 1);
  }

  /// Transformed coordinates (Plan section "M17", eq:ct_deriv): A, G, the conversion matrices of
  /// eq:ct_rayleigh / eq:ct_inc and numerical bases of the half spaces (in X).
  static void transform(Setup& s) {
    FG_PROF_SCOPE("solve.transform");
    s.transformed = true;
    const int nf = s.map.size();
    auto t = detail::transformedCoordinates(s.map, s.tr, s.pmlX.get(), s.pmlZ.get(), s.alpha, s.gamma,
                                            s.inc.alpha0 * s.k0, s.inc.gamma0 * s.k0);
    s.A = std::move(t.A);
    s.G = std::move(t.G);
    s.toPhys = std::move(t.toPhys);
    s.fromPhys = std::move(t.fromPhys);
    s.refPhys = std::move(s.ref);
    s.subPhys = std::move(s.subModes);
    s.ref = numericBasis(buildPropagationMatrix(homogeneousIsotropicQ(s.epsSup, nf), s.A, s.G));
    s.subModes = numericBasis(buildPropagationMatrix(
        s.sub->isotropic ? homogeneousIsotropicQ(s.sub->eps, nf) : homogeneousTensorQ(s.sub->epsTensor, nf), s.A, s.G));
  }

  /// True if the problem decouples exactly into TE and TM (Plan section "Spezialfaelle"): 1D,
  /// gamma = 0 for all orders, isotropic substrate and layers (the superstrate is always isotropic).
  bool decoupled(const Stack& stack, const Setup& s) const {
    if (!opt_.splitPolarizations || stack.is2D() || !s.gamma.isZero(0.0) || !s.sub->isotropic || s.transformed ||
        s.scattered)
      return false;
    for (const auto& le : stack.layers) {
      const auto names = layerMaterials(le.spec);
      if (!names) return false;
      for (const auto& n : *names)
        if (!stack.material(n).isotropic) return false;
    }
    return true;
  }

  // ------------------------------------------------------------------ assembly
  /// Thread-local settings of this solve, re-applied on worker threads (backend, FFT sampling,
  /// precision, refinement counter).
  struct Settings {
    const SolverOptions* opt;
    linalg::RefinementCounter* counter;
    void apply() const {
      detail::applyThreadSettings(*opt);
      linalg::setRefinementCounter(counter);
    }
  };
  Settings currentSettings() const { return {&opt_, linalg::refinementCounter()}; }

  /// One unit of work of the assembly: the element(s) of a homogeneous layer, an RCWA slice or an
  /// RK stability slice (bisection may give several), in stack order.
  using Task = std::function<void(std::vector<Element>&, Diagnostics&)>;

  /// Threads for the layers/slices of one computation (Plan section "M15", step 2).
  int sliceThreads() const {
    return opt_.threadsSlices > 0 ? opt_.threadsSlices : threadPlan(opt_, 1, availableThreads()).slices;
  }

  /// Builds the total S-matrix. With record = true every element is kept (S, modes, RK transfer
  /// matrices at the requested normalized heights 'sampleY', sorted).
  /// Reduced systems (sys.polarization() set) only for the total S-matrix (record = false).
  /// The elements are independent (Plan section "M15", step 3): they are computed in blocks of
  /// sliceThreads() tasks in parallel and chained in stack order (sequential or fixed-shape tree).
  Assembly assemble(const Setup& s, const detail::ModalBasis& sys, const std::vector<const detail::StackLayer*>& layers,
                    bool record, const std::vector<double>& sampleY, Diagnostics& diag) {
    FG_PROF_SCOPE("solve.assemble");
    FG_REQUIRE(!record || !sys.polarization(), "field records need the full 4 n_f system");
    // resident chain (Plan section "M15", step 5): slice modes, slice S-matrices and star products
    // stay on the device; only the total S-matrix is downloaded (not for near-field records)
    const bool res = !record && !s.scattered && linalg::gpuSelected(sys.ref().Psi.rows()) && linalg::gpu().resident();
    std::shared_ptr<const linalg::DeviceMatrix> psiInvD;
    if (res) {
      linalg::Backend& d = linalg::gpu();
      psiInvD = std::make_shared<const linalg::DeviceMatrix>(d.upload(sys.ref().PsiInv));
      d.synchronize();  // read by the task threads (their own streams)
    }
    std::vector<Task> tasks;
    double y = 0;
    for (const detail::StackLayer* mp : layers) {
      const detail::StackLayer& m = *mp;
      const double H = m.thickness_um() * s.k0;
      const double y0 = y;
      if (m.homogeneous()) {
        tasks.push_back([&s, &sys, &m, y0, H, record, res, psiInvD](std::vector<Element>& out, Diagnostics&) {
          auto modes = [&] { return sys.homogeneousModes(m); };
          Element e;
          e.kind = Element::Kind::Modes;
          e.y0 = y0;
          e.y1 = y0 + H;
          e.model = dynamic_cast<const LayerModel*>(&m);
          e.layerY0 = y0;
          if (res) {
            linalg::Backend& d = linalg::gpu();
            e.dS = resident::layerSMatrix(d, resident::upload(d, modes()), H, *psiInvD);
            d.synchronize();
            out.push_back(std::move(e));
            return;
          }
          RCWASliceRecord rec;
          e.S = stableLayerSMatrix(modes(), H, sys.ref(), record ? &rec : nullptr);
          e.modes = std::move(rec.modes);
          e.G = std::move(rec.G);
          out.push_back(std::move(e));
        });
      } else if (opt_.kind == SolverKind::RCWA) {
        const int K = m.yInvariant() ? 1 : m.rcwaSlices;
        for (int k = 0; k < K; ++k)
          tasks.push_back([this, &s, &sys, &m, y0, H, K, k, record, res, psiInvD](std::vector<Element>& out, Diagnostics& dg) {
            const double yq = (k + 0.5) * m.thickness_um() / K;  // staircase: Q at slice centre
            const MatC M = sys.M(m, yq);
            if (s.scattered) {
              scatteredSlices(s, sys, m, M, y0, y0 + H * k / K, y0 + H * (k + 1) / K, yq, record, out, dg);
              return;
            }
            Element e;
            e.kind = Element::Kind::Modes;
            e.y0 = y0 + H * k / K;
            e.y1 = y0 + H * (k + 1) / K;
            e.model = dynamic_cast<const LayerModel*>(&m);
            e.layerY0 = y0;
            e.qHeight_um = yq;
            if (res) {
              linalg::Backend& d = linalg::gpu();
              e.dS = resident::layerSMatrix(d, resident::layerModes(d, M, sys.polarization()), H / K, *psiInvD);
              d.synchronize();
              out.push_back(std::move(e));
              return;
            }
            RCWASliceRecord rec;
            e.S = stableLayerSMatrix(sys.modes(M), H / K, sys.ref(), record ? &rec : nullptr);
            e.modes = std::move(rec.modes);
            e.G = std::move(rec.G);
            out.push_back(std::move(e));
          });
      } else {
        rkLayer(s, sys, m, y0, H, record, sampleY, tasks);
      }
      y += H;
    }

    Assembly as;
    as.yTop = y;
    // with the GPU: at most deviceStreams() threads (one stream and its workspace each)
    const int threads = res ? std::min(sliceThreads(), linalg::deviceStreams()) : sliceThreads();
    const Settings settings = currentSettings();
    BasicStarChain<SMatrix> chain(opt_.treeReduction, threads, [settings](const SMatrix& a, const SMatrix& b) {
      settings.apply();  // products of the tree run on worker threads
      return star(a, b);
    });
    BasicStarChain<resident::SMatrixD> dchain(opt_.treeReduction, threads,
                                              [settings](const resident::SMatrixD& a, const resident::SMatrixD& b) {
                                                settings.apply();
                                                linalg::Backend& d = linalg::gpu();
                                                resident::SMatrixD c = resident::star(d, a, b);
                                                d.synchronize();
                                                return c;
                                              });
    BasicStarChain<SourcedSMatrix> schain(opt_.treeReduction, threads,
                                          [settings](const SourcedSMatrix& a, const SourcedSMatrix& b) {
                                            settings.apply();
                                            return star(a, b);
                                          });
    auto push = [&](Element& e) {
      if (s.scattered) {
        const int nh = e.S.nh();
        schain.push({e.S, e.su.size() ? e.su : VecC(VecC::Zero(nh)), e.sv.size() ? e.sv : VecC(VecC::Zero(nh))});
        return;
      }
      if (!res) {
        chain.push(record ? e.S : std::move(e.S));
        return;
      }
      if (e.dS.S11.empty()) {  // host element (substrate, RK slice)
        linalg::Backend& d = linalg::gpu();
        e.dS = resident::upload(d, e.S);
        d.synchronize();
      }
      dchain.push(std::move(e.dS));
    };
    {
      Element e0;
      e0.S = fromTransfer(linalg::product(sys.ref().PsiInv, sys.sub().Psi));  // substrate -> reference (Plan eq:Stot)
      push(e0);
      if (record) as.elems.push_back(std::move(e0));
    }
    // blocks of `threads` tasks bound the memory to O(threads) S-matrices
    const int nt = static_cast<int>(tasks.size());
    for (int b0 = 0; b0 < nt; b0 += threads) {
      const int nb = std::min(threads, nt - b0);
      std::vector<std::vector<Element>> out(static_cast<std::size_t>(nb));
      std::vector<Diagnostics> dg(static_cast<std::size_t>(nb));
      parallelFor(nb, threads, [&](int i) {
        settings.apply();  // thread-local backend / FFT / precision settings
        tasks[static_cast<std::size_t>(b0 + i)](out[static_cast<std::size_t>(i)], dg[static_cast<std::size_t>(i)]);
      });
      for (int i = 0; i < nb; ++i) {
        diag.rkSteps += dg[static_cast<std::size_t>(i)].rkSteps;
        for (Element& e : out[static_cast<std::size_t>(i)]) {
          push(e);
          if (record) as.elems.push_back(std::move(e));
        }
      }
    }
    if (s.scattered) {
      as.totalSrc = schain.result();
      as.total = as.totalSrc.S;
    } else if (res) {
      linalg::Backend& d = linalg::gpu();
      as.total = resident::download(d, dchain.result());
    } else {
      as.total = chain.result();
    }
    return as;
  }

  /// RK treatment of one layer: stability sub-slices (Plan eq:slicecrit) with the spectral radius
  /// of the truncated M, a posteriori bisection, transfer matrices at recorded heights. Adds one
  /// task per stability slice.
  void rkLayer(const Setup& s, const detail::ModalBasis& sys, const detail::StackLayer& m, double y0, double H, bool record,
               const std::vector<double>& sampleY, std::vector<Task>& tasks) {
    const double a = m.thickness_um();
    std::function<MatC(double)> Mfun;
    if (m.yInvariant()) {
      auto Mconst = std::make_shared<const MatC>(sys.M(m, 0.0));
      Mfun = [Mconst](double) { return *Mconst; };
    } else {
      Mfun = [&s, &sys, &m](double yl) { return sys.M(m, yl / s.k0); };
    }
    double rho = sys.rhoEstimate(m);
    if (m.yInvariant()) {
      rho = std::max(rho, spectralRadius(Mfun(0.0)));
    } else {
      const double ts[5] = {0.02, 0.25, 0.5, 0.75, 0.98};
      double r[5];
      const Settings settings = currentSettings();
      parallelFor(5, sliceThreads(), [&](int i) {
        settings.apply();
        r[i] = spectralRadius(Mfun(ts[i] * a * s.k0));
      });
      for (double v : r) rho = std::max(rho, v);
    }
    const int K = std::max(1, static_cast<int>(std::ceil(rho * H / std::log(growthLimit()))));
    const auto tab = std::make_shared<const ButcherTableau>(ButcherTableau::byName(opt_.rkTableau));
    for (int k = 0; k < K; ++k) {
      // scattered-field formulation: also cut at the background layer boundaries (profile heights)
      const std::vector<double> cut = bgCuts(s, y0 + H * k / K, y0 + H * (k + 1) / K);
      for (std::size_t c = 0; c + 1 < cut.size(); ++c) {
        const double ya = cut[c] - y0, yb = cut[c + 1] - y0;
        tasks.push_back([this, &s, &sys, &m, &sampleY, Mfun, tab, y0, ya, yb, rho, record](std::vector<Element>& out,
                                                                                         Diagnostics& diag) {
          rkSlice(s, sys, m, Mfun, m.yInvariant(), y0, ya, yb, rho, *tab, record, sampleY, out, diag, 0);
        });
      }
    }
  }

  /// Field transfer over [ya, yb] (layer-local normalized y), recording the transfer matrices at
  /// the requested heights; bisects while the growth exceeds growthLimit.
  void rkSlice(const Setup& s, const detail::ModalBasis& sys, const detail::StackLayer& m, const std::function<MatC(double)>& Mfun, bool constant,
               double y0, double ya, double yb, double rho, const ButcherTableau& tab, bool record,
               const std::vector<double>& sampleY, std::vector<Element>& out, Diagnostics& diag, int depth) const {
    if (s.scattered) {
      const bool bloch = s.blochSource;
      auto sigma = [&](double yl, const MatC& Mm) {
        const Setup::BgLayer& L = bgLayer(s, 0.5 * (ya + yb) + y0);
        const MatC MB = bloch ? buildPropagationMatrix(cartesian(m).Q(yl / s.k0), s.alphaB, s.gammaB) : MatC();
        return sourceOf(s, m, Mm, bloch ? &MB : nullptr, L, bgField(L, y0 + yl));
      };
      Element e;
      const Eigen::Index n = sys.ref().Psi.rows();
      MatC T = MatC::Identity(n, n);
      VecC f = VecC::Zero(n);
      double prev = ya;
      auto advance = [&](double to) {  // (T, f) <- (T2 T, T2 f + f2)
        MatC T2;
        VecC f2;
        affineTransfer(Mfun, sigma, prev, to, tab, rho, T2, f2, diag);
        T = T2 * T;
        f = T2 * f + f2;
        prev = to;
      };
      if (record)
        for (double ys : sampleY) {
          const double yl = ys - y0;
          if (yl > ya && yl < yb) {
            advance(yl);
            e.ys.push_back(ys);
            e.T.push_back(T);
            e.f.push_back(f);
          }
        }
      advance(yb);
      if ((!T.allFinite() || T.cwiseAbs().maxCoeff() > growthLimit()) && depth < 12) {
        const double ym = 0.5 * (ya + yb);
        rkSlice(s, sys, m, Mfun, constant, y0, ya, ym, rho, tab, record, sampleY, out, diag, depth + 1);
        rkSlice(s, sys, m, Mfun, constant, y0, ym, yb, rho, tab, record, sampleY, out, diag, depth + 1);
        return;
      }
      FG_REQUIRE(T.allFinite() && f.allFinite(), "RK transfer matrix not finite");
      e.kind = Element::Kind::RK;
      e.y0 = y0 + ya;
      e.y1 = y0 + yb;
      e.model = dynamic_cast<const LayerModel*>(&m);
      e.layerY0 = y0;
      e.S = rkSliceSMatrix(T, sys.ref());
      affineSources(e.S, f, sys.ref(), e.su, e.sv);
      out.push_back(std::move(e));
      return;
    }
    auto transfer = [&](double from, double to) {
      RKStats st;
      MatC T = constant  ? rkTransferConstant(Mfun(0), to - from, tab, rkTol(), rho, st)
               : sys.polarization() ? rkTransfer(Mfun, from, to, tab, rkTol(), rho, st)
                         : rkTransferAuto(Mfun, from, to, tab, rkTol(), rho, st);
      diag.rkSteps += st.steps;
      return T;
    };
    Element e;
    e.kind = Element::Kind::RK;
    e.y0 = y0 + ya;
    e.y1 = y0 + yb;
    e.model = dynamic_cast<const LayerModel*>(&m);
    e.layerY0 = y0;
    const Eigen::Index n = sys.ref().Psi.rows();
    MatC T = MatC::Identity(n, n);
    double prev = ya;
    if (record)
      for (double ys : sampleY) {
        const double yl = ys - y0;
        if (yl > ya && yl < yb) {
          T = transfer(prev, yl) * T;
          prev = yl;
          e.ys.push_back(ys);
          e.T.push_back(T);
        }
      }
    T = transfer(prev, yb) * T;
    if ((!T.allFinite() || T.cwiseAbs().maxCoeff() > growthLimit()) && depth < 12) {
      const double ym = 0.5 * (ya + yb);
      rkSlice(s, sys, m, Mfun, constant, y0, ya, ym, rho, tab, record, sampleY, out, diag, depth + 1);
      rkSlice(s, sys, m, Mfun, constant, y0, ym, yb, rho, tab, record, sampleY, out, diag, depth + 1);
      return;
    }
    FG_REQUIRE(T.allFinite(), "RK transfer matrix not finite");
    e.S = rkSliceSMatrix(T, sys.ref());
    out.push_back(std::move(e));
  }

  /// Modal slice [ya, yb] of a structured layer in the scattered-field formulation: cut at the
  /// background layer boundaries, S and the source amplitudes per piece; a piece whose background
  /// exponent resonates with a mode of M is integrated with RK (warning).
  void scatteredSlices(const Setup& s, const detail::ModalBasis& sys, const detail::StackLayer& m, const MatC& M, double y0,
                       double ya, double yb, double yq, bool record, std::vector<Element>& out, Diagnostics& diag) const {
    const LayerModes modes = sys.modes(M);
    const MatC MB = s.blochSource ? buildPropagationMatrix(cartesian(m).Q(yq), s.alphaB, s.gammaB) : MatC();
    const std::vector<double> cut = bgCuts(s, ya, yb);
    for (std::size_t c = 0; c + 1 < cut.size(); ++c) {
      Element e;
      e.kind = Element::Kind::Modes;
      e.y0 = cut[c];
      e.y1 = cut[c + 1];
      e.model = dynamic_cast<const LayerModel*>(&m);
      e.layerY0 = y0;
      e.qHeight_um = yq;
      RCWASliceRecord rec;
      e.S = stableLayerSMatrix(modes, e.y1 - e.y0, sys.ref(), record ? &rec : nullptr);
      e.modes = std::move(rec.modes);
      e.G = std::move(rec.G);
      if (!modalSources(s, m, M, s.blochSource ? &MB : nullptr, modes, e.S, e.y0, e.y1, sys.ref(), e.su, e.sv,
                        record ? &e : nullptr)) {
        FG_REQUIRE(!record, "scattered field: resonant slice (RK fallback) in a near-field computation");
        diag.warnings.push_back("scattered field: background exponent resonates with a slice mode, slice integrated with RK");
        const ButcherTableau tab = ButcherTableau::byName(opt_.rkTableau);
        const Setup::BgLayer& L = bgLayer(s, 0.5 * (e.y0 + e.y1));
        MatC T;
        VecC f;
        const double rho = std::max(modes.lu.cwiseAbs().maxCoeff(), modes.ld.cwiseAbs().maxCoeff());
        affineTransfer([&M](double) { return M; },
                       [&](double y, const MatC& Mm) { return sourceOf(s, m, Mm, s.blochSource ? &MB : nullptr, L, bgField(L, y)); },
                       e.y0, e.y1, tab, rho, T, f, diag);
        e.kind = Element::Kind::RK;
        e.S = rkSliceSMatrix(T, sys.ref());
        affineSources(e.S, f, sys.ref(), e.su, e.sv);
      }
      out.push_back(std::move(e));
    }
  }

  // ------------------------------------------------------------------ efficiencies
  void efficiencies(const Setup& s, const SMatrix& S, Result& res) const {
    const VecC a_in = incidentAmplitudes(s);
    efficienciesFrom(s, S.S12 * a_in, S.S22 * a_in, res);
  }
  /// Efficiencies from the up-going amplitudes at the top and the down-going ones in the substrate.
  void efficienciesFrom(const Setup& s, VecC b_top, VecC a_sub, Result& res) const {
    FG_PROF_SCOPE("solve.efficiencies");
    const int nf = s.map.size(), nh = 2 * nf;
    const int c = s.map.center();
    const ModeBasis& rp = s.physRef();
    const ModeBasis& sp = s.physSub();
    if (s.transformed) {
      // physical Rayleigh amplitudes from the tangential fields in X (Plan eq:ct_rayleigh)
      b_top = (rp.PsiInv * componentwise(s.toPhys, s.ref.Psi.rightCols(nh) * b_top)).tail(nh);
      a_sub = (sp.PsiInv * componentwise(s.toPhys, s.subModes.Psi.leftCols(nh) * a_sub)).head(nh);
    }
    const VecC a_ph = physicalIncidence(s);
    const auto incOrder = planeWaveOrder(s.epsSup, s.alpha(c), s.gamma(c), rp.beta(c), -1, a_ph(c), a_ph(nf + c));
    const double Sinc = std::abs(incOrder.Sy);
    const Material& sub = *s.sub;
    const bool subLossless = sub.isotropic ? (std::abs(sub.eps.imag()) < 1e-14 && sub.eps.real() > 0)
                                           : (sub.epsTensor - sub.epsTensor.adjoint()).cwiseAbs().maxCoeff() < 1e-14;
    const VecC Fsub = sp.Psi.leftCols(nh) * a_sub;
    auto split = [](OrderEfficiency& o, double e) {
      o.efficiency = e;
      const double ws = std::norm(o.amp_s), wp = std::norm(o.amp_p);
      if (ws + wp > 0) {
        o.eff_s = e * ws / (ws + wp);
        o.eff_p = e * wp / (ws + wp);
      }
    };
    double sumR = 0, sumT = 0;
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      const auto ro = planeWaveOrder(s.epsSup, s.alpha(r), s.gamma(r), rp.beta(r), +1, b_top(r), b_top(nf + r));
      if (ro.propagating) {
        OrderEfficiency o{n, m, 0, ro.amp_s, ro.amp_p};
        split(o, ro.Sy / Sinc);
        res.reflected.push_back(o);
        sumR += ro.Sy / Sinc;
      }
      if (!subLossless) continue;
      if (sub.isotropic) {
        const auto to = planeWaveOrder(sub.eps, s.alpha(r), s.gamma(r), sp.beta(r), -1, a_sub(r), a_sub(nf + r));
        if (to.propagating) {
          OrderEfficiency o{n, m, 0, to.amp_s, to.amp_p};
          split(o, -to.Sy / Sinc);
          res.transmitted.push_back(o);
          sumT += -to.Sy / Sinc;
        }
      } else {
        const double sy = 0.5 * std::real(Fsub(nf + r) * std::conj(Fsub(2 * nf + r)) - Fsub(r) * std::conj(Fsub(3 * nf + r)));
        if (std::abs(sy) > 1e-14 * Sinc) {
          res.transmitted.push_back({n, m, -sy / Sinc, 0, 0});
          sumT += -sy / Sinc;
        }
      }
    }
    res.absorption = 1.0 - sumR - sumT;
    res.diag.energyBalance = sumR + sumT;
  }

  // ------------------------------------------------------------------ aperiodic structures
  /// Physical region [x0, x1] x [z0, z1] (um); 1D: the z integrals are per unit length.
  struct Window {
    double x0, x1, z0, z1;
    bool twoD;
  };
  static Window window(const Stack& stack, const Setup& s) {
    return {s.pmlX ? s.pmlX->physicalBegin() : 0.0, s.pmlX ? s.pmlX->physicalEnd() : stack.dx_um,
            s.pmlZ ? s.pmlZ->physicalBegin() : 0.0, s.pmlZ ? s.pmlZ->physicalEnd() : stack.dz_um, stack.is2D()};
  }
  static cplx axisInt(double kappa, double a, double b) {
    if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
    return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
  }
  /// Field on a plane as coefficient vectors (Ex, Ey, Ez, Hx, Hy, Hz) over lateral wave numbers kx, kz
  /// (1/um). The Bloch-free scattered field of oblique incidence takes the background as one more
  /// "order" (alpha0, gamma0); otherwise the background belongs to the centre order.
  struct PlaneField {
    std::array<VecC, 6> c;
    std::vector<double> kx, kz;
  };
  /// F (4 n_f), E_y, H_y (may be empty: zero) plus sign * background (6 components, may be empty).
  static PlaneField planeField(const Setup& s, const VecC& F, const VecC& Ey, const VecC& Hy, const VecC& bg,
                               double sign = 1.0) {
    const int nf = s.map.size(), c = s.map.center();
    const bool ext = s.blochSource;
    const int n = nf + (ext ? 1 : 0);
    PlaneField P;
    const VecC zero = VecC::Zero(nf);
    const VecC* src[6] = {nullptr, Ey.size() ? &Ey : &zero, nullptr, nullptr, Hy.size() ? &Hy : &zero, nullptr};
    for (int k = 0; k < 6; ++k) {
      P.c[static_cast<std::size_t>(k)] = VecC::Zero(n);
      VecC& v = P.c[static_cast<std::size_t>(k)];
      if (src[k]) v.head(nf) = *src[k];
      else if (F.size()) v.head(nf) = F.segment((k == 0 ? 0 : k == 2 ? 1 : k == 3 ? 2 : 3) * nf, nf);
      if (bg.size()) v(ext ? nf : c) += sign * bg(k);
    }
    for (int r = 0; r < nf; ++r) {
      P.kx.push_back(s.k0 * s.alpha(r).real());
      P.kz.push_back(s.k0 * s.gamma(r).real());
    }
    if (ext) {
      P.kx.push_back(s.k0 * s.inc.alpha0);
      P.kz.push_back(s.k0 * s.inc.gamma0);
    }
    return P;
  }
  /// Upward power through y = const over the window: 1/2 Re(Hx^H W Ez - Hz^H W Ex),
  /// W(r', r) = int exp(i ((kx_r - kx_r') x + (kz_r - kz_r') z)).
  static double planeFlux(const Window& w, const PlaneField& P) {
    const auto n = static_cast<Eigen::Index>(P.kx.size());
    MatC W(n, n);
    for (Eigen::Index r = 0; r < n; ++r)
      for (Eigen::Index c = 0; c < n; ++c) {
        cplx v = axisInt(P.kx[static_cast<std::size_t>(c)] - P.kx[static_cast<std::size_t>(r)], w.x0, w.x1);
        if (w.twoD) v *= axisInt(P.kz[static_cast<std::size_t>(c)] - P.kz[static_cast<std::size_t>(r)], w.z0, w.z1);
        W(r, c) = v;
      }
    return 0.5 * std::real(P.c[3].dot(W * P.c[2]) - P.c[5].dot(W * P.c[0]));  // dot conjugates its first argument
  }
  /// Power in +x through x = xs (axis 0) or in +z through z = zs (axis 1), integrated over the other
  /// lateral axis of the window: S_x = 1/2 Re(Ey Hz* - Ez Hy*), S_z = 1/2 Re(Ex Hy* - Ey Hx*).
  static double sideFlux(const Window& w, const PlaneField& P, int axis, double pos) {
    const auto n = static_cast<Eigen::Index>(P.kx.size());
    MatC W(n, n);
    for (Eigen::Index r = 0; r < n; ++r)
      for (Eigen::Index c = 0; c < n; ++c) {
        const auto ur = static_cast<std::size_t>(r), uc = static_cast<std::size_t>(c);
        W(r, c) = axis == 0 ? std::exp(I * ((P.kx[uc] - P.kx[ur]) * pos)) *
                                  (w.twoD ? axisInt(P.kz[uc] - P.kz[ur], w.z0, w.z1) : cplx(1.0))
                            : std::exp(I * ((P.kz[uc] - P.kz[ur]) * pos)) * axisInt(P.kx[uc] - P.kx[ur], w.x0, w.x1);
      }
    return axis == 0 ? 0.5 * std::real(P.c[5].dot(W * P.c[1]) - P.c[4].dot(W * P.c[2]))
                     : 0.5 * std::real(P.c[4].dot(W * P.c[0]) - P.c[3].dot(W * P.c[1]));
  }
  /// Outward power through the lateral sides of the window (PML axes).
  static double sidesFlux(const Setup& s, const Window& w, const PlaneField& P) {
    double out = 0;
    if (s.pmlX) out += sideFlux(w, P, 0, w.x1) - sideFlux(w, P, 0, w.x0);
    if (s.pmlZ) out += sideFlux(w, P, 1, w.z1) - sideFlux(w, P, 1, w.z0);
    return out;
  }

  /// Heights (um) and weights of the lateral-side quadrature: Gauss nodes per layer, and in a 1D profile
  /// layer also split at its height in the zones (E_y, E_z jump on the sides there).
  static void sideNodes(const Stack& stack, const Setup& s, const detail::LayerModels& models, std::vector<double>& ys,
                        std::vector<double>& wy) {
    std::vector<double> br{0.0};
    for (std::size_t k = 0; k < stack.layers.size(); ++k) {
      if (const auto* p = std::get_if<Profile1DLayer>(&stack.layers[k].spec)) {
        const double g0 = profileValue(*p, 0.0, stack.dx_um);
        if (g0 > 0 && g0 < p->depth_um) br.push_back(br.back() + g0);
      }
      br.push_back(br.back() + models[k]->thickness_um());
    }
    const auto& g = quadrature::gaussLegendre(8);
    const double hmax = 0.1 * 2 * PI / s.k0;  // panels <= lambda / 10
    for (std::size_t k = 0; k + 1 < br.size(); ++k) {
      const double a = br[k], b = br[k + 1];
      const int m = std::max(1, static_cast<int>(std::ceil((b - a) / hmax)));
      for (int p = 0; p < m; ++p)
        for (std::size_t q = 0; q < g.x.size(); ++q) {
          ys.push_back(a + (b - a) * (p + 0.5 * (1 + g.x[q])) / m);
          wy.push_back(0.5 * (b - a) / m * g.w[q]);
        }
    }
  }

  /// Scattering outputs (Plan section "M18: Ausgaben"), total-field formulation (normal incidence): the
  /// background stack (same transformation, same N, so the discretization errors cancel) gives E_bg;
  /// fluxes over the physical region. The absorption is the inward total-field flux through the closed
  /// box physical region x [0, top] (top, bottom and the lateral sides), minus that of the background.
  void scattering(const Stack& stack, const Setup& s, const detail::LayerModels& models, Result& res) {
    FG_PROF_SCOPE("solve.scattering");
    const Stack bgStack = backgroundStack(stack);
    const auto bgModels = buildModels(bgStack, s);
    const SMatrix S = assemble(s, CartesianBasis(s), layerPtrs(models), false, {}, res.diag).total;
    const SMatrix Sb = assemble(s, CartesianBasis(s), layerPtrs(bgModels), false, {}, res.diag).total;
    const int nf = s.map.size(), nh = 2 * nf;
    const VecC a_in = incidentAmplitudes(s);
    const VecC b = S.S12 * a_in, bb = Sb.S12 * a_in, t = S.S22 * a_in, tb = Sb.S22 * a_in;
    const Window w = window(stack, s);
    const VecC none;
    auto flux = [&](const VecC& F) { return planeFlux(w, planeField(s, F, none, none, none)); };
    const MatC& up = s.ref.Psi;
    const auto down = s.subModes.Psi.leftCols(nh);
    const double area = (w.x1 - w.x0) * (w.twoD ? w.z1 - w.z0 : 1.0);
    const double Sinc = std::abs(flux(up.leftCols(nh) * a_in)) / area;  // incident intensity (normal incidence)
    VecC top(2 * nh), topb(2 * nh);
    top << a_in, b;
    topb << a_in, bb;
    ScatteringResult& sc = res.scattering;
    sc.present = true;
    sc.up = flux(up.rightCols(nh) * (b - bb)) / Sinc;
    sc.down = -flux(down * (t - tb)) / Sinc;
    // lateral sides: the background has no lateral flux (laterally uniform)
    std::vector<double> ys, wy;
    sideNodes(stack, s, models, ys, wy);
    double side = 0, sideSca = 0;
    if (!ys.empty()) {
      const FieldSampleSet fb = samples(s, bgModels, ys, nullptr);
      const FieldSampleSet fs = samples(s, models, ys, nullptr);
      for (std::size_t i = 0; i < ys.size(); ++i) {
        const FieldSample &f = fs.samples[i], &g = fb.samples[i];
        side += wy[i] * sidesFlux(s, w, planeField(s, f.F, f.Ey, f.Hy, none));
        sideSca += wy[i] * sidesFlux(s, w, planeField(s, f.F - g.F, f.Ey - g.Ey, f.Hy - g.Hy, none));
      }
    }
    sc.lateral = sideSca / Sinc;
    sc.absorbed = ((flux(down * t) - flux(up * top) - side) - (flux(down * tb) - flux(up * topb))) / Sinc;
    res.diag.energyBalance = 1.0;  // no diffraction orders; the M16 safeguard does not apply
  }

  /// Scattering outputs of the scattered-field formulation: field samples just above the stack, just
  /// below it and on the lateral sides; total = scattered + background, the absorption is the inward
  /// total-field flux through the closed box minus that of the background alone.
  void scatteredOutputs(const Stack& stack, const Setup& s, const detail::LayerModels& models, Result& res) {
    FG_PROF_SCOPE("solve.scattering");
    double top_um = 0;
    for (const auto& m : models) top_um += m->thickness_um();
    std::vector<double> ys{top_um + 1e-9, -1e-9}, wy{0.0, 0.0};
    sideNodes(stack, s, models, ys, wy);
    const FieldSampleSet fs = samples(s, models, ys, nullptr);
    const Window w = window(stack, s);
    const VecC none;
    // total, background and scattered field of sample i
    auto tot = [&](const FieldSample& f) {
      return planeField(s, f.F, f.Ey, f.Hy, fs.backgroundIncluded ? none : f.bg);
    };
    auto bgo = [&](const FieldSample& f) { return planeField(s, none, none, none, f.bg); };
    auto sca = [&](const FieldSample& f) {
      return planeField(s, f.F, f.Ey, f.Hy, fs.backgroundIncluded ? f.bg : none, -1.0);
    };
    const FieldSample &ft = fs.samples[0], &fb = fs.samples[1];
    double sideTot = 0, sideBg = 0, sideSca = 0;
    for (std::size_t i = 2; i < ys.size(); ++i) {
      sideTot += wy[i] * sidesFlux(s, w, tot(fs.samples[i]));
      sideBg += wy[i] * sidesFlux(s, w, bgo(fs.samples[i]));
      sideSca += wy[i] * sidesFlux(s, w, sca(fs.samples[i]));
    }
    const double I0 = s.incIntensity;
    ScatteringResult& sc = res.scattering;
    sc.present = true;
    sc.up = planeFlux(w, sca(ft)) / I0;
    sc.down = -planeFlux(w, sca(fb)) / I0;
    sc.lateral = sideSca / I0;
    const double netTot = planeFlux(w, tot(fb)) - planeFlux(w, tot(ft)) - sideTot;
    const double netBg = planeFlux(w, bgo(fb)) - planeFlux(w, bgo(ft)) - sideBg;
    sc.absorbed = (netTot - netBg) / I0;
    res.diag.energyBalance = 1.0;
  }

  /// Guided modes of a half-space basis among the columns [first, first + n_h): n_eff = |Im lambda| above the
  /// cladding index, |Re lambda| < 1e-3 n_eff (the PML perturbs guided modes by up to 1e-6, radiation modes
  /// are strongly damped) and >= 90 % of int |E|^2 dX inside the physical region (PML modes with nearly real
  /// lambda live in the zones); sorted by descending n_eff.
  static std::vector<int> guidedModes(const Stack& stack, const Setup& s, const ModeBasis& b, int first, int nh,
                                      double nclad) {
    const Window w = window(stack, s);
    const VecC none;
    const double cell = stack.dx_um * (stack.is2D() ? stack.dz_um : 1.0);
    auto inside = [&](const VecC& F) {  // int_window |E_t|^2 / int_cell |E_t|^2 (Parseval in X)
      PlaneField P = planeField(s, F, none, none, none);
      const auto n = static_cast<Eigen::Index>(P.kx.size());
      MatC W(n, n);
      for (Eigen::Index r = 0; r < n; ++r)
        for (Eigen::Index c = 0; c < n; ++c) {
          cplx v = axisInt(P.kx[static_cast<std::size_t>(c)] - P.kx[static_cast<std::size_t>(r)], w.x0, w.x1);
          if (w.twoD) v *= axisInt(P.kz[static_cast<std::size_t>(c)] - P.kz[static_cast<std::size_t>(r)], w.z0, w.z1);
          W(r, c) = v;
        }
      const double in = std::real(P.c[0].dot(W * P.c[0]) + P.c[2].dot(W * P.c[2]));
      return in / (cell * (P.c[0].squaredNorm() + P.c[2].squaredNorm()));
    };
    std::vector<std::pair<double, int>> g;
    for (int k = first; k < first + nh; ++k) {
      const cplx l = b.lambda(k);
      if (std::abs(l.real()) <= 1e-3 * std::abs(l.imag()) && std::abs(l.imag()) > nclad * (1 + 1e-9) &&
          inside(b.Psi.col(k)) >= 0.9)
        g.push_back({std::abs(l.imag()), k});
    }
    std::sort(g.begin(), g.end(), [](const auto& a, const auto& c) { return a.first > c.first; });
    std::vector<int> out;
    for (const auto& [n, k] : g) out.push_back(k);
    return out;
  }

  /// Modal excitation (Plan "Einfall (b)", "Ausgaben"): a guided mode of the structured superstrate comes
  /// in; reflected guided modes (superstrate) and transmitted ones (structured substrate), powers from the
  /// flux over the physical region, the rest is radiated.
  void modalOutputs(const Stack& stack, const Setup& s, const detail::LayerModels& models, Result& res) {
    FG_PROF_SCOPE("solve.modes");
    const AperiodicSpec& ap = stack.aperiodic;
    const SMatrix S = assemble(s, CartesianBasis(s), layerPtrs(models), false, {}, res.diag).total;
    const int nf = s.map.size(), nh = 2 * nf;
    const Window w = window(stack, s);
    const VecC none;
    auto up = [&](const VecC& F) { return planeFlux(w, planeField(s, F, none, none, none)); };
    auto nclad = [&](const std::string& name) { return std::sqrt(std::max(stack.material(name).tensor()(0, 0).real(), 0.0)); };
    auto describe = [&](const ModeBasis& b, int col, int index) {
      ModeCoefficient m;
      m.index = index;
      m.nEff = cplx(std::abs(b.lambda(col).imag()), std::abs(b.lambda(col).real()));  // n_eff + i damping
      const VecC E = b.Psi.col(col);
      const double ex = E.segment(0, nf).squaredNorm(), ez = E.segment(nf, nf).squaredNorm();
      m.teFraction = ex + ez > 0 ? ez / (ex + ez) : 0.0;
      return m;
    };
    const std::vector<int> down = guidedModes(stack, s, s.ref, 0, nh, nclad(stack.superstrate));
    FG_REQUIRE(!down.empty(), "modal excitation: the structured superstrate has no guided mode");
    int j = 0;
    if (ap.nEff > 0) {
      for (int k = 1; k < static_cast<int>(down.size()); ++k)
        if (std::abs(std::abs(s.ref.lambda(down[static_cast<std::size_t>(k)]).imag()) - ap.nEff) <
            std::abs(std::abs(s.ref.lambda(down[static_cast<std::size_t>(j)]).imag()) - ap.nEff))
          j = k;
    } else {
      FG_REQUIRE(ap.mode >= 0 && ap.mode < static_cast<int>(down.size()),
                 "modal excitation: mode " + std::to_string(ap.mode) + " not guided (" + std::to_string(down.size()) +
                     " guided modes)");
      j = ap.mode;
    }
    const int cin = down[static_cast<std::size_t>(j)];
    VecC a_in = VecC::Zero(nh);
    a_in(cin) = 1.0;
    const double Pin = -up(s.ref.Psi.col(cin));  // down-going power
    FG_REQUIRE(Pin > 0, "modal excitation: the incident mode carries no power through the physical region");
    ModalResult& mr = res.modal;
    mr.present = true;
    mr.incident = describe(s.ref, cin, j);
    mr.incident.amplitude = 1.0;
    mr.incident.power = 1.0;
    const VecC b = S.S12 * a_in, t = S.S22 * a_in;
    double sum = 0;
    const std::vector<int> rUp = guidedModes(stack, s, s.ref, nh, nh, nclad(stack.superstrate));
    for (std::size_t k = 0; k < rUp.size(); ++k) {
      ModeCoefficient m = describe(s.ref, rUp[k], static_cast<int>(k));
      const double Pk = up(s.ref.Psi.col(rUp[k]));
      m.amplitude = b(rUp[k] - nh) * std::sqrt(std::abs(Pk) / Pin);
      m.power = std::norm(b(rUp[k] - nh)) * Pk / Pin;
      sum += m.power;
      mr.reflected.push_back(m);
    }
    if (stack.substrateStructure) {
      const std::vector<int> tDown = guidedModes(stack, s, s.subModes, 0, nh, nclad(stack.substrate));
      for (std::size_t k = 0; k < tDown.size(); ++k) {
        ModeCoefficient m = describe(s.subModes, tDown[k], static_cast<int>(k));
        const double Pk = -up(s.subModes.Psi.col(tDown[k]));
        m.amplitude = t(tDown[k]) * std::sqrt(std::abs(Pk) / Pin);
        m.power = std::norm(t(tDown[k])) * Pk / Pin;
        sum += m.power;
        mr.transmitted.push_back(m);
      }
    }
    mr.radiated = 1.0 - sum;
    res.diag.energyBalance = 1.0;
  }

  /// Amplitudes (a_x, a_z) of the incident wave in the analytic basis of the physical orders.
  static VecC physicalIncidence(const Setup& s) {
    const int nf = s.map.size();
    VecC a_in = VecC::Zero(2 * nf);
    a_in(s.map.center()) = s.inc.A(0);
    a_in(nf + s.map.center()) = s.inc.A(2);
    return a_in;
  }
  /// Incident amplitudes in the reference basis; transformed coordinates: the incident field expanded
  /// in X (Plan eq:ct_inc) and projected onto the down-going numerical modes.
  static VecC incidentAmplitudes(const Setup& s) {
    const VecC a = physicalIncidence(s);
    if (!s.transformed) return a;
    const Eigen::Index nh = a.size();
    VecC b = (s.ref.PsiInv * componentwise(s.fromPhys, s.refPhys.Psi.leftCols(nh) * a)).head(nh);
    // a field incident from y -> infinity contains no modes growing upwards: the projection onto the
    // evanescent down modes is truncation noise, which exp(lambda y) amplifies above the stack
    for (Eigen::Index c = 0; c < nh; ++c)
      if (s.ref.lambda(c).real() > 1e-8) b(c) = 0;
    return b;
  }

  // ------------------------------------------------------------------ field coefficients
  FieldSampleSet samples(const Setup& s, const std::vector<std::unique_ptr<LayerModel>>& models,
                         const std::vector<double>& y_um, std::vector<const Element*>* where) {
    FG_PROF_SCOPE("nearfield.coefficients");
    std::vector<double> yn;
    for (double y : y_um) yn.push_back(y * s.k0);
    std::vector<double> sorted = yn;
    std::sort(sorted.begin(), sorted.end());
    Diagnostics diag;
    assembly_ = assemble(s, CartesianBasis(s), layerPtrs(models), true, sorted, diag);
    const auto& el = assembly_.elems;
    const std::size_t K = el.size() - 1;
    const int nf = s.map.size(), nh = 2 * nf;
    const VecC b_in = VecC::Zero(nh);
    std::vector<VecC> amp(K + 1);  // (a_q; b_q) in the reference basis (Plan eq:internal)
    VecC a_in, a_sub;
    if (s.scattered) {
      // scattered field: nothing comes in, the sources of the elements drive it
      auto src = [&](std::size_t q) {
        const Element& e = el[q];
        return SourcedSMatrix{e.S, e.su.size() ? e.su : VecC(VecC::Zero(nh)), e.sv.size() ? e.sv : VecC(VecC::Zero(nh))};
      };
      std::vector<SourcedSMatrix> SL(K + 1), SU(K + 1);
      SL[0] = src(0);
      for (std::size_t q = 1; q <= K; ++q) SL[q] = star(SL[q - 1], src(q));
      SU[K] = SourcedSMatrix::identity(nh);
      for (std::size_t q = K; q-- > 0;) SU[q] = star(src(q + 1), SU[q + 1]);
      a_in = VecC::Zero(nh);
      for (std::size_t q = 0; q <= K; ++q) amp[q] = interfaceAmplitudes(SL[q], SU[q], b_in, a_in);
      a_sub = SL[K].v;
    } else {
      // cumulative S-matrices below (S_L) and above (S_U) every interface q (top of element q)
      std::vector<SMatrix> SL(K + 1), SU(K + 1);
      SL[0] = el[0].S;
      for (std::size_t q = 1; q <= K; ++q) SL[q] = star(SL[q - 1], el[q].S);
      SU[K] = SMatrix::identity(nh);
      for (std::size_t q = K; q-- > 0;) SU[q] = star(el[q + 1].S, SU[q + 1]);
      a_in = incidentAmplitudes(s);
      for (std::size_t q = 0; q <= K; ++q) amp[q] = interfaceAmplitudes(SL[q], SU[q], b_in, a_in);
      a_sub = SL[K].S22 * a_in;  // substrate amplitudes at y = 0 (total = SL[K])
    }

    FieldSampleSet out;
    out.k0 = s.k0;
    out.alpha = s.alpha;
    out.gamma = s.gamma;
    out.alpha0 = s.inc.alpha0;
    out.gamma0 = s.inc.gamma0;
    out.backgroundIncluded = s.scattered && !s.blochSource;
    if (where) where->assign(y_um.size(), nullptr);
    for (std::size_t i = 0; i < yn.size(); ++i) {
      const double y = yn[i];
      FieldSample fsm;
      fsm.y_um = y_um[i];
      QBlocks Q;
      if (y < 0) {
        fsm.region = -1;
        VecC ph(nh);
        for (int c = 0; c < nh; ++c) ph(c) = std::exp(s.subModes.lambda(c) * y);
        fsm.F = s.subModes.Psi.leftCols(nh) * ph.cwiseProduct(a_sub);
        Q = s.subQ                ? *s.subQ
            : s.sub->isotropic ? homogeneousIsotropicQ(s.sub->eps, nf)
                               : homogeneousTensorQ(s.sub->epsTensor, nf);
      } else if (y > assembly_.yTop) {
        fsm.region = 1;
        const double d = y - assembly_.yTop;
        VecC v(2 * nh);
        for (int c = 0; c < nh; ++c) {
          v(c) = a_in(c) * std::exp(s.ref.lambda(c) * d);
          v(nh + c) = amp[K](nh + c) * std::exp(s.ref.lambda(nh + c) * d);
        }
        fsm.F = s.ref.Psi * v;
        Q = s.supQ ? *s.supQ : homogeneousIsotropicQ(s.epsSup, nf);
      } else {
        std::size_t k = 1;
        while (k < K && y > el[k].y1) ++k;
        const Element& e = el[k];
        if (where) (*where)[i] = &e;
        const VecC& lo = amp[k - 1];  // (a_{k-1}; b_{k-1}) at the bottom
        const VecC& hi = amp[k];      // (a_k; b_k) at the top
        if (e.kind == Element::Kind::Modes) {
          const Eigen::Index n = e.modes.lu.size();
          VecC rhs(2 * n);
          rhs << lo.tail(n), hi.head(n);  // (b_bottom; a_top)
          if (!e.pVec.empty()) {  // homogeneous part = total - particular solution
            const VecC g0 = s.ref.PsiInv * e.particular(e.y0), g1 = s.ref.PsiInv * e.particular(e.y1);
            rhs.head(n) -= g0.tail(n);
            rhs.tail(n) -= g1.head(n);
          }
          const VecC cc = e.G.partialPivLu().solve(rhs);
          VecC cu = cc.head(n), cd = cc.tail(n);
          for (Eigen::Index j = 0; j < n; ++j) {
            cu(j) *= std::exp(e.modes.lu(j) * (y - e.y0));
            cd(j) *= std::exp(e.modes.ld(j) * (y - e.y1));
          }
          fsm.F = e.modes.Wu * cu + e.modes.Wd * cd;
          if (!e.pVec.empty()) fsm.F += e.particular(y);
          Q = e.model->Q(e.qHeight_um >= 0 ? e.qHeight_um : 0.0);
        } else {
          const VecC F0 = s.ref.Psi * lo;
          fsm.F = F0;
          for (std::size_t j = 0; j < e.ys.size(); ++j)
            if (std::abs(e.ys[j] - y) <= 1e-14 * std::max(1.0, std::abs(y)))
              fsm.F = e.f.empty() ? VecC(e.T[j] * F0) : VecC(e.T[j] * F0 + e.f[j]);
          if (y >= e.y1 - 1e-14 * std::max(1.0, std::abs(y))) fsm.F = s.ref.Psi * hi;
          Q = e.model->Q((y - e.layerY0) / s.k0);
        }
      }
      VecC eyCorr;
      if (s.scattered) {
        // background field (incidence order): superstrate incident + reflected, substrate transmitted
        const VecC a0 = VecC::Constant(1, s.inc.alpha0), g0 = VecC::Constant(1, s.inc.gamma0);
        VecC Fb(4);
        const QBlocks* Qb1 = nullptr;
        QBlocks Qh;
        if (fsm.region > 0) {
          const ModeBasis sup = isotropicModes(s.epsSup, a0, g0);
          const double d = y - assembly_.yTop;
          VecC v(4);
          v << s.inc.A(0) * std::exp(sup.lambda(0) * d), s.inc.A(2) * std::exp(sup.lambda(1) * d),
              s.bgTop(0) * std::exp(sup.lambda(2) * d), s.bgTop(1) * std::exp(sup.lambda(3) * d);
          Fb = sup.Psi * v;
          Qh = homogeneousIsotropicQ(s.epsSup, 1);
          Qb1 = &Qh;
        } else if (fsm.region < 0) {
          const ModeBasis sb = s.sub->isotropic ? isotropicModes(s.sub->eps, a0, g0) : anisotropicModes(s.sub->epsTensor, a0, g0);
          VecC v = VecC::Zero(4);
          v(0) = s.bgSub(0) * std::exp(sb.lambda(0) * y);
          v(1) = s.bgSub(1) * std::exp(sb.lambda(1) * y);
          Fb = sb.Psi * v;
          Qh = s.sub->isotropic ? homogeneousIsotropicQ(s.sub->eps, 1) : homogeneousTensorQ(s.sub->epsTensor, 1);
          Qb1 = &Qh;
        } else {
          const Setup::BgLayer& L = bgLayer(s, y);
          Fb = bgField(L, y);
          Qb1 = &L.Q1;
          if (s.blochSource) {
            // E_y of the scattered field: Q_yy^-1 (... - P_y) = E_y(Q, F_s) + [E_y(Q, F_bg) - E_y(Q_bg, F_bg)]
            VecC FB = VecC::Zero(4 * nf);
            for (int c = 0; c < 4; ++c) FB(c * nf + s.map.center()) = Fb(c);
            eyCorr = s.Tsrc * (computeEy(Q, s.alphaB, s.gammaB, FB) - computeEy(L.Q, s.alphaB, s.gammaB, FB));
          }
        }
        fsm.bg.resize(6);
        fsm.bg << Fb(0), computeEy(*Qb1, a0, g0, Fb)(0), Fb(1), Fb(2), computeHy(a0, g0, Fb)(0), Fb(3);
        if (!s.blochSource)
          for (int c = 0; c < 4; ++c) fsm.F(c * nf + s.map.center()) += Fb(c);  // total field
      }
      fsm.Ey = s.transformed ? computeEy(Q, s.A, s.G, fsm.F) : computeEy(Q, s.alpha, s.gamma, fsm.F);
      if (eyCorr.size()) fsm.Ey += eyCorr;
      fsm.Hy = s.transformed ? computeHy(s.A, s.G, fsm.F) : computeHy(s.alpha, s.gamma, fsm.F);
      out.samples.push_back(std::move(fsm));
    }
    return out;
  }

  // ------------------------------------------------------------------ synthesis
  /// Real-space field on the points (x, z) of every sample i0..i1-1 (Plan section "Nahfeld");
  /// row-major with the outer index running over samples (vertical cut) or z (horizontal cut).
  FieldCut synthesize(const Stack& stack, const Setup& s, const FieldSampleSet& fs,
                      const std::vector<const Element*>& where, const std::vector<double>& xs,
                      const std::vector<double>& zs, int i0, int i1, const NearFieldRequest& req,
                      const std::string& name) const {
    FieldCut cut;
    cut.cut = name;
    cut.lambda_um = 2 * PI / s.k0;
    cut.axis1_um = xs;
    const int nf = s.map.size();
    for (int i = i0; i < i1; ++i) cut.axis2_um.push_back(fs.samples[static_cast<std::size_t>(i)].y_um);
    // phase factors exp(i alpha_r X), exp(i gamma_r Z) at X = f^{-1}(x), Z = g^{-1}(z) (Plan section "M17")
    std::vector<VecC> px, pz;
    for (double x_um : xs) {
      const double X = s.tr.x ? s.tr.x->inverse(x_um) : x_um;
      VecC v(nf);
      for (int r = 0; r < nf; ++r) v(r) = std::exp(I * (s.alpha(r) * X * s.k0));
      px.push_back(v);
    }
    for (double z_um : zs) {
      const double Z = s.tr.z ? s.tr.z->inverse(z_um) : z_um;
      VecC v(nf);
      for (int r = 0; r < nf; ++r) v(r) = std::exp(I * (s.gamma(r) * Z * s.k0));
      pz.push_back(v);
    }
    for (int i = i0; i < i1; ++i) {
      const FieldSample& f = fs.samples[static_cast<std::size_t>(i)];
      const Element* e = where[static_cast<std::size_t>(i)];
      // the Bloch-free scattered field is not D = Q E (sources): direct reconstruction
      const bool fff = req.reconstruction == Reconstruction::FFF && e && e->model && !e->model->homogeneous() &&
                       !(f.bg.size() && !fs.backgroundIncluded);
      // scattered-field formulation: add (total) or remove (scattered) the background of the incidence order
      const double bgSign = f.bg.size() == 0 ? 0.0
                            : fs.backgroundIncluded ? (req.scattered ? -1.0 : 0.0)
                                                    : (req.scattered ? 0.0 : 1.0);
      // Fourier vectors of E and H
      VecC E[3] = {f.F.segment(0, nf), f.Ey, f.F.segment(nf, nf)};
      VecC Hc[3] = {f.F.segment(2 * nf, nf), f.Hy, f.F.segment(3 * nf, nf)};
      VecC et[3], dn[3];
      const double yq = e ? (e->qHeight_um >= 0 ? e->qHeight_um : f.y_um - e->layerY0 / s.k0) : 0.0;
      if (fff) {
        // FFF-consistent reconstruction (Plan eq:recon): e_t = P_T E, d_n = P_N D continuous
        const QBlocks Q = e->model->Q(yq);
        VecC D[3];
        for (int a = 0; a < 3; ++a) D[a] = Q(a, 0) * E[0] + Q(a, 1) * E[1] + Q(a, 2) * E[2];
        const auto& nn = e->model->NN(yq);
        for (int a = 0; a < 3; ++a) {
          dn[a] = VecC::Zero(nf);
          et[a] = E[a];
          for (int b = 0; b < 3; ++b) {
            const auto& g = nn.NN[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            if (g.Px() < 2 * s.map.Nx() || g.Pz() < 2 * s.map.Nz()) continue;  // zero entry
            const MatC P = fourier::toeplitz(s.map, g);
            dn[a] += P * D[b];
            et[a] -= P * E[b];
          }
        }
      }
      for (std::size_t iz = 0; iz < zs.size(); ++iz)
        for (std::size_t ix = 0; ix < xs.size(); ++ix) {
          const VecC ph = px[ix].cwiseProduct(pz[iz]);
          const double x_um = xs[ix], z_um = zs[iz];
          for (int c = 0; c < 3; ++c) cut.H[c].push_back(ph.transpose() * Hc[c]);
          const cplx bph = bgSign == 0.0 ? cplx(0.0) : bgSign * std::exp(I * (s.k0 * (fs.alpha0 * x_um + fs.gamma0 * z_um)));
          if (bgSign != 0.0)
            for (int c = 0; c < 3; ++c) cut.H[c].back() += bph * f.bg(3 + c);
          if (fff) {
            Vec3C etx, dnx;
            for (int a = 0; a < 3; ++a) {
              etx(a) = ph.transpose() * et[a];
              dnx(a) = ph.transpose() * dn[a];
            }
            const Mat3C eps = e->model->epsAt(x_um, yq, z_um);
            Vec3C Ex;
            if (e->model->rank1Projector()) {
              const ProjectorSet ps = projectorMatrices(eps, e->model->normalAt(x_um, yq, z_um));
              Ex = ps.H1 * etx + ps.invXi * dnx;
            } else {
              // isotropic, general projector field P: E = e_t + d_n / eps holds for any P
              Ex = etx + dnx / eps(0, 0);
            }
            for (int c = 0; c < 3; ++c) cut.E[c].push_back(Ex(c));
          } else {
            for (int c = 0; c < 3; ++c) cut.E[c].push_back(ph.transpose() * E[c]);
          }
          if (bgSign != 0.0)
            for (int c = 0; c < 3; ++c) cut.E[c].back() += bph * f.bg(c);
          cut.material.push_back(materialAt(stack, f, e, x_um, yq, z_um).real());
        }
    }
    return cut;
  }

  static cplx materialAt(const Stack& stack, const FieldSample& f, const Element* e, double x_um, double yq,
                         double z_um) {
    if (f.region < 0) return stack.material(stack.substrate).tensor()(0, 0);
    if (f.region > 0 || !e || !e->model) return stack.material(stack.superstrate).tensor()(0, 0);
    return e->model->epsAt(x_um, yq, z_um)(0, 0);
  }

  SolverOptions opt_;
  Assembly assembly_;  // last recorded assembly (referenced by 'where' pointers)
};

}  // namespace

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt) { return std::make_unique<StackSolver>(opt); }

namespace detail {

linalg::Precision resolvePrecision(const SolverOptions& opt) {
  std::string p = opt.precision;
  if (const char* env = std::getenv("FG_PRECISION"); env && *env) p = env;
  if (p != "auto") return linalg::parsePrecision(p);
  // auto (Plan section "M16: Standardwahl"): CPU always FP64; GPU by device policy
  const Eigen::Index n = 4 * static_cast<Eigen::Index>(2 * opt.Nx + 1) * (2 * opt.Nz + 1);
  return linalg::gpuSelected(n) ? linalg::gpuAutoPrecision() : linalg::Precision::FP64;
}

void applyThreadSettings(const SolverOptions& opt) {
  linalg::setPreference(opt.backend, opt.backendThreshold);  // thread-local (M10)
  fourier::setFFTSampling(opt.fftSampling);                   // thread-local
  linalg::setPrecision(resolvePrecision(opt));                // thread-local (M16)
}

Result solveWithModels(const SolverOptions& opt, const Stack& stack, const Wavelength& wl, const LayerModels& models) {
  return StackSolver(opt).solveWithModels(stack, wl, models);
}

}  // namespace detail

ThreadPlan threadPlan(const SolverOptions& opt, std::size_t scanPoints, int cores) {
  ThreadPlan p;
  cores = std::max(1, cores);
  p.blas = std::clamp(opt.threadsBlas, 1, cores);
  const int avail = std::max(1, cores / p.blas);
  const int points = static_cast<int>(std::clamp<std::size_t>(scanPoints, 1, 1 << 30));
  p.scan = std::min(opt.threadsScan > 0 ? std::min(opt.threadsScan, avail) : avail, points);
  const int rest = std::max(1, avail / p.scan);
  p.slices = opt.threadsSlices > 0 ? std::min(opt.threadsSlices, rest) : rest;
  return p;
}

std::vector<Result> solveScan(const SolverOptions& opt, const Stack& stack, const std::vector<Wavelength>& scan) {
  FG_PROF_SCOPE("scan");
  const ThreadPlan plan = threadPlan(opt, scan.size(), availableThreads());
  SolverOptions o = opt;
  o.threadsScan = plan.scan;
  o.threadsSlices = plan.slices;
  o.threadsBlas = plan.blas;
  setEigenThreads(plan);
  std::vector<Result> out(scan.size());
  // one solver per point (with the CUDA backend: one stream per thread, M10)
  parallelFor(static_cast<int>(scan.size()), plan.scan, [&](int i) {
    out[static_cast<std::size_t>(i)] = makeSolver(o)->solve(stack, scan[static_cast<std::size_t>(i)]);
  });
  return out;
}

}  // namespace fg
