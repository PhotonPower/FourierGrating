// FourierGrating - stack assembly shared by the RK and RCWA solvers, efficiencies and near fields
// (Plan sections "S-Matrix-Formalismus", "Loeser A/B", "Nahfeld").
#include <algorithm>
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
#include "fg/physics/modes.hpp"
#include "fg/solvers/rcwa.hpp"
#include "fg/solvers/resident.hpp"
#include "fg/solvers/rk.hpp"
#include "fg/solvers/solver.hpp"

#include "layer_model.hpp"

namespace fg {

namespace {

using detail::LayerModel;

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
  const ModeBasis& physRef() const { return transformed ? refPhys : ref; }
  const ModeBasis& physSub() const { return transformed ? subPhys : subModes; }
};

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
  const std::vector<Eigen::Index> idx = modeOrder(lam, [&](Eigen::Index k) { return modeFlux(W.col(k)); });
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

/// Field system of an assembly: the full 4 n_f problem, or one polarization of the decoupled 1D
/// in-plane problem (2 n_f, Plan eq:TE / eq:TM).
struct System {
  std::optional<Polarization> pol;
  ModeBasis ref, sub;

  MatC M(const QBlocks& Q, const Setup& s) const {
    if (pol) return buildPropagationMatrix(Q, s.alpha, *pol);
    return s.transformed ? buildPropagationMatrix(Q, s.A, s.G) : buildPropagationMatrix(Q, s.alpha, s.gamma);
  }
  LayerModes modes(const MatC& M) const { return pol ? layerModesPolarization(M, *pol) : layerModesAuto(M); }
  ModeBasis basis(const ModeBasis& full) const { return pol ? polarizationBasis(full, *pol) : full; }
};

System fullSystem(const Setup& s) { return {std::nullopt, s.ref, s.subModes}; }
System polarizationSystem(const Setup& s, Polarization p) {
  return {p, polarizationBasis(s.ref, p), polarizationBasis(s.subModes, p)};
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
};

struct Assembly {
  SMatrix total;
  std::vector<Element> elems;  // only filled when recording
  double yTop = 0;             // normalized
};

class StackSolver final : public Solver {
 public:
  explicit StackSolver(SolverOptions o) : opt_(std::move(o)) {}

  // ------------------------------------------------------------------ efficiencies
  Result solve(const Stack& stack, const Wavelength& wl) override {
    const Setup s = setup(stack, wl);
    return solveModels(stack, s, buildModels(stack, s));
  }

  /// solve() with prebuilt layer models (wavelength independent; Simulation, M12).
  Result solveWithModels(const Stack& stack, const Wavelength& wl, const detail::LayerModels& models) {
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
        const SMatrix Sp = assemble(s, polarizationSystem(s, p), models, false, {}, res.diag).total;
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
    const Assembly as = assemble(s, fullSystem(s), models, false, {}, res.diag);
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
    const Setup s = setup(stack, wl);
    const auto models = buildModels(stack, s);
    return samples(s, models, y_um, nullptr);
  }

  // ------------------------------------------------------------------ near field cuts
  NearFieldResult nearField(const Stack& stack, const Wavelength& wl, const NearFieldRequest& req) override {
    FG_PROF_SCOPE("nearfield");
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
    s.alpha.resize(nf);
    s.gamma.resize(nf);
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = s.map.orders(r);
      s.alpha(r) = s.inc.alpha0 + n * Kx;
      s.gamma(r) = s.inc.gamma0 + m * Kz;
    }
    s.ref = isotropicModes(s.epsSup, s.alpha, s.gamma);
    s.sub = &stack.material(stack.substrate);
    s.subModes = s.sub->isotropic ? isotropicModes(s.sub->eps, s.alpha, s.gamma)
                                  : anisotropicModes(s.sub->epsTensor, s.alpha, s.gamma);
    s.tr = resolveTransform(stack);
    if (stack.aperiodic.enabled) {
      const AperiodicSpec& ap = stack.aperiodic;
      FG_REQUIRE(!s.tr.x && !s.tr.z, "aperiodic: PML and coordinate_transform (ASR) cannot be combined");
      FG_REQUIRE(!ap.z.active() || stack.is2D(), "aperiodic: pml.z needs a 2D grating");
      checkPmlZones(stack);
      // total-field formulation (Plan "Einfall (a)"): the incident wave would be damped in the zones
      FG_REQUIRE(std::abs(s.inc.alpha0) < 1e-12 && std::abs(s.inc.gamma0) < 1e-12,
                 "aperiodic: plane-wave excitation with PML needs normal incidence (phi = 0)");
      s.aperiodic = true;
      if (ap.x.active()) s.pmlX = std::make_shared<const PmlMap>(stack.dx_um, ap.x.width_um, ap.x.beta);
      if (ap.z.active()) s.pmlZ = std::make_shared<const PmlMap>(stack.dz_um, ap.z.width_um, ap.z.beta);
    }
    if (s.tr.active() || s.aperiodic) transform(s);
    return s;
  }

  /// Transformed coordinates (Plan section "M17", eq:ct_deriv): A, G, the conversion matrices of
  /// eq:ct_rayleigh / eq:ct_inc and numerical bases of the half spaces (in X).
  static void transform(Setup& s) {
    FG_PROF_SCOPE("solve.transform");
    s.transformed = true;
    const int nf = s.map.size(), Nx = s.map.Nx(), Nz = s.map.Nz();
    std::vector<cplx> cx, cz;
    if (s.tr.x) cx = s.tr.x->inverseDerivativeCoefficients(2 * Nx);
    if (s.tr.z) cz = s.tr.z->inverseDerivativeCoefficients(2 * Nz);
    if (s.pmlX) cx = s.pmlX->inverseDerivativeCoefficients(2 * Nx);  // Plan eq:pml
    if (s.pmlZ) cz = s.pmlZ->inverseDerivativeCoefficients(2 * Nz);
    s.A = MatC::Zero(nf, nf);
    s.G = MatC::Zero(nf, nf);
    for (int r = 0; r < nf; ++r)
      for (int c = 0; c < nf; ++c) {
        const auto [n, m] = s.map.orders(r);
        const auto [n2, m2] = s.map.orders(c);
        // [[1/f']] is Toeplitz in n and the identity in m (and vice versa)
        if (m == m2) s.A(r, c) = (!cx.empty() ? cx[static_cast<std::size_t>(n - n2 + 2 * Nx)] : cplx(n == n2 ? 1.0 : 0.0)) * s.alpha(c);
        if (n == n2) s.G(r, c) = (!cz.empty() ? cz[static_cast<std::size_t>(m - m2 + 2 * Nz)] : cplx(m == m2 ? 1.0 : 0.0)) * s.gamma(c);
      }
    MatC Cx, Cz, Dx, Dz;
    if (s.tr.x) {
      Cx = s.tr.x->toPhysical(s.inc.alpha0 * s.k0, Nx);
      Dx = s.tr.x->fromPhysical(s.inc.alpha0 * s.k0, Nx);
    }
    if (s.tr.z) {
      Cz = s.tr.z->toPhysical(s.inc.gamma0 * s.k0, Nz);
      Dz = s.tr.z->fromPhysical(s.inc.gamma0 * s.k0, Nz);
    }
    s.toPhys = separable(s.map, s.tr.x ? &Cx : nullptr, s.tr.z ? &Cz : nullptr);
    s.fromPhys = separable(s.map, s.tr.x ? &Dx : nullptr, s.tr.z ? &Dz : nullptr);
    s.refPhys = std::move(s.ref);
    s.subPhys = std::move(s.subModes);
    s.ref = numericBasis(buildPropagationMatrix(homogeneousIsotropicQ(s.epsSup, nf), s.A, s.G));
    s.subModes = numericBasis(buildPropagationMatrix(
        s.sub->isotropic ? homogeneousIsotropicQ(s.sub->eps, nf) : homogeneousTensorQ(s.sub->epsTensor, nf), s.A, s.G));
  }

  /// True if the problem decouples exactly into TE and TM (Plan section "Spezialfaelle"): 1D,
  /// gamma = 0 for all orders, isotropic substrate and layers (the superstrate is always isotropic).
  bool decoupled(const Stack& stack, const Setup& s) const {
    if (!opt_.splitPolarizations || stack.is2D() || !s.gamma.isZero(0.0) || !s.sub->isotropic || s.transformed)
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
  /// Reduced systems (sys.pol set) only for the total S-matrix (record = false).
  /// The elements are independent (Plan section "M15", step 3): they are computed in blocks of
  /// sliceThreads() tasks in parallel and chained in stack order (sequential or fixed-shape tree).
  Assembly assemble(const Setup& s, const System& sys, const std::vector<std::unique_ptr<LayerModel>>& models,
                    bool record, const std::vector<double>& sampleY, Diagnostics& diag) {
    FG_PROF_SCOPE("solve.assemble");
    FG_REQUIRE(!record || !sys.pol, "field records need the full 4 n_f system");
    // resident chain (Plan section "M15", step 5): slice modes, slice S-matrices and star products
    // stay on the device; only the total S-matrix is downloaded (not for near-field records)
    const bool res = !record && linalg::gpuSelected(sys.ref.Psi.rows()) && linalg::gpu().resident();
    std::shared_ptr<const linalg::DeviceMatrix> psiInvD;
    if (res) {
      linalg::Backend& d = linalg::gpu();
      psiInvD = std::make_shared<const linalg::DeviceMatrix>(d.upload(sys.ref.PsiInv));
      d.synchronize();  // read by the task threads (their own streams)
    }
    std::vector<Task> tasks;
    double y = 0;
    for (const auto& mp : models) {
      const LayerModel& m = *mp;
      const double H = m.thickness_um() * s.k0;
      const double y0 = y;
      if (m.homogeneous()) {
        tasks.push_back([&s, &sys, &m, y0, H, record, res, psiInvD](std::vector<Element>& out, Diagnostics&) {
          const Material& mat = m.material();
          // analytic modes; in transformed coordinates numerical modes of M (Plan section "M17")
          auto modes = [&] {
            return s.transformed ? homogeneousModes(sys.M(m.Q(0.0), s))
                                 : layerModesHomogeneous(sys.basis(mat.isotropic
                                                                       ? isotropicModes(mat.eps, s.alpha, s.gamma)
                                                                       : anisotropicModes(mat.epsTensor, s.alpha, s.gamma)));
          };
          Element e;
          e.kind = Element::Kind::Modes;
          e.y0 = y0;
          e.y1 = y0 + H;
          e.model = &m;
          e.layerY0 = y0;
          if (res) {
            linalg::Backend& d = linalg::gpu();
            e.dS = resident::layerSMatrix(d, resident::upload(d, modes()), H, *psiInvD);
            d.synchronize();
            out.push_back(std::move(e));
            return;
          }
          RCWASliceRecord rec;
          e.S = stableLayerSMatrix(modes(), H, sys.ref, record ? &rec : nullptr);
          e.modes = std::move(rec.modes);
          e.G = std::move(rec.G);
          out.push_back(std::move(e));
        });
      } else if (opt_.kind == SolverKind::RCWA) {
        const int K = m.yInvariant() ? 1 : m.rcwaSlices;
        for (int k = 0; k < K; ++k)
          tasks.push_back([&s, &sys, &m, y0, H, K, k, record, res, psiInvD](std::vector<Element>& out, Diagnostics&) {
            const double yq = (k + 0.5) * m.thickness_um() / K;  // staircase: Q at slice centre
            const MatC M = sys.M(m.Q(yq), s);
            Element e;
            e.kind = Element::Kind::Modes;
            e.y0 = y0 + H * k / K;
            e.y1 = y0 + H * (k + 1) / K;
            e.model = &m;
            e.layerY0 = y0;
            e.qHeight_um = yq;
            if (res) {
              linalg::Backend& d = linalg::gpu();
              e.dS = resident::layerSMatrix(d, resident::layerModes(d, M, sys.pol), H / K, *psiInvD);
              d.synchronize();
              out.push_back(std::move(e));
              return;
            }
            RCWASliceRecord rec;
            e.S = stableLayerSMatrix(sys.modes(M), H / K, sys.ref, record ? &rec : nullptr);
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
    auto push = [&](Element& e) {
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
      e0.S = fromTransfer(linalg::product(sys.ref.PsiInv, sys.sub.Psi));  // substrate -> reference (Plan eq:Stot)
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
    if (res) {
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
  void rkLayer(const Setup& s, const System& sys, const LayerModel& m, double y0, double H, bool record,
               const std::vector<double>& sampleY, std::vector<Task>& tasks) {
    const double a = m.thickness_um();
    std::function<MatC(double)> Mfun;
    if (m.yInvariant()) {
      auto Mconst = std::make_shared<const MatC>(sys.M(m.Q(0.0), s));
      Mfun = [Mconst](double) { return *Mconst; };
    } else {
      Mfun = [&s, &sys, &m](double yl) { return sys.M(m.Q(yl / s.k0), s); };
    }
    double rho = rhoEstimate(s, m.epsList);
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
    for (int k = 0; k < K; ++k)
      tasks.push_back([this, &s, &sys, &m, &sampleY, Mfun, tab, y0, H, K, k, rho, record](std::vector<Element>& out,
                                                                                         Diagnostics& diag) {
        rkSlice(s, sys, m, Mfun, m.yInvariant(), y0, H * k / K, H * (k + 1) / K, rho, *tab, record, sampleY, out, diag,
                0);
      });
  }

  /// Field transfer over [ya, yb] (layer-local normalized y), recording the transfer matrices at
  /// the requested heights; bisects while the growth exceeds growthLimit.
  void rkSlice(const Setup& s, const System& sys, const LayerModel& m, const std::function<MatC(double)>& Mfun, bool constant,
               double y0, double ya, double yb, double rho, const ButcherTableau& tab, bool record,
               const std::vector<double>& sampleY, std::vector<Element>& out, Diagnostics& diag, int depth) const {
    auto transfer = [&](double from, double to) {
      RKStats st;
      MatC T = constant  ? rkTransferConstant(Mfun(0), to - from, tab, rkTol(), rho, st)
               : sys.pol ? rkTransfer(Mfun, from, to, tab, rkTol(), rho, st)
                         : rkTransferAuto(Mfun, from, to, tab, rkTol(), rho, st);
      diag.rkSteps += st.steps;
      return T;
    };
    Element e;
    e.kind = Element::Kind::RK;
    e.y0 = y0 + ya;
    e.y1 = y0 + yb;
    e.model = &m;
    e.layerY0 = y0;
    const Eigen::Index n = sys.ref.Psi.rows();
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
    e.S = rkSliceSMatrix(T, sys.ref);
    out.push_back(std::move(e));
  }

  // ------------------------------------------------------------------ efficiencies
  void efficiencies(const Setup& s, const SMatrix& S, Result& res) const {
    FG_PROF_SCOPE("solve.efficiencies");
    const int nf = s.map.size(), nh = 2 * nf;
    const VecC a_in = incidentAmplitudes(s);
    const int c = s.map.center();
    VecC b_top = S.S12 * a_in;
    VecC a_sub = S.S22 * a_in;
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
  /// W(r', r) = int over the physical region of exp(i k0 ((alpha_r - alpha_r') x + (gamma_r - gamma_r') z))
  /// (x, z in um): the flux of a field through y = const is 1/2 Re(Hx^H W Ez - Hz^H W Ex).
  static MatC physicalWeights(const Stack& stack, const Setup& s) {
    auto axis = [](double kappa, double a, double b) -> cplx {
      if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
      return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
    };
    const double x0 = s.pmlX ? s.pmlX->physicalBegin() : 0.0, x1 = s.pmlX ? s.pmlX->physicalEnd() : stack.dx_um;
    const double z0 = s.pmlZ ? s.pmlZ->physicalBegin() : 0.0, z1 = s.pmlZ ? s.pmlZ->physicalEnd() : stack.dz_um;
    const int nf = s.map.size();
    MatC W(nf, nf);
    for (int r = 0; r < nf; ++r)
      for (int c = 0; c < nf; ++c) {
        cplx w = axis(s.k0 * (s.alpha(c) - s.alpha(r)).real(), x0, x1);
        if (stack.is2D()) w *= axis(s.k0 * (s.gamma(c) - s.gamma(r)).real(), z0, z1);
        W(r, c) = w;
      }
    return W;
  }

  /// Outward power of the total field through the lateral sides of the box physical region x [0, top]
  /// (Poynting vector S_x, S_z from the field coefficients at Gauss nodes in y; per layer, and in a 1D
  /// profile layer also split at its height in the zones, where E_y, E_z jump on the sides).
  double lateralFlux(const Stack& stack, const Setup& s, const detail::LayerModels& models) {
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
    std::vector<double> ys, wy;
    for (std::size_t k = 0; k + 1 < br.size(); ++k) {
      const double a = br[k], b = br[k + 1];
      const int m = std::max(1, static_cast<int>(std::ceil((b - a) / hmax)));
      for (int p = 0; p < m; ++p)
        for (std::size_t q = 0; q < g.x.size(); ++q) {
          ys.push_back(a + (b - a) * (p + 0.5 * (1 + g.x[q])) / m);
          wy.push_back(0.5 * (b - a) / m * g.w[q]);
        }
    }
    if (ys.empty()) return 0.0;
    const FieldSampleSet fs = samples(s, models, ys, nullptr);
    const int nf = s.map.size();
    // side weights: phase of the orders at the side times the integral along the other lateral axis
    auto axisInt = [](double kappa, double a, double b) -> cplx {
      if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
      return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
    };
    const double x0 = s.pmlX ? s.pmlX->physicalBegin() : 0.0, x1 = s.pmlX ? s.pmlX->physicalEnd() : stack.dx_um;
    const double z0 = s.pmlZ ? s.pmlZ->physicalBegin() : 0.0, z1 = s.pmlZ ? s.pmlZ->physicalEnd() : stack.dz_um;
    auto sideX = [&](double xs) {  // W(r', r) for the plane x = xs
      MatC W(nf, nf);
      for (int r = 0; r < nf; ++r)
        for (int c = 0; c < nf; ++c) {
          cplx w = std::exp(I * (s.k0 * (s.alpha(c) - s.alpha(r)).real() * xs));
          if (stack.is2D()) w *= axisInt(s.k0 * (s.gamma(c) - s.gamma(r)).real(), z0, z1);
          W(r, c) = w;
        }
      return W;
    };
    auto sideZ = [&](double zs) {
      MatC W(nf, nf);
      for (int r = 0; r < nf; ++r)
        for (int c = 0; c < nf; ++c)
          W(r, c) = std::exp(I * (s.k0 * (s.gamma(c) - s.gamma(r)).real() * zs)) *
                    axisInt(s.k0 * (s.alpha(c) - s.alpha(r)).real(), x0, x1);
      return W;
    };
    struct Side {
      MatC W;
      int axis;     // 0: S_x, 1: S_z
      double sign;  // outward normal
    };
    std::vector<Side> sides;
    if (s.pmlX) {
      sides.push_back({sideX(x0), 0, -1.0});
      sides.push_back({sideX(x1), 0, +1.0});
    }
    if (s.pmlZ) {
      sides.push_back({sideZ(z0), 1, -1.0});
      sides.push_back({sideZ(z1), 1, +1.0});
    }
    double P = 0;
    for (std::size_t i = 0; i < ys.size(); ++i) {
      const FieldSample& f = fs.samples[i];
      const VecC Ex = f.F.segment(0, nf), Ez = f.F.segment(nf, nf), Hx = f.F.segment(2 * nf, nf),
                 Hz = f.F.segment(3 * nf, nf);
      for (const Side& sd : sides) {
        // S_x = 1/2 Re(Ey Hz* - Ez Hy*), S_z = 1/2 Re(Ex Hy* - Ey Hx*)
        const double sx = sd.axis == 0 ? 0.5 * std::real(Hz.dot(sd.W * f.Ey) - f.Hy.dot(sd.W * Ez))
                                       : 0.5 * std::real(f.Hy.dot(sd.W * Ex) - Hx.dot(sd.W * f.Ey));
        P += sd.sign * wy[i] * sx;
      }
    }
    return P;
  }

  /// Scattering outputs (Plan section "M18: Ausgaben"): the background stack (same transformation,
  /// same N, so the discretization errors cancel) gives E_bg; fluxes over the physical region. The
  /// absorption is the inward total-field flux through the closed box physical region x [0, top] (top,
  /// bottom and the lateral sides), minus that of the background.
  void scattering(const Stack& stack, const Setup& s, const detail::LayerModels& models, Result& res) {
    FG_PROF_SCOPE("solve.scattering");
    const Stack bgStack = backgroundStack(stack);
    const auto bgModels = buildModels(bgStack, s);
    const SMatrix S = assemble(s, fullSystem(s), models, false, {}, res.diag).total;
    const SMatrix Sb = assemble(s, fullSystem(s), bgModels, false, {}, res.diag).total;
    const int nf = s.map.size(), nh = 2 * nf;
    const VecC a_in = incidentAmplitudes(s);
    const VecC b = S.S12 * a_in, bb = Sb.S12 * a_in, t = S.S22 * a_in, tb = Sb.S22 * a_in;
    const MatC W = physicalWeights(stack, s);
    auto flux = [&](const VecC& F) {  // upward power through y = const over the physical region
      const VecC Ex = F.segment(0, nf), Ez = F.segment(nf, nf), Hx = F.segment(2 * nf, nf), Hz = F.segment(3 * nf, nf);
      return 0.5 * std::real(Hx.dot(W * Ez) - Hz.dot(W * Ex));  // dot: conjugates the first argument
    };
    const MatC& up = s.ref.Psi;
    const auto down = s.subModes.Psi.leftCols(nh);
    const cplx area = W(s.map.center(), s.map.center());
    const double Sinc = std::abs(flux(up.leftCols(nh) * a_in)) / area.real();  // incident intensity
    VecC top(2 * nh), topb(2 * nh);
    top << a_in, b;
    topb << a_in, bb;
    ScatteringResult& sc = res.scattering;
    sc.present = true;
    sc.up = flux(up.rightCols(nh) * (b - bb)) / Sinc;
    sc.down = -flux(down * (t - tb)) / Sinc;
    const double side = lateralFlux(stack, s, models);  // the background has none (laterally uniform)
    sc.absorbed = ((flux(down * t) - flux(up * top) - side) - (flux(down * tb) - flux(up * topb))) / Sinc;
    res.diag.energyBalance = 1.0;  // no diffraction orders; the M16 safeguard does not apply
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
    assembly_ = assemble(s, fullSystem(s), models, true, sorted, diag);
    const auto& el = assembly_.elems;
    const std::size_t K = el.size() - 1;
    const int nf = s.map.size(), nh = 2 * nf;
    // cumulative S-matrices below (S_L) and above (S_U) every interface q (top of element q)
    std::vector<SMatrix> SL(K + 1), SU(K + 1);
    SL[0] = el[0].S;
    for (std::size_t q = 1; q <= K; ++q) SL[q] = star(SL[q - 1], el[q].S);
    SU[K] = SMatrix::identity(nh);
    for (std::size_t q = K; q-- > 0;) SU[q] = star(el[q + 1].S, SU[q + 1]);
    const VecC a_in = incidentAmplitudes(s);
    const VecC b_in = VecC::Zero(nh);
    std::vector<VecC> amp(K + 1);  // (a_q; b_q) in the reference basis (Plan eq:internal)
    for (std::size_t q = 0; q <= K; ++q) amp[q] = interfaceAmplitudes(SL[q], SU[q], b_in, a_in);
    const VecC a_sub = SL[K].S22 * a_in;  // substrate amplitudes at y = 0 (total = SL[K])

    FieldSampleSet out;
    out.k0 = s.k0;
    out.alpha = s.alpha;
    out.gamma = s.gamma;
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
        Q = s.sub->isotropic ? homogeneousIsotropicQ(s.sub->eps, nf) : homogeneousTensorQ(s.sub->epsTensor, nf);
      } else if (y > assembly_.yTop) {
        fsm.region = 1;
        const double d = y - assembly_.yTop;
        VecC v(2 * nh);
        for (int c = 0; c < nh; ++c) {
          v(c) = a_in(c) * std::exp(s.ref.lambda(c) * d);
          v(nh + c) = amp[K](nh + c) * std::exp(s.ref.lambda(nh + c) * d);
        }
        fsm.F = s.ref.Psi * v;
        Q = homogeneousIsotropicQ(s.epsSup, nf);
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
          const VecC cc = e.G.partialPivLu().solve(rhs);
          VecC cu = cc.head(n), cd = cc.tail(n);
          for (Eigen::Index j = 0; j < n; ++j) {
            cu(j) *= std::exp(e.modes.lu(j) * (y - e.y0));
            cd(j) *= std::exp(e.modes.ld(j) * (y - e.y1));
          }
          fsm.F = e.modes.Wu * cu + e.modes.Wd * cd;
          Q = e.model->Q(e.qHeight_um >= 0 ? e.qHeight_um : 0.0);
        } else {
          const VecC F0 = s.ref.Psi * lo;
          fsm.F = F0;
          for (std::size_t j = 0; j < e.ys.size(); ++j)
            if (std::abs(e.ys[j] - y) <= 1e-14 * std::max(1.0, std::abs(y))) fsm.F = e.T[j] * F0;
          if (y >= e.y1 - 1e-14 * std::max(1.0, std::abs(y))) fsm.F = s.ref.Psi * hi;
          Q = e.model->Q((y - e.layerY0) / s.k0);
        }
      }
      fsm.Ey = s.transformed ? computeEy(Q, s.A, s.G, fsm.F) : computeEy(Q, s.alpha, s.gamma, fsm.F);
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
      const bool fff = req.reconstruction == Reconstruction::FFF && e && e->model && !e->model->homogeneous();
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
