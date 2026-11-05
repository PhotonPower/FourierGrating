// FourierGrating - stack assembly shared by the RK and RCWA solvers, efficiencies and near fields
// (Plan sections "S-Matrix-Formalismus", "Loeser A/B", "Nahfeld").
#include <algorithm>
#include <exception>
#include <cmath>
#include <functional>
#include <optional>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/fourier/sampling.hpp"
#include "fg/fourier/toeplitz.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/physics/anisotropic.hpp"
#include "fg/physics/efficiency.hpp"
#include "fg/physics/modes.hpp"
#include "fg/solvers/rcwa.hpp"
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
};

/// Field system of an assembly: the full 4 n_f problem, or one polarization of the decoupled 1D
/// in-plane problem (2 n_f, Plan eq:TE / eq:TM).
struct System {
  std::optional<Polarization> pol;
  ModeBasis ref, sub;

  MatC M(const QBlocks& Q, const Setup& s) const {
    return pol ? buildPropagationMatrix(Q, s.alpha, *pol) : buildPropagationMatrix(Q, s.alpha, s.gamma);
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

double rhoEstimate(const Setup& s, const std::vector<cplx>& eps) {
  double rho = 0;
  for (Eigen::Index r = 0; r < s.alpha.size(); ++r)
    for (const cplx e : eps)
      rho = std::max(rho, std::abs(std::sqrt(s.alpha(r) * s.alpha(r) + s.gamma(r) * s.gamma(r) - e)));
  return rho;
}

/// Spectral radius for the slice criterion (Plan eq:slicecrit). Exact eigenvalues up to
/// n = 800 (all 1D cases); above (2D) a deterministic power iteration, which is cheap and errs on
/// the safe side (transient growth); the a-posteriori growth check bisects slices if needed.
double spectralRadius(const MatC& M) {
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
    return solveModels(stack, s, detail::buildLayerModels(stack, s.map, opt_.rcwaSlices));
  }

  /// solve() with prebuilt layer models (wavelength independent; Simulation, M12).
  Result solveWithModels(const Stack& stack, const Wavelength& wl, const detail::LayerModels& models) {
    return solveModels(stack, setup(stack, wl), models);
  }

  Result solveModels(const Stack& stack, const Setup& s, const detail::LayerModels& models) {
    Result res;
    for (const auto& m : models) res.diag.warnings.insert(res.diag.warnings.end(), m->warnings.begin(), m->warnings.end());
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
      return res;
    }
    const Assembly as = assemble(s, fullSystem(s), models, false, {}, res.diag);
    efficiencies(s, as.total, res);
    return res;
  }

  // ------------------------------------------------------------------ field coefficients
  FieldSampleSet fieldCoefficients(const Stack& stack, const Wavelength& wl, const std::vector<double>& y_um) override {
    const Setup s = setup(stack, wl);
    const auto models = detail::buildLayerModels(stack, s.map, opt_.rcwaSlices);
    return samples(s, models, y_um, nullptr);
  }

  // ------------------------------------------------------------------ near field cuts
  NearFieldResult nearField(const Stack& stack, const Wavelength& wl, const NearFieldRequest& req) override {
    const Setup s = setup(stack, wl);
    const auto models = detail::buildLayerModels(stack, s.map, opt_.rcwaSlices);
    std::vector<double> ys;
    for (int k = 0; k < req.ny; ++k)
      ys.push_back(req.ny == 1 ? req.y_min_um : req.y_min_um + (req.y_max_um - req.y_min_um) * k / (req.ny - 1));
    ys.push_back(req.horizontal_y_um);
    std::vector<const Element*> where;
    const FieldSampleSet fs = samples(s, models, ys, &where);
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
  // ------------------------------------------------------------------ setup
  Setup setup(const Stack& stack, const Wavelength& wl) const {
    detail::applyThreadSettings(opt_);
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
    return s;
  }

  /// True if the problem decouples exactly into TE and TM (Plan section "Spezialfaelle"): 1D,
  /// gamma = 0 for all orders, isotropic substrate and layers (the superstrate is always isotropic).
  bool decoupled(const Stack& stack, const Setup& s) const {
    if (!opt_.splitPolarizations || stack.is2D() || !s.gamma.isZero(0.0) || !s.sub->isotropic) return false;
    for (const auto& le : stack.layers) {
      const auto names = layerMaterials(le.spec);
      if (!names) return false;
      for (const auto& n : *names)
        if (!stack.material(n).isotropic) return false;
    }
    return true;
  }

  // ------------------------------------------------------------------ assembly
  /// Builds the total S-matrix. With record = true every element is kept (S, modes, RK transfer
  /// matrices at the requested normalized heights 'sampleY', sorted).
  /// Reduced systems (sys.pol set) only for the total S-matrix (record = false).
  Assembly assemble(const Setup& s, const System& sys, const std::vector<std::unique_ptr<LayerModel>>& models,
                    bool record, const std::vector<double>& sampleY, Diagnostics& diag) {
    FG_REQUIRE(!record || !sys.pol, "field records need the full 4 n_f system");
    Assembly as;
    Element e0;
    e0.S = fromTransfer(sys.ref.PsiInv * sys.sub.Psi);  // substrate -> reference (Plan eq:Stot)
    as.total = e0.S;
    if (record) as.elems.push_back(std::move(e0));
    double y = 0;
    for (const auto& mp : models) {
      const LayerModel& m = *mp;
      const double H = m.thickness_um() * s.k0;
      const double y0 = y;
      if (m.homogeneous()) {
        const Material& mat = m.material();
        const ModeBasis mb = sys.basis(mat.isotropic ? isotropicModes(mat.eps, s.alpha, s.gamma)
                                                     : anisotropicModes(mat.epsTensor, s.alpha, s.gamma));
        Element e;
        e.kind = Element::Kind::Modes;
        e.y0 = y0;
        e.y1 = y0 + H;
        e.model = &m;
        e.layerY0 = y0;
        RCWASliceRecord rec;
        e.S = stableLayerSMatrix(layerModesHomogeneous(mb), H, sys.ref, &rec);
        e.modes = rec.modes;
        e.G = rec.G;
        addElement(as, std::move(e), record);
      } else if (opt_.kind == SolverKind::RCWA) {
        const int K = m.yInvariant() ? 1 : m.rcwaSlices;
        for (int k = 0; k < K; ++k) {
          const double yq = (k + 0.5) * m.thickness_um() / K;  // staircase: Q at slice centre
          const MatC M = sys.M(m.Q(yq), s);
          Element e;
          e.kind = Element::Kind::Modes;
          e.y0 = y0 + H * k / K;
          e.y1 = y0 + H * (k + 1) / K;
          e.model = &m;
          e.layerY0 = y0;
          e.qHeight_um = yq;
          RCWASliceRecord rec;
          e.S = stableLayerSMatrix(sys.modes(M), H / K, sys.ref, &rec);
          e.modes = rec.modes;
          e.G = rec.G;
          addElement(as, std::move(e), record);
        }
      } else {
        rkLayer(s, sys, m, y0, H, record, sampleY, as, diag);
      }
      y += H;
    }
    as.yTop = y;
    return as;
  }

  void addElement(Assembly& as, Element&& e, bool record) {
    as.total = star(as.total, e.S);
    if (record) as.elems.push_back(std::move(e));
  }

  /// RK treatment of one layer: stability sub-slices (Plan eq:slicecrit) with the spectral radius
  /// of the truncated M, a posteriori bisection, transfer matrices at recorded heights.
  void rkLayer(const Setup& s, const System& sys, const LayerModel& m, double y0, double H, bool record,
               const std::vector<double>& sampleY, Assembly& as, Diagnostics& diag) {
    const double a = m.thickness_um();
    std::function<MatC(double)> Mfun;
    MatC Mconst;
    if (m.yInvariant()) {
      Mconst = sys.M(m.Q(0.0), s);
      Mfun = [&Mconst](double) { return Mconst; };
    } else {
      Mfun = [&](double yl) { return sys.M(m.Q(yl / s.k0), s); };
    }
    double rho = rhoEstimate(s, m.epsList);
    if (m.yInvariant())
      rho = std::max(rho, spectralRadius(Mconst));
    else
      for (double t : {0.02, 0.25, 0.5, 0.75, 0.98}) rho = std::max(rho, spectralRadius(Mfun(t * a * s.k0)));
    const int K = std::max(1, static_cast<int>(std::ceil(rho * H / std::log(opt_.growthLimit))));
    const auto tab = ButcherTableau::byName(opt_.rkTableau);
    for (int k = 0; k < K; ++k)
      rkSlice(s, sys, m, Mfun, m.yInvariant(), y0, H * k / K, H * (k + 1) / K, rho, tab, record, sampleY, as, diag, 0);
  }

  /// Field transfer over [ya, yb] (layer-local normalized y), recording the transfer matrices at
  /// the requested heights; bisects while the growth exceeds growthLimit.
  void rkSlice(const Setup& s, const System& sys, const LayerModel& m, const std::function<MatC(double)>& Mfun, bool constant,
               double y0, double ya, double yb, double rho, const ButcherTableau& tab, bool record,
               const std::vector<double>& sampleY, Assembly& as, Diagnostics& diag, int depth) {
    auto transfer = [&](double from, double to) {
      RKStats st;
      MatC T = constant  ? rkTransferConstant(Mfun(0), to - from, tab, opt_.rkTol, rho, st)
               : sys.pol ? rkTransfer(Mfun, from, to, tab, opt_.rkTol, rho, st)
                         : rkTransferAuto(Mfun, from, to, tab, opt_.rkTol, rho, st);
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
    if ((!T.allFinite() || T.cwiseAbs().maxCoeff() > opt_.growthLimit) && depth < 12) {
      const double ym = 0.5 * (ya + yb);
      rkSlice(s, sys, m, Mfun, constant, y0, ya, ym, rho, tab, record, sampleY, as, diag, depth + 1);
      rkSlice(s, sys, m, Mfun, constant, y0, ym, yb, rho, tab, record, sampleY, as, diag, depth + 1);
      return;
    }
    FG_REQUIRE(T.allFinite(), "RK transfer matrix not finite");
    e.S = rkSliceSMatrix(T, sys.ref);
    addElement(as, std::move(e), record);
  }

  // ------------------------------------------------------------------ efficiencies
  void efficiencies(const Setup& s, const SMatrix& S, Result& res) const {
    const int nf = s.map.size(), nh = 2 * nf;
    const VecC a_in = incidentAmplitudes(s);
    const int c = s.map.center();
    const VecC b_top = S.S12 * a_in;
    const VecC a_sub = S.S22 * a_in;
    const auto incOrder = planeWaveOrder(s.epsSup, s.alpha(c), s.gamma(c), s.ref.beta(c), -1, a_in(c), a_in(nf + c));
    const double Sinc = std::abs(incOrder.Sy);
    const Material& sub = *s.sub;
    const bool subLossless = sub.isotropic ? (std::abs(sub.eps.imag()) < 1e-14 && sub.eps.real() > 0)
                                           : (sub.epsTensor - sub.epsTensor.adjoint()).cwiseAbs().maxCoeff() < 1e-14;
    const VecC Fsub = s.subModes.Psi.leftCols(nh) * a_sub;
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
      const auto ro = planeWaveOrder(s.epsSup, s.alpha(r), s.gamma(r), s.ref.beta(r), +1, b_top(r), b_top(nf + r));
      if (ro.propagating) {
        OrderEfficiency o{n, m, 0, ro.amp_s, ro.amp_p};
        split(o, ro.Sy / Sinc);
        res.reflected.push_back(o);
        sumR += ro.Sy / Sinc;
      }
      if (!subLossless) continue;
      if (sub.isotropic) {
        const auto to = planeWaveOrder(sub.eps, s.alpha(r), s.gamma(r), s.subModes.beta(r), -1, a_sub(r), a_sub(nf + r));
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

  static VecC incidentAmplitudes(const Setup& s) {
    const int nf = s.map.size();
    VecC a_in = VecC::Zero(2 * nf);
    a_in(s.map.center()) = s.inc.A(0);
    a_in(nf + s.map.center()) = s.inc.A(2);
    return a_in;
  }

  // ------------------------------------------------------------------ field coefficients
  FieldSampleSet samples(const Setup& s, const std::vector<std::unique_ptr<LayerModel>>& models,
                         const std::vector<double>& y_um, std::vector<const Element*>* where) {
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
      fsm.Ey = computeEy(Q, s.alpha, s.gamma, fsm.F);
      fsm.Hy = computeHy(s.alpha, s.gamma, fsm.F);
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
    // phase factors exp(i alpha_r x), exp(i gamma_r z)
    std::vector<VecC> px, pz;
    for (double x_um : xs) {
      VecC v(nf);
      for (int r = 0; r < nf; ++r) v(r) = std::exp(I * (s.alpha(r) * x_um * s.k0));
      px.push_back(v);
    }
    for (double z_um : zs) {
      VecC v(nf);
      for (int r = 0; r < nf; ++r) v(r) = std::exp(I * (s.gamma(r) * z_um * s.k0));
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
        const auto& nn = e->model->NN();
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
              const ProjectorSet ps = projectorMatrices(eps, e->model->normalAt(x_um, z_um));
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

void applyThreadSettings(const SolverOptions& opt) {
  linalg::setPreference(opt.backend, opt.backendThreshold);  // thread-local (M10)
  fourier::setFFTSampling(opt.fftSampling);                   // thread-local
}

Result solveWithModels(const SolverOptions& opt, const Stack& stack, const Wavelength& wl, const LayerModels& models) {
  return StackSolver(opt).solveWithModels(stack, wl, models);
}

}  // namespace detail

std::vector<Result> solveScan(const SolverOptions& opt, const Stack& stack, const std::vector<Wavelength>& scan) {
  std::vector<Result> out(scan.size());
  std::vector<std::exception_ptr> err(scan.size());
  const long n = static_cast<long>(scan.size());
#pragma omp parallel
  {
    auto solver = makeSolver(opt);
#pragma omp for schedule(dynamic, 1)
    for (long i = 0; i < n; ++i) {
      try {
        out[static_cast<std::size_t>(i)] = solver->solve(stack, scan[static_cast<std::size_t>(i)]);
      } catch (...) {
        err[static_cast<std::size_t>(i)] = std::current_exception();
      }
    }
  }
  for (const auto& e : err)
    if (e) std::rethrow_exception(e);
  return out;
}

}  // namespace fg
