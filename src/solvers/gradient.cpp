// FourierGrating - discrete adjoint of the RCWA computation (Plan section "M19: Methode"): the forward
// chain is recorded (slice matrices, modes, S-matrices, prefix products) and differentiated in reverse:
// efficiencies -> amplitudes -> star products -> slice S-matrices -> eigen-decompositions -> M -> Q ->
// coefficient grids -> parameters (thickness, permittivity, breakpoints, shapes, profiles). Phase 2: wavelength,
// angles and the superstrate/substrate permittivity through alpha, gamma, k0, the half-space bases, the
// incident vector and the efficiency normalization.
#include "fg/adjoint/gradient.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <variant>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/adjoint/primitives.hpp"
#include "fg/adjoint/shape.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/io/config.hpp"
#include "fg/physics/efficiency.hpp"
#include "fg/physics/incidence.hpp"
#include "fg/physics/modes.hpp"
#include "fg/solvers/butcher.hpp"
#include "fg/solvers/rcwa.hpp"
#include "fg/solvers/rk.hpp"

#include "layer_model.hpp"

namespace fg {

namespace {

using json = nlohmann::json;
using Grid = fourier::CoefficientGrid;
constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

struct Param {
  enum class Kind { Thickness, Material, Breakpoint, Shape, ProfileDepth, ProfileShift, Lambda, Phi, Psi, Delta };
  Kind kind = Kind::Thickness;
  std::size_t layer = kNone;
  std::string material;
  bool viaN = false;
  int comp = 0;    ///< material: 0 = real part, 1 = imaginary part
  int index = 0;   ///< breakpoint index
  // shape: the structure at theta +- delta; the layers whose regions may move
  std::shared_ptr<const Stack> plus, minus;
  double delta = 0;
  std::vector<std::size_t> layers;
  /// profiles: d [[N_i N_j]] / d theta, d [[N_i]] / d theta (the analytic profile normal field is smooth and
  /// differentiated; the raster field of Algorithm NF of 2D cross-sections stays frozen)
  bool hasDnf = false;
  NormalFieldCoefficients dnf;
};

std::vector<std::string> tokens(const std::string& ptr) {
  std::vector<std::string> t;
  std::stringstream ss(ptr);
  std::string s;
  while (std::getline(ss, s, '/'))
    if (!s.empty()) t.push_back(s);
  return t;
}

/// Regions of a polygon2d / contour2d layer (normalized as in the forward computation).
std::vector<geometry::Region> regionsOf(const LayerSpec& spec) {
  if (const auto* c = std::get_if<Contour2DLayer>(&spec)) return c->regions;
  if (const auto* p = std::get_if<Polygon2DLayer>(&spec)) {
    std::vector<geometry::Region> out;
    for (const auto& inc : p->inclusions) {
      geometry::Region r;
      r.outer = geometry::Contour2D::polygon(inc.vertices);
      r.material = inc.material;
      r.normalize(1e-9);
      out.push_back(std::move(r));
    }
    return out;
  }
  throw Error("gradient: shape parameters need a polygon2d or contour2d layer");
}

bool sameControlPoints(const geometry::Region& a, const geometry::Region& b) {
  auto same = [](const geometry::Contour2D& x, const geometry::Contour2D& y) {
    if (x.segments.size() != y.segments.size()) return false;
    for (std::size_t s = 0; s < x.segments.size(); ++s)
      if (x.segments[s].P != y.segments[s].P || x.segments[s].w != y.segments[s].w) return false;
    return true;
  };
  if (!same(a.outer, b.outer) || a.holes.size() != b.holes.size()) return false;
  for (std::size_t h = 0; h < a.holes.size(); ++h)
    if (!same(a.holes[h], b.holes[h])) return false;
  return true;
}

Param resolve(const json& j, const io::RunConfig& cfg, const std::string& ptr, const IndexMap& map) {
  const json::json_pointer jp(ptr);
  const auto t = tokens(ptr);
  Param p;
  if (t.size() == 2 && t[0] == "incidence") {  // M19 Phase 2: at the current scan point (lists and ranges too)
    FG_REQUIRE(j.contains(jp), "gradient: '" + ptr + "' is not in the configuration");
    const std::string& f = t[1];
    if (f == "lambda_um") p.kind = Param::Kind::Lambda;
    else if (f == "phi_deg") p.kind = Param::Kind::Phi;
    else if (f == "psi_deg") p.kind = Param::Kind::Psi;
    else if (f == "delta_deg") p.kind = Param::Kind::Delta;
    else throw Error("gradient: unsupported incidence parameter '" + ptr + "'");
    return p;
  }
  FG_REQUIRE(j.contains(jp) && j.at(jp).is_number(), "gradient: '" + ptr + "' is not a number of the configuration");
  if (t.size() >= 3 && t[0] == "materials") {
    p.kind = Param::Kind::Material;
    p.material = t[1];
    FG_REQUIRE(t[2] == "n" || t[2] == "eps", "gradient: material parameters are n or eps ('" + ptr + "')");
    p.viaN = t[2] == "n";
    p.comp = t.size() >= 4 ? std::stoi(t[3]) : 0;
    FG_REQUIRE(p.material != cfg.stack.superstrate || p.comp == 0,
               "gradient: the superstrate stays lossless (no imaginary part as parameter): '" + ptr + "'");
    FG_REQUIRE(cfg.stack.material(p.material).isotropic, "gradient: anisotropic materials are not available: '" + ptr + "'");
    return p;
  }
  FG_REQUIRE(t.size() >= 3 && t[0] == "layers", "gradient: unsupported parameter '" + ptr + "'");
  for (const auto& l : j["layers"])
    FG_REQUIRE(l.value("type", std::string()) != "cad", "gradient: layer parameters with cad layers (M19 Phase 2)");
  p.layer = static_cast<std::size_t>(std::stoul(t[1]));
  FG_REQUIRE(p.layer < cfg.stack.layers.size(), "gradient: layer index out of range in '" + ptr + "'");
  const auto& spec = cfg.stack.layers[p.layer].spec;
  const std::string& f = t[2];
  if (f == "thickness_um") {
    p.kind = Param::Kind::Thickness;
  } else if (f == "depth_um" || f == "x0_um") {
    const auto* pr = std::get_if<Profile1DLayer>(&spec);
    FG_REQUIRE(pr && pr->shape == Profile1DLayer::Shape::Sinus, "gradient: depth_um/x0_um need a sinus profile1d layer");
    p.kind = f == "depth_um" ? Param::Kind::ProfileDepth : Param::Kind::ProfileShift;
    const double v = j.at(jp).get<double>(), delta = 1e-6 * std::max(1.0, std::abs(v));
    Profile1DLayer a = *pr, b = *pr;
    (f == "depth_um" ? a.depth_um : a.x0_um) += delta;
    (f == "depth_um" ? b.depth_um : b.x0_um) -= delta;
    const auto na = normalField(a, cfg.stack.dx_um, map), nb = normalField(b, cfg.stack.dx_um, map);
    p.dnf = na;
    auto diff = [&](Grid& d, const Grid& x, const Grid& y) {
      for (int pp = -x.Px(); pp <= x.Px(); ++pp)
        for (int qq = -x.Pz(); qq <= x.Pz(); ++qq) d(pp, qq) = (x(pp, qq) - y(pp, qq)) / (2 * delta);
    };
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) diff(p.dnf.NN[i][k], na.NN[i][k], nb.NN[i][k]);
      diff(p.dnf.N[i], na.N[i], nb.N[i]);
    }
    p.hasDnf = true;
  } else if (f == "breakpoints_um" && t.size() == 4) {
    const auto* l = std::get_if<Lamellar1DLayer>(&spec);
    FG_REQUIRE(l != nullptr, "gradient: breakpoints_um needs a lamellar1d layer");
    p.kind = Param::Kind::Breakpoint;
    p.index = std::stoi(t[3]);
    FG_REQUIRE(p.index >= 1 && p.index + 1 < static_cast<int>(l->breakpoints_um.size()),
               "gradient: only inner breakpoints are parameters (0 and the period are fixed)");
  } else if (f == "regions" || f == "inclusions") {
    p.kind = Param::Kind::Shape;
    (void)regionsOf(spec);  // polygon2d / contour2d
    const double v = j.at(jp).get<double>();
    p.delta = 1e-6 * std::max(1.0, std::abs(v));
    json jp1 = j, jm1 = j;
    jp1[jp] = v + p.delta;
    jm1[jp] = v - p.delta;
    p.plus = std::make_shared<const Stack>(io::parseConfig(jp1.dump()).stack);
    p.minus = std::make_shared<const Stack>(io::parseConfig(jm1.dump()).stack);
    p.layers = {p.layer};
  } else {
    throw Error("gradient: unsupported layer parameter '" + ptr + "'");
  }
  return p;
}

/// One recorded element of the chain.
struct Elem {
  std::size_t layer = kNone;
  bool homogeneous = false;
  double y_um = 0;  ///< slice centre (layer-relative)
  int K = 1;        ///< slices of the layer
  double h = 0;     ///< normalized slice thickness
  QBlocks Q;
  LayerModes modes;
  SMatrix S;
  // RK slices: transfer matrix, fixed step sequence (layer-local normalized y), constant M if y-invariant
  bool rk = false, constant = false;
  MatC T, Mc;
  RKSteps steps;
  double ya = 0, yb = 0, rho = 0;  ///< RK: slice range (normalized) and spectral bound, to recompute it
  bool stored = true;              ///< false: slice data dropped (checkpoint_every > 1), recomputed backwards
};

/// Spectral radius for the slice criterion (as the solver: eigenvalues up to n = 800, else power iteration).
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

/// Hermitian 2 x 2 form G of S_y(B) = B^H G B of an order (polarization identity on planeWaveOrder).
Eigen::Matrix2cd fluxForm(cplx eps, cplx al, cplx ga, cplx beta, int dir) {
  auto sy = [&](cplx bx, cplx bz) { return planeWaveOrder(eps, al, ga, beta, dir, bx, bz).Sy; };
  const double g11 = sy(1, 0), g22 = sy(0, 1);
  const double re = 0.5 * (sy(1, 1) - g11 - g22), im = -0.5 * (sy(1, I) - g11 - g22);
  Eigen::Matrix2cd G;
  G << g11, cplx(re, im), cplx(re, -im), g22;
  return G;
}

}  // namespace

GradientResult gradient(const std::string& configJson, const Objective& objective, const std::vector<std::string>& params) {
  const io::RunConfig cfg = io::parseConfig(configJson);
  FG_REQUIRE(!cfg.scan.empty(), "gradient: the configuration has no incidence");
  return gradient(configJson, cfg.scan.front(), objective, params);
}

GradientResult gradient(const std::string& configJson, const Wavelength& wl, const Objective& objective,
                        const std::vector<std::string>& params) {
  return detail::gradientWithModels(configJson, io::parseConfig(configJson), wl, nullptr, objective, params);
}

namespace {
GradientResult gradientCore(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                            const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps);
}  // namespace

GradientResult gradient(const SolverOptions& opt, const Stack& stack, const Wavelength& wl, const Objective& objective,
                        const std::vector<StackParameter>& params) {
  std::vector<Param> ps;
  for (const auto& sp : params) {
    Param p;
    p.kind = Param::Kind::Shape;
    p.delta = sp.delta;
    p.plus = std::make_shared<const Stack>(sp.at(sp.delta));
    p.minus = std::make_shared<const Stack>(sp.at(-sp.delta));
    FG_REQUIRE(p.plus->layers.size() == stack.layers.size() && p.minus->layers.size() == stack.layers.size(),
               "gradient: parameter '" + sp.name + "' changes the layer structure (critical heights, Plan M19 Phase 2)");
    for (std::size_t k = 0; k < stack.layers.size(); ++k) {
      const auto& s = stack.layers[k].spec;
      if (std::holds_alternative<Contour2DLayer>(s) || std::holds_alternative<Polygon2DLayer>(s) ||
          std::holds_alternative<SectionedLayer>(s))
        p.layers.push_back(k);
      else
        FG_REQUIRE(p.plus->layers[k].spec.index() == s.index(), "gradient: parameter '" + sp.name + "' changes a layer type");
    }
    ps.push_back(std::move(p));
  }
  return gradientCore(stack, opt, wl, nullptr, objective, std::move(ps));
}

namespace detail {

GradientResult gradientWithModels(const std::string& configJson, const io::RunConfig& cfg, const Wavelength& wl,
                                  const LayerModels* prebuilt, const Objective& objective,
                                  const std::vector<std::string>& params) {
  const json j = json::parse(configJson);
  std::vector<Param> ps;
  const IndexMap map(cfg.solver.Nx, cfg.stack.is2D() ? cfg.solver.Nz : 0);
  for (const auto& ptr : params) ps.push_back(resolve(j, cfg, ptr, map));
  return gradientCore(cfg.stack, cfg.solver, wl, prebuilt, objective, std::move(ps));
}

}  // namespace detail

namespace {

GradientResult gradientCore(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                            const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps) {
  FG_PROF_SCOPE("gradient");
  FG_REQUIRE(!st.aperiodic.enabled && !opt.scatteredField, "gradient: aperiodic structures are not available");
  FG_REQUIRE(!resolveTransform(st).active(), "gradient: coordinate transformations are not available");
  detail::applyThreadSettings(opt);
  const IndexMap map(opt.Nx, st.is2D() ? opt.Nz : 0);
  const int nf = map.size(), nh = 2 * nf, c = map.center();
  const double k0 = 2 * PI / wl.lambda_um;
  const Material& sup = st.material(st.superstrate);
  const Material& sub = st.material(st.substrate);
  FG_REQUIRE(sup.isotropic && std::abs(sup.eps.imag()) < 1e-14 && sub.isotropic,
             "gradient: isotropic substrate and lossless superstrate required");
  const IncidentWave inc = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
  const double Kx = 2 * PI / (k0 * st.dx_um), Kz = st.is2D() ? 2 * PI / (k0 * st.dz_um) : 0.0;
  VecC al(nf), ga(nf);
  for (int r = 0; r < nf; ++r) {
    const auto [n, m] = map.orders(r);
    al(r) = inc.alpha0 + n * Kx;
    ga(r) = inc.gamma0 + m * Kz;
  }
  const ModeBasis ref = isotropicModes(sup.eps, al, ga), sm = isotropicModes(sub.eps, al, ga);
  VecC a_in = VecC::Zero(nh);
  a_in(c) = inc.A(0);
  a_in(nf + c) = inc.A(2);

  // ---------------------------------------------------------------- forward (recorded)
  detail::LayerModels built;
  if (!prebuilt) built = detail::buildLayerModels(st, map, opt.rcwaSlices);
  const detail::LayerModels& models = prebuilt ? *prebuilt : built;
  const ButcherTableau tab = ButcherTableau::byName(opt.rkTableau);
  const std::size_t ckpt = static_cast<std::size_t>(std::max(1, opt.checkpointEvery));
  /// slice data of a recorded element from its description (forward pass and checkpoint recomputation)
  auto rebuild = [&](Elem& e) {
    const auto& m = *models[e.layer];
    if (e.homogeneous) {
      const Material& mat = m.material();
      e.Q = m.Q(0.0);
      e.modes = layerModesHomogeneous(mat.isotropic ? isotropicModes(mat.eps, al, ga) : anisotropicModes(mat.epsTensor, al, ga));
      e.S = stableLayerSMatrix(e.modes, e.h, ref);
    } else if (e.rk) {
      auto Mfun = [&m, &al, &ga, k0](double yl) { return buildPropagationMatrix(m.Q(yl / k0), al, ga); };
      RKStats rs;
      e.steps = RKSteps{};
      if (e.constant) {
        e.Mc = Mfun(0.0);
        e.T = rkTransferConstant(e.Mc, e.yb - e.ya, tab, opt.rkTol, e.rho, rs, &e.steps);
        for (double& y : e.steps.y) y += e.ya;
      } else {
        e.T = rkTransfer(Mfun, e.ya, e.yb, tab, opt.rkTol, e.rho, rs, &e.steps);
      }
      e.S = fromTransfer(ref.PsiInv * e.T * ref.Psi);
    } else {
      e.Q = m.Q(e.y_um);
      e.modes = layerModesAuto(buildPropagationMatrix(e.Q, al, ga));
      e.S = stableLayerSMatrix(e.modes, e.h, ref);
    }
    e.stored = true;
  };
  auto drop = [](Elem& e) {
    e.Q = QBlocks{};
    e.modes = LayerModes{};
    e.S = SMatrix{};
    e.T = MatC();
    e.Mc = MatC();
    e.steps = RKSteps{};
    e.stored = false;
  };
  std::vector<Elem> el;
  std::vector<SMatrix> P;  // prefix products; with checkpoints only every ckpt-th (and the last)
  SMatrix Plast;
  auto push = [&](Elem&& e) {
    const std::size_t i = el.size();
    Plast = i == 0 ? e.S : star(Plast, e.S);
    P.push_back(i % ckpt == 0 ? Plast : SMatrix{});
    el.push_back(std::move(e));
    if (ckpt > 1 && i > 0) drop(el.back());
  };
  {
    Elem e0;
    e0.S = fromTransfer(ref.PsiInv * sm.Psi);
    push(std::move(e0));
  }
  for (std::size_t k = 0; k < models.size(); ++k) {
    const auto& m = *models[k];
    const double H = m.thickness_um() * k0;
    if (m.homogeneous()) {
      const Material& mat = m.material();
      (void)mat;
      Elem e;
      e.layer = k;
      e.homogeneous = true;
      e.h = H;
      rebuild(e);
      push(std::move(e));
      continue;
    }
    if (opt.kind == SolverKind::RK) {
      // stability slices (Plan eq:slicecrit) and adaptive steps as in the solver; the steps are recorded
      auto Mfun = [&m, &al, &ga, k0](double yl) { return buildPropagationMatrix(m.Q(yl / k0), al, ga); };
      double rho = 0;
      for (int r = 0; r < nf; ++r)
        for (const cplx eps : m.epsList) rho = std::max(rho, std::abs(std::sqrt(al(r) * al(r) + ga(r) * ga(r) - eps)));
      MatC Mc;
      if (m.yInvariant()) {
        Mc = Mfun(0.0);
        rho = std::max(rho, spectralRadius(Mc));
      } else {
        for (double tt : {0.02, 0.25, 0.5, 0.75, 0.98}) rho = std::max(rho, spectralRadius(Mfun(tt * H)));
      }
      const double growth = opt.growthLimit > 0 ? opt.growthLimit : 1e8;
      const int K = std::max(1, static_cast<int>(std::ceil(rho * H / std::log(growth))));
      std::function<void(double, double, int)> slice = [&](double ya, double yb, int depth) {
        Elem e;
        e.layer = k;
        e.rk = true;
        e.constant = m.yInvariant();
        e.y_um = 0.5 * (ya + yb) / k0;
        e.ya = ya;
        e.yb = yb;
        e.rho = rho;
        rebuild(e);
        if ((!e.T.allFinite() || e.T.cwiseAbs().maxCoeff() > growth) && depth < 12) {
          slice(ya, 0.5 * (ya + yb), depth + 1);
          slice(0.5 * (ya + yb), yb, depth + 1);
          return;
        }
        push(std::move(e));
      };
      for (int s = 0; s < K; ++s) slice(H * s / K, H * (s + 1) / K, 0);
      continue;
    }
    const int K = m.yInvariant() ? 1 : m.rcwaSlices;
    for (int s = 0; s < K; ++s) {
      Elem e;
      e.layer = k;
      e.K = K;
      e.y_um = (s + 0.5) * m.thickness_um() / K;
      e.h = H / K;
      rebuild(e);
      push(std::move(e));
    }
  }
  const VecC b = Plast.S12 * a_in, t = Plast.S22 * a_in;

  // ---------------------------------------------------------------- efficiencies and objective
  GradientResult out;
  Result& res = out.result;
  const auto incOrder = planeWaveOrder(sup.eps, al(c), ga(c), ref.beta(c), -1, a_in(c), a_in(nf + c));
  const double Sinc = std::abs(incOrder.Sy);
  const bool subLossless = std::abs(sub.eps.imag()) < 1e-14 && sub.eps.real() > 0;
  VecC bbar = VecC::Zero(nh), tbar = VecC::Zero(nh);
  double sumR = 0, sumT = 0;
  auto weight = [&](char side, int n, int m) {
    double w = 0;
    for (const auto& [o, v] : objective.weights)
      if (o.side == side && o.n == n && o.m == m) w += v;
    return w;
  };
  for (int r = 0; r < nf; ++r) {
    const auto [n, m] = map.orders(r);
    const auto ro = planeWaveOrder(sup.eps, al(r), ga(r), ref.beta(r), +1, b(r), b(nf + r));
    if (ro.propagating) {
      const double e = ro.Sy / Sinc;
      res.reflected.push_back({n, m, e, ro.amp_s, ro.amp_p});
      sumR += e;
      const double w = weight('R', n, m);
      out.value += w * e;
      if (w != 0) {  // e = B^H G B / Sinc -> Bbar = 2 w G B / Sinc
        const Eigen::Matrix2cd G = fluxForm(sup.eps, al(r), ga(r), ref.beta(r), +1);
        const Eigen::Vector2cd B(b(r), b(nf + r)), Bb = 2 * w / Sinc * (G * B);
        bbar(r) += Bb(0);
        bbar(nf + r) += Bb(1);
      }
    }
    if (!subLossless) continue;
    const auto to = planeWaveOrder(sub.eps, al(r), ga(r), sm.beta(r), -1, t(r), t(nf + r));
    if (to.propagating) {
      const double e = -to.Sy / Sinc;
      res.transmitted.push_back({n, m, e, to.amp_s, to.amp_p});
      sumT += e;
      const double w = weight('T', n, m);
      out.value += w * e;
      if (w != 0) {
        const Eigen::Matrix2cd G = fluxForm(sub.eps, al(r), ga(r), sm.beta(r), -1);
        const Eigen::Vector2cd B(t(r), t(nf + r)), Bb = -2 * w / Sinc * (G * B);
        tbar(r) += Bb(0);
        tbar(nf + r) += Bb(1);
      }
    }
  }
  res.absorption = 1.0 - sumR - sumT;
  res.diag.energyBalance = sumR + sumT;

  // ---------------------------------------------------------------- reverse sweep
  out.grad.assign(ps.size(), 0.0);
  // M19 Phase 2: global parameters (wavelength, angles, superstrate/substrate) act on every element
  bool global = false;
  for (const auto& p : ps)
    global = global || p.kind == Param::Kind::Lambda || p.kind == Param::Kind::Phi || p.kind == Param::Kind::Psi ||
             p.kind == Param::Kind::Delta ||
             (p.kind == Param::Kind::Material && (p.material == st.superstrate || p.material == st.substrate));
  VecC alb = VecC::Zero(nf), gab = VecC::Zero(nf);
  MatC refPsib, refPsiInvb, smPsib;
  double dLdlnk0 = 0;  // all normalized thicknesses scale with k0
  // material cotangents (complex, Plan convention) per material name
  std::map<std::string, cplx> epsBar;
  auto onLayer = [](const Param& p, std::size_t k) {
    return p.layer == k || std::find(p.layers.begin(), p.layers.end(), k) != p.layers.end();
  };
  auto layerHasParams = [&](std::size_t k, const std::vector<std::string>& mats) {
    for (const auto& p : ps) {
      if (onLayer(p, k)) return true;
      if (p.kind == Param::Kind::Material && std::find(mats.begin(), mats.end(), p.material) != mats.end()) return true;
    }
    return false;
  };
  auto thicknessParam = [&](std::size_t k, double dLdH) {  // dL / d thickness (um)
    for (std::size_t q = 0; q < ps.size(); ++q)
      if (ps[q].layer == k && (ps[q].kind == Param::Kind::Thickness || ps[q].kind == Param::Kind::ProfileDepth))
        out.grad[q] += dLdH;
  };
  const int Px = 2 * map.Nx(), Pz = 2 * map.Nz();
  auto gridDot = [&](const Grid& gb, const Grid& d) {  // Re sum conj(gbar) d
    double s = 0;
    if (gb.Px() < Px || gb.Pz() < Pz) return 0.0;
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q) s += std::real(std::conj(gb(p, q)) * d(p, q));
    return s;
  };
  auto gridInner = [&](const Grid& gb, const Grid& chi) {  // sum gbar conj(chi)
    cplx s = 0;
    if (gb.Px() < Px || gb.Pz() < Pz) return s;
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q) s += gb(p, q) * std::conj(chi(p, q));
    return s;
  };
  auto epsOf = [&](const std::string& name) { return st.material(name).eps; };
  /// dL/dQ of the slice of layer k at height y_um -> parameters (materials, geometry, profiles)
  auto fromQbar = [&](std::size_t layer, bool homogeneous, double y_um, const QBlocks& Qb) {
    const auto& model = *models[layer];
    if (homogeneous) {
      const Material& mat = model.material();
      FG_REQUIRE(mat.isotropic, "gradient: anisotropic materials are not available");
      epsBar[std::get<HomogeneousLayer>(st.layers[layer].spec).material] += Qb(0, 0).trace() + Qb(1, 1).trace() + Qb(2, 2).trace();
      return;
    }
    detail::AdjointSlice a;
    FG_REQUIRE(model.adjointSlice(y_um, a),
               "gradient: layer " + std::to_string(layer) + " has no adjoint (anisotropic, FFF_B or transformed)");
    adjoint::IsoGridBar g;
    switch (a.fac) {
      case Factorization::Laurent: adjoint::qLaurent(a.sc, map, Qb, g); break;
      case Factorization::LiLamellar: adjoint::qLiLamellar(a.sc, map, Qb, g); break;
      case Factorization::FFF_iso: adjoint::qFFFIsotropic(a.sc, map, false, Qb, g); break;
      case Factorization::FFF_sym: adjoint::qFFFIsotropic(a.sc, map, true, Qb, g); break;
      case Factorization::FFF_rig: adjoint::qFFFIsotropicRigorous(a.sc, map, Qb, g); break;
      default: throw Error("gradient: factorization without adjoint");
    }
    // materials: eps grid = sum eps_k chi_k, 1/eps grid = sum chi_k / eps_k
    for (std::size_t k = 0; k < a.materials.size(); ++k) {
      const cplx ek = epsOf(a.materials[k]);
      epsBar[a.materials[k]] += gridInner(g.eps, a.chi[k]) + std::conj(-1.0 / (ek * ek)) * gridInner(g.invEps, a.chi[k]);
    }
    // a moving jump at x (um) between the materials left and right: d chi = (1/d) e^{-i p K x} (left grows)
    auto jump = [&](double x, const std::string& left, const std::string& right) {
      const cplx eL = epsOf(left), eR = epsOf(right);
      Grid de(Px, Pz), di(Px, Pz);
      for (int p = -Px; p <= Px; ++p) {
        const cplx d = std::exp(-I * (2 * PI * p * x / st.dx_um)) / st.dx_um;
        de(p, 0) = d * (eL - eR);
        di(p, 0) = d * (1.0 / eL - 1.0 / eR);
      }
      return gridDot(g.eps, de) + gridDot(g.invEps, di);
    };
    for (std::size_t q = 0; q < ps.size(); ++q) {
      const Param& p = ps[q];
      if (!onLayer(p, layer)) continue;
      if (p.kind == Param::Kind::Breakpoint) {
        const auto j0 = static_cast<std::size_t>(p.index);
        out.grad[q] += jump(a.breaks[j0], a.materials[j0 - 1], a.materials[j0]);
      }
      if (p.hasDnf) {  // the profile normal field moves with the parameter
        for (std::size_t ii = 0; ii < 3; ++ii) {
          for (std::size_t kk = 0; kk < 3; ++kk) out.grad[q] += gridDot(g.NN[ii][kk], p.dnf.NN[ii][kk]);
          out.grad[q] += gridDot(g.N[ii], p.dnf.N[ii]);
        }
      }
      if (p.kind == Param::Kind::ProfileShift) {
        // the crossings move with the profile: dx_j / dx0 = 1 (the period ends 0 and d are no crossings)
        for (std::size_t jj = 1; jj + 1 < a.breaks.size(); ++jj)
          out.grad[q] += jump(a.breaks[jj], a.materials[jj - 1], a.materials[jj]);
      } else if (p.kind == Param::Kind::Shape) {
        // regions of this layer at this height for theta and theta +- delta (Plan eq:shape2d / eq:shape3d)
        auto regionsAt = [&](const Stack& s, const detail::LayerModel* mdl) {
          if (std::holds_alternative<SectionedLayer>(s.layers[layer].spec)) {
            const auto& sl = std::get<SectionedLayer>(s.layers[layer].spec);
            return sl.provider->section(std::clamp(y_um, 1e-14 * sl.provider->thickness_um(),
                                                   sl.provider->thickness_um() * (1 - 1e-14))).regions;
          }
          (void)mdl;
          return regionsOf(s.layers[layer].spec);
        };
        const auto R0 = regionsAt(st, &model), Rp = regionsAt(*p.plus, nullptr), Rm = regionsAt(*p.minus, nullptr);
        FG_REQUIRE(R0.size() == Rp.size() && R0.size() == Rm.size() && R0.size() + 1 == a.materials.size(),
                   "gradient: the parameter changes the number of regions");
        const cplx eb = epsOf(a.materials[0]);
        for (std::size_t r = 0; r < R0.size(); ++r) {
          if (sameControlPoints(R0[r], Rp[r]) && sameControlPoints(R0[r], Rm[r])) continue;
          const std::size_t idx = r + 1;
          const Grid dchi = adjoint::shapeDerivative(R0[r], Rp[r], Rm[r], p.delta, st.dx_um, st.dz_um, Px, Pz);
          const cplx er = epsOf(a.materials[idx]);
          Grid de = dchi, di = dchi;
          for (int pp = -Px; pp <= Px; ++pp)
            for (int qq = -Pz; qq <= Pz; ++qq) {
              de(pp, qq) *= er - eb;
              di(pp, qq) *= 1.0 / er - 1.0 / eb;
            }
          out.grad[q] += gridDot(g.eps, de) + gridDot(g.invEps, di);
        }
      }
    }
  };
  /// backward of element i >= 1 from the cotangent Eb of its S-matrix
  auto backward = [&](std::size_t i, const SMatrix& Eb) {
    const Elem& e = el[i];
    const auto& model = *models[e.layer];
    std::vector<std::string> mats;
    if (e.homogeneous) {
      mats = {std::get<HomogeneousLayer>(st.layers[e.layer].spec).material};
    } else {
      detail::AdjointSlice a;
      if (model.adjointSlice(e.y_um, a)) mats = a.materials;
    }
    const bool own = layerHasParams(e.layer, mats);
    if (!own && !global) return;
    // eigen-decompositions / stage matrices only if Q matters (materials, geometry, profile normal field)
    bool needQ = false;
    for (const auto& p : ps)
      needQ = needQ || (p.kind == Param::Kind::Material && std::find(mats.begin(), mats.end(), p.material) != mats.end()) ||
              (onLayer(p, e.layer) && (p.hasDnf || (p.kind != Param::Kind::Thickness && p.kind != Param::Kind::ProfileDepth)));
    needQ = needQ || global;
    // dL/dM of a node or slice -> Q (own parameters) and alpha, gamma (global parameters)
    auto fromMbar = [&](const QBlocks& Q, const MatC& Mb, double y_um, bool homogeneous) {
      QBlocks Qb;
      if (global) adjoint::propagationMatrix(Q, al, ga, Mb, Qb, alb, gab);
      else adjoint::propagationMatrix(Q, al, ga, Mb, Qb);
      if (own) fromQbar(e.layer, homogeneous, y_um, Qb);
    };
    if (e.rk) {
      // S = fromTransfer(Psi^-1 T Psi), T = P_m ... P_1 (fixed step sequence, Plan "RK-Propagator")
      const MatC Tref = ref.PsiInv * e.T * ref.Psi;
      MatC Trefb;
      adjoint::fromTransfer(Tref, Eb, Trefb);
      if (global) {
        adjoint::accumulate(refPsiInvb, MatC(Trefb * (e.T * ref.Psi).adjoint()));
        adjoint::accumulate(refPsib, MatC((ref.PsiInv * e.T).adjoint() * Trefb));
      }
      MatC G = ref.PsiInv.adjoint() * Trefb * ref.Psi.adjoint();
      const std::size_t ms = e.steps.h.size();
      const Eigen::Index nn = e.T.rows();
      auto nodesAt = [&](std::size_t k) {
        std::vector<MatC> nodes;
        for (int s = 0; s < tab.stages(); ++s)
          nodes.push_back(e.constant ? e.Mc : MatC(buildPropagationMatrix(model.Q((e.steps.y[k] + tab.c(s) * e.steps.h[k]) / k0), al, ga)));
        return nodes;
      };
      std::vector<MatC> Pk(ms), Lk(ms);  // step propagators and prefix products L_k = P_k ... P_1
      for (std::size_t k = 0; k < ms; ++k) {
        Pk[k] = rkStepPropagator(nodesAt(k), e.steps.h[k], tab);
        Lk[k] = k == 0 ? Pk[k] : MatC(Pk[k] * Lk[k - 1]);
      }
      double dLdH = 0;
      MatC Msum = MatC::Zero(nn, nn);
      for (std::size_t k = ms; k-- > 0;) {
        const MatC Pb = k == 0 ? G : MatC(G * Lk[k - 1].adjoint());
        std::vector<MatC> nb;
        double hb = 0;
        const auto nodes = nodesAt(k);
        adjoint::rkStep(nodes, e.steps.h[k], tab, Pb, nb, &hb);
        dLdH += hb * e.steps.h[k];  // h_k proportional to the thickness
        if (needQ) {
          if (e.constant) {
            for (const auto& m : nb) Msum += m;
          } else {
            for (int s = 0; s < tab.stages(); ++s) {
              const double y_um = (e.steps.y[k] + tab.c(s) * e.steps.h[k]) / k0;
              fromMbar(model.Q(y_um), nb[static_cast<std::size_t>(s)], y_um, false);
            }
          }
        }
        G = Pk[k].adjoint() * G;
      }
      thicknessParam(e.layer, dLdH / model.thickness_um());
      dLdlnk0 += dLdH;
      if (needQ && e.constant) fromMbar(model.Q(0.0), Msum, 0.0, false);
      return;
    }
    MatC Wub, Wdb;
    VecC lub, ldb;
    double hb = 0;
    adjoint::layerSMatrix(e.modes, e.h, ref, Eb, Wub, Wdb, lub, ldb, hb, global ? &refPsiInvb : nullptr);
    thicknessParam(e.layer, hb * k0 / e.K);
    dLdlnk0 += hb * e.h;
    if (!needQ) return;
    const Eigen::Index n = e.modes.lu.size();
    MatC W(4 * nf, 2 * n), Wb(4 * nf, 2 * n);
    W << e.modes.Wu, e.modes.Wd;
    Wb << Wub, Wdb;
    VecC lam(2 * n), lamb(2 * n);
    lam << e.modes.lu, e.modes.ld;
    lamb << lub, ldb;
    std::vector<int> cls(static_cast<std::size_t>(2 * n), 0);
    for (Eigen::Index q = n; q < 2 * n; ++q) cls[static_cast<std::size_t>(q)] = 1;
    const MatC Mb = adjoint::eig(W, lam, Wb, lamb, cls);
    fromMbar(e.Q, Mb, e.y_um, e.homogeneous);
  };
  // reverse sweep over the star chain, per checkpoint segment (recompute the dropped slices and prefixes)
  SMatrix E0bar;
  {
    SMatrix Pbar;
    Pbar.S11 = MatC::Zero(nh, nh);
    Pbar.S21 = MatC::Zero(nh, nh);
    Pbar.S12 = bbar * a_in.adjoint();
    Pbar.S22 = tbar * a_in.adjoint();
    std::size_t i = el.size() - 1;
    while (i >= 1) {
      const std::size_t s0 = (i - 1) / ckpt * ckpt;  // P[s0] is stored
      for (std::size_t j = s0 + 1; j <= i; ++j) {
        if (!el[j].stored) rebuild(el[j]);
        if (j < i && P[j].S11.size() == 0) P[j] = star(P[j - 1], el[j].S);
      }
      for (std::size_t j = i; j > s0; --j) {
        SMatrix Ab, Bb;
        adjoint::star(P[j - 1], el[j].S, Pbar, Ab, Bb);
        Pbar = std::move(Ab);
        backward(j, Bb);
        if (ckpt > 1) {
          drop(el[j]);
          if (j % ckpt != 0) P[j] = SMatrix{};
        }
      }
      i = s0;
    }
    E0bar = std::move(Pbar);
  }
  if (global) {
    // first element S0 = fromTransfer(PsiInv_sup Psi_sub)
    MatC Xb;
    adjoint::fromTransfer(MatC(ref.PsiInv * sm.Psi), E0bar, Xb);
    adjoint::accumulate(refPsiInvb, MatC(Xb * sm.Psi.adjoint()));
    smPsib = ref.PsiInv.adjoint() * Xb;
    cplx epsSupb = 0, epsSubb = 0;
    adjoint::isotropicModes(sup.eps, al, ga, ref, refPsib, refPsiInvb, VecC(), epsSupb, alb, gab);
    adjoint::isotropicModes(sub.eps, al, ga, sm, smPsib, MatC(), VecC(), epsSubb, alb, gab);
    // amplitudes b = S12 a_in, t = S22 a_in
    const VecC ainb = Plast.S12.adjoint() * bbar + Plast.S22.adjoint() * tbar;
    // directional derivatives w.r.t. the scalars the bases, M and the amplitudes depend on
    double dAlpha0 = alb.real().sum(), dGamma0 = gab.real().sum(), dKx = 0, dKz = 0;
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = map.orders(r);
      dKx += n * alb(r).real();
      dKz += m * gab(r).real();
    }
    double dA0 = ainb(c).real(), dA2 = ainb(nf + c).real(), dEpsSup = epsSupb.real();
    cplx dEpsSub = epsSubb;
    // efficiency normalization and flux forms (b, t fixed): central differences of the objective
    struct Glob {
      double epsSup, alpha0, gamma0, Kx, Kz, A0, A2;
      cplx epsSub;
    };
    const Glob g0{sup.eps.real(), inc.alpha0, inc.gamma0, Kx, Kz, inc.A(0), inc.A(2), sub.eps};
    auto objectiveAt = [&](const Glob& g) {
      double v = 0;
      VecC a(nf), gg(nf);
      for (int r = 0; r < nf; ++r) {
        const auto [n, m] = map.orders(r);
        a(r) = g.alpha0 + n * g.Kx;
        gg(r) = g.gamma0 + m * g.Kz;
      }
      const auto io = planeWaveOrder(g.epsSup, a(c), gg(c), betaBranch(g.epsSup, a(c), gg(c)), -1, g.A0, g.A2);
      const double S = std::abs(io.Sy);
      for (int r = 0; r < nf; ++r) {
        const auto [n, m] = map.orders(r);
        if (const double w = weight('R', n, m); w != 0) {
          const auto ro = planeWaveOrder(g.epsSup, a(r), gg(r), betaBranch(g.epsSup, a(r), gg(r)), +1, b(r), b(nf + r));
          if (ro.propagating) v += w * ro.Sy / S;
        }
        if (!subLossless) continue;
        if (const double w = weight('T', n, m); w != 0) {
          const auto to = planeWaveOrder(g.epsSub, a(r), gg(r), betaBranch(g.epsSub, a(r), gg(r)), -1, t(r), t(nf + r));
          if (to.propagating) v -= w * to.Sy / S;
        }
      }
      return v;
    };
    auto central = [&](auto&& shift) {
      constexpr double h = 1e-6;
      Glob gp = g0, gm = g0;
      shift(gp, h);
      shift(gm, -h);
      return (objectiveAt(gp) - objectiveAt(gm)) / (2 * h);
    };
    dAlpha0 += central([](Glob& g, double h) { g.alpha0 += h; });
    dGamma0 += central([](Glob& g, double h) { g.gamma0 += h; });
    dKx += central([](Glob& g, double h) { g.Kx += h; });
    dKz += central([](Glob& g, double h) { g.Kz += h; });
    dA0 += central([](Glob& g, double h) { g.A0 += h; });
    dA2 += central([](Glob& g, double h) { g.A2 += h; });
    dEpsSup += central([](Glob& g, double h) { g.epsSup += h; });
    if (subLossless) {
      dEpsSub += central([](Glob& g, double h) { g.epsSub += h; });  // the substrate stays lossless
    } else {
      dEpsSub += cplx(central([](Glob& g, double h) { g.epsSub += h; }), central([](Glob& g, double h) { g.epsSub += cplx(0, h); }));
    }
    // incident wave (alpha0, gamma0, A) as a function of nu_sup and the angles
    auto incidenceDot = [&](double dnu, double dphi, double dpsi, double ddelta) {
      constexpr double h = 1e-6;
      const double nu = std::sqrt(sup.eps.real());
      const IncidentWave p = makeIncidence(nu + h * dnu, wl.phi + h * dphi, wl.psi + h * dpsi, wl.delta + h * ddelta);
      const IncidentWave m = makeIncidence(nu - h * dnu, wl.phi - h * dphi, wl.psi - h * dpsi, wl.delta - h * ddelta);
      return (dAlpha0 * (p.alpha0 - m.alpha0) + dGamma0 * (p.gamma0 - m.gamma0) + dA0 * (p.A(0) - m.A(0)) +
              dA2 * (p.A(2) - m.A(2))) / (2 * h);
    };
    epsBar[st.superstrate] += dEpsSup + incidenceDot(0.5 / std::sqrt(sup.eps.real()), 0, 0, 0);
    epsBar[st.substrate] += dEpsSub;
    constexpr double deg = PI / 180;
    for (std::size_t q = 0; q < ps.size(); ++q) {
      switch (ps[q].kind) {
        // Kx = lambda / dx, Kz = lambda / dz, normalized thicknesses ~ k0 = 2 pi / lambda
        case Param::Kind::Lambda: out.grad[q] = (dKx * Kx + dKz * Kz - dLdlnk0) / wl.lambda_um; break;
        case Param::Kind::Phi: out.grad[q] = incidenceDot(0, 1, 0, 0) * deg; break;
        case Param::Kind::Psi: out.grad[q] = incidenceDot(0, 0, 1, 0) * deg; break;
        case Param::Kind::Delta: out.grad[q] = incidenceDot(0, 0, 0, 1) * deg; break;
        default: break;
      }
    }
  }
  for (std::size_t q = 0; q < ps.size(); ++q) {
    const Param& p = ps[q];
    if (p.kind != Param::Kind::Material) continue;
    const cplx eb = epsBar[p.material];
    if (!p.viaN) {
      out.grad[q] = p.comp == 0 ? eb.real() : eb.imag();
    } else {
      const cplx nidx = std::sqrt(st.material(p.material).eps);  // eps = n^2
      const cplx deps = p.comp == 0 ? 2.0 * nidx : 2.0 * nidx * I;
      out.grad[q] = std::real(std::conj(eb) * deps);
    }
  }
  out.notes = {"RCWA slices fixed (Plan M19 \"Methode\" 1.)", "normal field frozen (Plan M19 \"Methode\" 2.)"};
  return out;
}

}  // namespace

}  // namespace fg
