// FourierGrating - discrete adjoint of the RCWA computation (Plan section "M19: Methode"): the forward
// chain is recorded (slice matrices, modes, S-matrices, prefix products) and differentiated in reverse:
// efficiencies -> amplitudes -> star products -> slice S-matrices -> eigen-decompositions -> M -> Q ->
// coefficient grids -> parameters (thickness, permittivity, breakpoints, shapes, profiles).
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
  enum class Kind { Thickness, Material, Breakpoint, Shape, ProfileDepth, ProfileShift };
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
  FG_REQUIRE(j.contains(jp) && j.at(jp).is_number(), "gradient: '" + ptr + "' is not a number of the configuration");
  const auto t = tokens(ptr);
  Param p;
  if (t.size() >= 3 && t[0] == "materials") {
    p.kind = Param::Kind::Material;
    p.material = t[1];
    FG_REQUIRE(t[2] == "n" || t[2] == "eps", "gradient: material parameters are n or eps ('" + ptr + "')");
    p.viaN = t[2] == "n";
    p.comp = t.size() >= 4 ? std::stoi(t[3]) : 0;
    FG_REQUIRE(p.material != cfg.stack.superstrate && p.material != cfg.stack.substrate,
               "gradient: the permittivity of superstrate/substrate is not available (M19 Phase 2): '" + ptr + "'");
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
  std::vector<Elem> el;
  {
    Elem e0;
    e0.S = fromTransfer(ref.PsiInv * sm.Psi);
    el.push_back(std::move(e0));
  }
  for (std::size_t k = 0; k < models.size(); ++k) {
    const auto& m = *models[k];
    const double H = m.thickness_um() * k0;
    if (m.homogeneous()) {
      const Material& mat = m.material();
      Elem e;
      e.layer = k;
      e.homogeneous = true;
      e.h = H;
      e.Q = m.Q(0.0);
      e.modes = layerModesHomogeneous(mat.isotropic ? isotropicModes(mat.eps, al, ga) : anisotropicModes(mat.epsTensor, al, ga));
      e.S = stableLayerSMatrix(e.modes, H, ref);
      el.push_back(std::move(e));
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
        e.Mc = Mc;
        e.y_um = 0.5 * (ya + yb) / k0;
        RKStats rs;
        if (e.constant) {
          e.T = rkTransferConstant(Mc, yb - ya, tab, opt.rkTol, rho, rs, &e.steps);
          for (double& y : e.steps.y) y += ya;
        } else {
          e.T = rkTransfer(Mfun, ya, yb, tab, opt.rkTol, rho, rs, &e.steps);
        }
        if ((!e.T.allFinite() || e.T.cwiseAbs().maxCoeff() > growth) && depth < 12) {
          slice(ya, 0.5 * (ya + yb), depth + 1);
          slice(0.5 * (ya + yb), yb, depth + 1);
          return;
        }
        e.S = fromTransfer(ref.PsiInv * e.T * ref.Psi);
        el.push_back(std::move(e));
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
      e.Q = m.Q(e.y_um);
      e.modes = layerModesAuto(buildPropagationMatrix(e.Q, al, ga));
      e.S = stableLayerSMatrix(e.modes, e.h, ref);
      el.push_back(std::move(e));
    }
  }
  std::vector<SMatrix> P(el.size());
  P[0] = el[0].S;
  for (std::size_t i = 1; i < el.size(); ++i) P[i] = star(P[i - 1], el[i].S);
  const VecC b = P.back().S12 * a_in, t = P.back().S22 * a_in;

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
  std::vector<SMatrix> Ebar(el.size());
  {
    SMatrix Pbar;
    Pbar.S11 = MatC::Zero(nh, nh);
    Pbar.S21 = MatC::Zero(nh, nh);
    Pbar.S12 = bbar * a_in.adjoint();
    Pbar.S22 = tbar * a_in.adjoint();
    for (std::size_t i = el.size() - 1; i >= 1; --i) {
      SMatrix Ab, Bb;
      adjoint::star(P[i - 1], el[i].S, Pbar, Ab, Bb);
      Ebar[i] = std::move(Bb);
      Pbar = std::move(Ab);
    }
  }
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
  for (std::size_t i = 1; i < el.size(); ++i) {
    const Elem& e = el[i];
    const auto& model = *models[e.layer];
    std::vector<std::string> mats;
    if (e.homogeneous) {
      mats = {std::get<HomogeneousLayer>(st.layers[e.layer].spec).material};
    } else {
      detail::AdjointSlice a;
      if (model.adjointSlice(e.y_um, a)) mats = a.materials;
    }
    if (!layerHasParams(e.layer, mats)) continue;
    // eigen-decompositions / stage matrices only if Q matters (materials, geometry, profile normal field)
    bool needQ = false;
    for (const auto& p : ps)
      needQ = needQ || (p.kind == Param::Kind::Material && std::find(mats.begin(), mats.end(), p.material) != mats.end()) ||
              (onLayer(p, e.layer) && (p.hasDnf || (p.kind != Param::Kind::Thickness && p.kind != Param::Kind::ProfileDepth)));
    if (e.rk) {
      // S = fromTransfer(Psi^-1 T Psi), T = P_m ... P_1 (fixed step sequence, Plan "RK-Propagator")
      const MatC Tref = ref.PsiInv * e.T * ref.Psi;
      MatC Trefb;
      adjoint::fromTransfer(Tref, Ebar[i], Trefb);
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
              QBlocks Qb;
              adjoint::propagationMatrix(model.Q(y_um), al, ga, nb[static_cast<std::size_t>(s)], Qb);
              fromQbar(e.layer, false, y_um, Qb);
            }
          }
        }
        G = Pk[k].adjoint() * G;
      }
      thicknessParam(e.layer, dLdH / model.thickness_um());
      if (needQ && e.constant) {
        QBlocks Qb;
        adjoint::propagationMatrix(model.Q(0.0), al, ga, Msum, Qb);
        fromQbar(e.layer, false, 0.0, Qb);
      }
      continue;
    }
    MatC Wub, Wdb;
    VecC lub, ldb;
    double hb = 0;
    adjoint::layerSMatrix(e.modes, e.h, ref, Ebar[i], Wub, Wdb, lub, ldb, hb);
    thicknessParam(e.layer, hb * k0 / e.K);
    if (!needQ) continue;
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
    QBlocks Qb;
    adjoint::propagationMatrix(e.Q, al, ga, Mb, Qb);
    fromQbar(e.layer, e.homogeneous, e.y_um, Qb);
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
