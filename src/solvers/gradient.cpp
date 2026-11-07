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
#include <type_traits>
#include <variant>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/adjoint/primitives.hpp"
#include "fg/adjoint/shape.hpp"
#include "fg/fourier/contour.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/geometry/aperiodic.hpp"
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
  bool tensor = false;  ///< anisotropic material: d eps_tensor / d theta in dE (eps_tensor, eps_principal, euler_deg)
  Mat3C dE = Mat3C::Zero();
  int index = 0;   ///< breakpoint index
  // shape: the structure at theta + d (normal field difference quotient with the raster step, M19)
  std::function<Stack(double)> stackAt;
  // shape: the structure at theta +- delta; the layers whose regions may move
  std::shared_ptr<const Stack> plus, minus;
  double delta = 0;
  std::vector<std::size_t> layers;
  /// profiles: d [[N_i N_j]] / d theta, d [[N_i]] / d theta (the analytic profile normal field is smooth and
  /// differentiated; the raster field of Algorithm NF of 2D cross-sections: quotient over one raster cell, stackAt)
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
    if (!cfg.stack.material(p.material).isotropic) {  // M19 Phase 2: tensor entries, principal values, Euler angles
      FG_REQUIRE(t[2] == "eps_tensor" || t[2] == "eps_principal" || t[2] == "euler_deg",
                 "gradient: anisotropic material parameters are eps_tensor, eps_principal or euler_deg ('" + ptr + "')");
      FG_REQUIRE(p.material != cfg.stack.superstrate && p.material != cfg.stack.substrate,
                 "gradient: superstrate and substrate are isotropic");
      p.tensor = true;
      const double v = j.at(jp).get<double>(), h = 1e-6 * std::max(1.0, std::abs(v));
      json a = j, b = j;
      a[jp] = v + h;
      b[jp] = v - h;
      p.dE = (io::parseConfig(a.dump()).stack.material(p.material).epsTensor -
              io::parseConfig(b.dump()).stack.material(p.material).epsTensor) / (2 * h);
      return p;
    }
    FG_REQUIRE(t[2] == "n" || t[2] == "eps", "gradient: material parameters are n or eps ('" + ptr + "')");
    p.viaN = t[2] == "n";
    p.comp = t.size() >= 4 ? std::stoi(t[3]) : 0;
    FG_REQUIRE(p.material != cfg.stack.superstrate || p.comp == 0,
               "gradient: the superstrate stays lossless (no imaginary part as parameter): '" + ptr + "'");
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
  } else if (f == "depth_um" && std::holds_alternative<Profile2DLayer>(spec)) {
    // profile2d (egg_crate, sinus_x): the slices keep their relative height, the crossings stay; the analytic
    // normal field moves (isotropic: its coefficients, anisotropic: the coefficient grids of the slice)
    const auto& pr = std::get<Profile2DLayer>(spec);
    p.kind = Param::Kind::ProfileDepth;
    const double v = j.at(jp).get<double>(), delta = 1e-6 * std::max(1.0, std::abs(v));
    p.delta = delta;
    Profile2DLayer a = pr, b = pr;
    a.depth_um += delta;
    b.depth_um -= delta;
    const ResolvedTransform tr = resolveTransform(cfg.stack);
    const double dx = cfg.stack.dx_um, dz = cfg.stack.dz_um;
    const auto na = normalField(a, dx, dz, map, tr.x.get(), tr.z.get()), nb = normalField(b, dx, dz, map, tr.x.get(), tr.z.get());
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
  } else if ((f == "x0_um" || f == "z0_um") && std::holds_alternative<Profile2DLayer>(spec)) {
    // profile2d (egg_crate, sinus_x): translation of the whole slice, g_pq -> -i (p Kx | q Kz) g_pq
    const auto& pr = std::get<Profile2DLayer>(spec);
    FG_REQUIRE(pr.shape != Profile2DLayer::Shape::Function, "gradient: " + f + " needs an egg_crate or sinus_x layer");
    const ResolvedTransform tr = resolveTransform(cfg.stack);
    FG_REQUIRE(!(f == "x0_um" ? tr.x : tr.z), "gradient: " + f + " of profile2d with a coordinate transformation in " +
                                                  (f == "x0_um" ? "x" : "z"));
    p.kind = Param::Kind::ProfileShift;
    p.index = f == "x0_um" ? 0 : 1;
  } else if (f == "depth_um" || f == "x0_um") {
    const auto* pr = std::get_if<Profile1DLayer>(&spec);
    FG_REQUIRE(pr && pr->shape == Profile1DLayer::Shape::Sinus, "gradient: depth_um/x0_um need a sinus profile1d layer");
    p.kind = f == "depth_um" ? Param::Kind::ProfileDepth : Param::Kind::ProfileShift;
    const double v = j.at(jp).get<double>(), delta = 1e-6 * std::max(1.0, std::abs(v));
    p.delta = delta;
    Profile1DLayer a = *pr, b = *pr;
    (f == "depth_um" ? a.depth_um : a.x0_um) += delta;
    (f == "depth_um" ? b.depth_um : b.x0_um) -= delta;
    const ResolvedTransform tr = resolveTransform(cfg.stack);  // normals at x = f(X) (Plan section "M17")
    const auto na = normalField(a, cfg.stack.dx_um, map, tr.x.get()), nb = normalField(b, cfg.stack.dx_um, map, tr.x.get());
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
    p.stackAt = [j, jp, v](double d) {
      json a = j;
      a[jp] = v + d;
      return io::parseConfig(a.dump()).stack;
    };
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
  // scattered-field formulation (M18, Plan v2.13 Stufe B): slice [y0n, y1n] (normalized, absolute), background
  // layer, particular solution F_p = sum_k p_k e^{kappa_k (y - y_ref,k)} and the source amplitudes su, sv
  double y0n = 0, y1n = 0;
  std::size_t bgl = 0;
  std::vector<int> pk;      ///< background exponent index k of each term
  std::vector<VecC> pv;     ///< p_k = (kappa_k - M)^-1 sigma_k amp_k
  std::vector<cplx> e0, e1; ///< e^{kappa_k (y0n - y_ref)}, e^{kappa_k (y1n - y_ref)}
  VecC su, sv, F0, F1;
};

/// The structure and incidence at theta + h for a parameter acting on the background field (scattered-field
/// formulation: central differences of the N = 0 background solution, of Tsrc and of I_0; h in parameter units).
void perturbedSetup(const Stack& st, const Wavelength& wl, const Param& p, double h, Stack& s2, Wavelength& w2) {
  s2 = st;
  w2 = wl;
  constexpr double deg = PI / 180;
  switch (p.kind) {
    case Param::Kind::Material: {
      Material& m = s2.materials.at(p.material);
      FG_REQUIRE(m.isotropic && !p.tensor, "gradient: anisotropic background materials with the scattered-field formulation");
      cplx deps = p.comp == 0 ? cplx(1.0) : I;
      if (p.viaN) deps *= 2.0 * std::sqrt(m.eps);
      m.eps += h * deps;
      m.epsTensor = m.eps * Mat3C::Identity();
      break;
    }
    case Param::Kind::Thickness:
      std::visit(
          [&](auto& l) {
            using T = std::decay_t<decltype(l)>;
            if constexpr (std::is_same_v<T, HomogeneousLayer> || std::is_same_v<T, Lamellar1DLayer> ||
                          std::is_same_v<T, Polygon2DLayer> || std::is_same_v<T, Contour2DLayer>)
              l.thickness_um += h;
            else
              throw Error("gradient: thickness of this layer type with the scattered-field formulation");
          },
          s2.layers[p.layer].spec);
      break;
    case Param::Kind::Lambda: w2.lambda_um += h; break;
    case Param::Kind::Phi: w2.phi += h * deg; break;
    case Param::Kind::Psi: w2.psi += h * deg; break;
    case Param::Kind::Delta: w2.delta += h * deg; break;
    default: throw Error("gradient: internal: no background perturbation for this parameter");
  }
}

/// Background field of the scattered-field formulation (as the solver, Plan section "Streufeldformulierung"): the
/// background stack at the incidence order only (N = 0), per layer four exponential modes (normalized y).
struct BgLayerField {
  double y0 = 0, y1 = 0;
  MatC Psi;
  VecC lambda, amp;  ///< amplitudes: down referenced at y1, up at y0
  std::string material;
};
std::vector<BgLayerField> backgroundField(const Stack& st, const SolverOptions& opt, const Wavelength& wl, double alpha0,
                                          double gamma0) {
  Stack bg = backgroundStack(st);
  bg.aperiodic = {};
  bg.transform = {};
  SolverOptions o = opt;
  o.Nx = 0;
  o.Nz = 0;
  o.kind = SolverKind::RCWA;
  o.scatteredField = false;
  o.splitPolarizations = false;
  o.checkpointEvery = 1;
  const double k0 = 2 * PI / wl.lambda_um;
  std::vector<double> ys{0.0};
  for (const auto& le : bg.layers) ys.push_back(ys.back() + std::get<HomogeneousLayer>(le.spec).thickness_um);
  const FieldSampleSet fs = makeSolver(o)->fieldCoefficients(bg, wl, ys);
  detail::applyThreadSettings(opt);
  const VecC a0 = VecC::Constant(1, alpha0), g0 = VecC::Constant(1, gamma0);
  std::vector<BgLayerField> out;
  for (std::size_t k = 0; k < bg.layers.size(); ++k) {
    BgLayerField L;
    L.material = std::get<HomogeneousLayer>(bg.layers[k].spec).material;
    const Material& mat = bg.material(L.material);
    L.y0 = ys[k] * k0;
    L.y1 = ys[k + 1] * k0;
    const ModeBasis mb = mat.isotropic ? isotropicModes(mat.eps, a0, g0) : anisotropicModes(mat.epsTensor, a0, g0);
    L.Psi = mb.Psi;
    L.lambda = mb.lambda;
    const VecC c0 = mb.PsiInv * fs.samples[k].F, c1 = mb.PsiInv * fs.samples[k + 1].F;
    L.amp.resize(4);
    L.amp << c1.head(2), c0.tail(2);
    out.push_back(std::move(L));
  }
  return out;
}

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
/// Objectives on the amplitudes (aperiodic structures, Plan section "Gradienten aperiodischer Strukturen"):
/// the up-going amplitudes b above and the down-going t below the stack (in X) with the bases and the incident
/// amplitudes; returns L and the cotangents bbar, tbar and fills the outputs of res.
struct AmpContext {
  const ModeBasis& ref;
  const ModeBasis& sub;
  const VecC& aIn;
};
/// Cotangents of the objective through the bases and the incident amplitudes directly (global parameters).
struct AmpBar {
  MatC refPsi, subPsi;
  VecC aIn;
  VecC refLam, subLam;  ///< cotangents of the half-space eigenvalues
  double dlnk0 = 0;     ///< dL / d ln k0 at fixed normalized quantities (sampling heights in um)
};
using AmpObjective =
    std::function<double(const AmpContext&, const VecC& b, const VecC& t, VecC& bbar, VecC& tbar, AmpBar& ab, Result& res)>;
/// Incident amplitudes from the superstrate basis (modal excitation: a guided mode; fixed, not differentiated).
using IncidentHook = std::function<VecC(const ModeBasis& ref)>;
GradientResult gradientCore(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                            const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps,
                            const AmpObjective* ampObj = nullptr, const IncidentHook* incident = nullptr);
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

GradientResult aperiodicGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                                 const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps);
GradientResult modalGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                             const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps);
GradientResult scatteredGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                                 const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps,
                                 double wUp, double wDown);

GradientResult gradientCore(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                            const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps,
                            const AmpObjective* ampObj, const IncidentHook* incident) {
  if (st.aperiodic.enabled && !ampObj) return aperiodicGradient(st, opt, wl, prebuilt, objective, std::move(ps));
  FG_PROF_SCOPE("gradient");
  FG_REQUIRE(!opt.scatteredField || st.aperiodic.enabled,
             "gradient: the scattered-field formulation of periodic gratings is not available");
  FG_REQUIRE(ampObj || (objective.scattering.empty() && objective.modes.empty()),
             "gradient: scattering and mode objectives need an aperiodic structure");
  const ResolvedTransform tr = resolveTransform(st);
  // PML zones of aperiodic structures (Plan section "M18"): a complex coordinate stretching, f(X) = X where
  // structures lie
  std::shared_ptr<const PmlMap> pmlX, pmlZ;
  if (st.aperiodic.enabled) {
    if (st.aperiodic.x.active()) pmlX = std::make_shared<const PmlMap>(st.dx_um, st.aperiodic.x.width_um, st.aperiodic.x.beta);
    if (st.aperiodic.z.active()) pmlZ = std::make_shared<const PmlMap>(st.dz_um, st.aperiodic.z.width_um, st.aperiodic.z.beta);
  }
  const bool trans = tr.active() || pmlX || pmlZ;
  // structured superstrate/substrate (aperiodic, M18): numerical bases of their y-invariant layer models
  std::unique_ptr<detail::LayerModel> halfSup, halfSub;
  auto halfSpaceModel = [&](const LayerEntry& le) {
    Stack t = st;
    t.layers = {le};
    auto m = detail::buildLayerModel(t, 0, IndexMap(opt.Nx, st.is2D() ? opt.Nz : 0), opt.rcwaSlices);
    FG_REQUIRE(m->yInvariant(), "a structured superstrate/substrate must be y-invariant");
    return m;
  };
  if (st.superstrateStructure) halfSup = halfSpaceModel(*st.superstrateStructure);
  if (st.substrateStructure) halfSub = halfSpaceModel(*st.substrateStructure);
  const QBlocks Qsup = halfSup ? halfSup->Q(0.0) : QBlocks{}, Qsub = halfSub ? halfSub->Q(0.0) : QBlocks{};
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
  // scattered-field formulation (M18): oblique incidence on a PML axis expands the scattered field without
  // Bloch phase there; the sources are formed in the Bloch basis (alB, gaB) and re-expanded (Tsrc)
  const bool oblX = pmlX && std::abs(inc.alpha0) >= 1e-12, oblZ = pmlZ && std::abs(inc.gamma0) >= 1e-12;
  const bool scat = st.aperiodic.enabled && (opt.scatteredField || oblX || oblZ), bloch = oblX || oblZ;
  VecC al(nf), ga(nf), alB(nf), gaB(nf);
  for (int r = 0; r < nf; ++r) {
    const auto [n, m] = map.orders(r);
    al(r) = (oblX ? 0.0 : inc.alpha0) + n * Kx;
    ga(r) = (oblZ ? 0.0 : inc.gamma0) + m * Kz;
    alB(r) = inc.alpha0 + n * Kx;
    gaB(r) = inc.gamma0 + m * Kz;
  }
  ModeBasis ref = isotropicModes(sup.eps, al, ga), sm = isotropicModes(sub.eps, al, ga);
  VecC a_ph = VecC::Zero(nh);  // physical incident amplitudes
  a_ph(c) = inc.A(0);
  a_ph(nf + c) = inc.A(2);
  VecC a_in = a_ph;
  // coordinate transformation (Plan section "M17"): M with A, G, numerical half-space bases in X; the
  // efficiencies from the physical Rayleigh amplitudes (eq:ct_rayleigh), the incidence via eq:ct_inc
  detail::TransformedCoordinates tc;
  ModeBasis refPhys, subPhys;
  MatC Tb, Tt;  // physical amplitudes (up in the superstrate, down in the substrate) from those in X
  if (trans) {
    for (const auto& p : ps) {
      // the transformation must not depend on the parameter ("auto" breakpoints follow the geometry)
      auto sameMap = [&](const ResolvedTransform& o) {
        auto same = [](const std::shared_ptr<const CoordinateMap>& a, const std::shared_ptr<const CoordinateMap>& b) {
          return (!a && !b) || (a && b && a->breakpoints() == b->breakpoints());
        };
        return same(o.x, tr.x) && same(o.z, tr.z);
      };
      bool moves = false;
      if (p.kind == Param::Kind::Breakpoint) {
        Stack sp = st;
        std::get<Lamellar1DLayer>(sp.layers[p.layer].spec).breakpoints_um[static_cast<std::size_t>(p.index)] += 1e-6 * st.dx_um;
        moves = !sameMap(resolveTransform(sp));
      } else if (p.plus) {
        moves = !sameMap(resolveTransform(*p.plus)) || !sameMap(resolveTransform(*p.minus));
      }
      FG_REQUIRE(!moves, "gradient: a geometry parameter moves the breakpoints of coordinate_transform (\"auto\"); "
                         "give breakpoints_um explicitly");
    }
    tc = detail::transformedCoordinates(map, tr, pmlX.get(), pmlZ.get(), al, ga, inc.alpha0 * k0, inc.gamma0 * k0);
    refPhys = std::move(ref);
    subPhys = std::move(sm);
    ref = detail::numericBasis(buildPropagationMatrix(halfSup ? Qsup : homogeneousIsotropicQ(sup.eps, nf), tc.A, tc.G));
    sm = detail::numericBasis(buildPropagationMatrix(halfSub ? Qsub : homogeneousIsotropicQ(sub.eps, nf), tc.A, tc.G));
    a_in = (ref.PsiInv * detail::componentwise(tc.fromPhys, refPhys.Psi.leftCols(nh) * a_ph)).head(nh);
    for (int q = 0; q < nh; ++q)  // no modes growing upwards in the incident field (as in the solver)
      if (ref.lambda(q).real() > 1e-8) a_in(q) = 0;
    auto cw = [&](const MatC& X) {
      MatC Y(X.rows(), X.cols());
      for (int k = 0; k < 4; ++k) Y.middleRows(k * nf, nf) = tc.toPhys * X.middleRows(k * nf, nf);
      return Y;
    };
    Tb = (refPhys.PsiInv * cw(ref.Psi.rightCols(nh))).bottomRows(nh);
    Tt = (subPhys.PsiInv * cw(sm.Psi.leftCols(nh))).topRows(nh);
  }
  if (incident) a_in = (*incident)(ref);
  // scattered-field formulation: background field, background propagation matrices, the re-expansion Tsrc
  std::vector<BgLayerField> bgf;
  std::vector<MatC> bgM, bgMB;
  std::vector<std::size_t> bgOf(st.layers.size(), kNone);
  MatC Tsrc;
  /// T(m, n) = 1/d int_phys exp(i (kx0 + (n - m) 2 pi / d) x) dx on the oblique PML axes (as the solver)
  auto makeTsrc = [&](double kx0, double kz0) {
    auto restrictT = [&](int N, double d, double w, double kk0) {
      MatC T(2 * N + 1, 2 * N + 1);
      const double xa = w, xb = d - w;
      for (int m = -N; m <= N; ++m)
        for (int n2 = -N; n2 <= N; ++n2) {
          const double kap = kk0 + (n2 - m) * (2 * PI / d);
          T(m + N, n2 + N) = std::abs(kap) * (xb - xa) < 1e-12
                                 ? cplx((xb - xa) / d)
                                 : (std::exp(I * (kap * xb)) - std::exp(I * (kap * xa))) / (I * kap * d);
        }
      return T;
    };
    MatC Tx, Tz;
    if (oblX) Tx = restrictT(map.Nx(), st.dx_um, st.aperiodic.x.width_um, kx0);
    if (oblZ) Tz = restrictT(map.Nz(), st.dz_um, st.aperiodic.z.width_um, kz0);
    return detail::separable(map, oblX ? &Tx : nullptr, oblZ ? &Tz : nullptr);
  };
  if (scat) {
    FG_REQUIRE(opt.checkpointEvery <= 1, "gradient: checkpoint_every with the scattered-field formulation");
    a_in = VecC::Zero(nh);  // the sources drive the scattered field
    bgf = backgroundField(st, opt, wl, inc.alpha0, inc.gamma0);
    std::size_t n = 0;
    for (std::size_t k = 0; k < st.layers.size(); ++k) {
      FG_REQUIRE(!std::holds_alternative<Profile1DLayer>(st.layers[k].spec),
                 "gradient: profile layers with the scattered-field formulation");
      Stack one = st;
      one.layers = {st.layers[k]};
      const std::size_t c1 = backgroundStack(one).layers.size();
      if (c1) bgOf[k] = n;
      n += c1;
    }
    for (const auto& L : bgf) {
      const Material& mat = st.material(L.material);
      FG_REQUIRE(mat.isotropic, "gradient: anisotropic background materials with the scattered-field formulation");
      const QBlocks Q = mat.isotropic ? homogeneousIsotropicQ(mat.eps, nf) : homogeneousTensorQ(mat.epsTensor, nf);
      bgM.push_back(trans ? buildPropagationMatrix(Q, tc.A, tc.G) : buildPropagationMatrix(Q, al, ga));
      if (bloch) bgMB.push_back(buildPropagationMatrix(Q, alB, gaB));
    }
    if (bloch) Tsrc = makeTsrc(k0 * inc.alpha0, k0 * inc.gamma0);
  }
  /// sigma = (M - M_bg) F_bg for the single-order background vector w (oblique PML: in the Bloch basis, re-expanded)
  auto sourceOf = [&](const MatC& M, const MatC& MB, std::size_t l, const VecC& w) {
    VecC sig = VecC::Zero(4 * nf);
    for (int k = 0; k < 4; ++k)
      if (w(k) != 0.0)
        sig += w(k) * (bloch ? VecC(MB.col(k * nf + c) - bgMB[l].col(k * nf + c)) : VecC(M.col(k * nf + c) - bgM[l].col(k * nf + c)));
    return bloch ? detail::componentwise(Tsrc, sig) : sig;
  };
  /// background field (incidence order, 4 components) of background layer l at the normalized height y
  auto bgFieldAt = [&](std::size_t l, double y, VecC* Ek = nullptr) {
    const BgLayerField& L = bgf[l];
    VecC f = VecC::Zero(4);
    if (Ek) Ek->resize(4);
    for (int k = 0; k < 4; ++k) {
      const cplx E = std::exp(L.lambda(k) * (y - (k < 2 ? L.y1 : L.y0)));
      if (Ek) (*Ek)(k) = E;
      f += L.Psi.col(k) * (L.amp(k) * E);
    }
    return f;
  };
  const ModeBasis& rp = trans ? refPhys : ref;
  const ModeBasis& sp = trans ? subPhys : sm;
  auto Mof = [&](const QBlocks& Q) {
    return trans ? buildPropagationMatrix(Q, tc.A, tc.G) : buildPropagationMatrix(Q, al, ga);
  };

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
      e.modes = trans ? detail::homogeneousModes(Mof(e.Q))
                      : layerModesHomogeneous(mat.isotropic ? isotropicModes(mat.eps, al, ga) : anisotropicModes(mat.epsTensor, al, ga));
      e.S = stableLayerSMatrix(e.modes, e.h, ref);
    } else if (e.rk) {
      auto Mfun = [&m, &Mof, k0](double yl) { return Mof(m.Q(yl / k0)); };
      RKStats rs;
      e.steps = RKSteps{};
      if (scat) {  // affine transfer F(yb) = T F(ya) + f: the augmented linear system (F; 1) (as the solver)
        if (e.constant) e.Mc = Mfun(0.0);
        const Eigen::Index n = 4 * nf;
        auto aug = [&](double yl) {
          const MatC Mm = e.constant ? e.Mc : Mfun(yl);
          const MatC MB = bloch ? buildPropagationMatrix(m.Q(yl / k0), alB, gaB) : MatC();
          MatC A = MatC::Zero(n + 1, n + 1);
          A.topLeftCorner(n, n) = Mm;
          A.topRightCorner(n, 1) = sourceOf(Mm, MB, e.bgl, bgFieldAt(e.bgl, e.y0n + yl));
          return A;
        };
        const MatC Ta = rkTransfer(aug, e.ya, e.yb, tab, opt.rkTol, e.rho, rs, &e.steps);
        e.T = Ta.topLeftCorner(n, n);
        e.F0 = Ta.topRightCorner(n, 1);
        e.S = fromTransfer(ref.PsiInv * e.T * ref.Psi);
        const VecC g = ref.PsiInv * e.F0;
        e.su = g.tail(nh) - e.S.S12 * g.head(nh);
        e.sv = -(e.S.S22 * g.head(nh));
        e.stored = true;
        return;
      }
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
      const MatC M = Mof(e.Q);
      e.modes = layerModesAuto(M);
      e.S = stableLayerSMatrix(e.modes, e.h, ref);
      if (scat) {  // particular solution per background exponent (as the solver's modalSources)
        const BgLayerField& L = bgf[e.bgl];
        const MatC MB = bloch ? buildPropagationMatrix(e.Q, alB, gaB) : MatC();
        const double lmax = std::max(e.modes.lu.cwiseAbs().maxCoeff(), e.modes.ld.cwiseAbs().maxCoeff());
        const Eigen::Index N = M.rows();
        e.F0 = VecC::Zero(N);
        e.F1 = VecC::Zero(N);
        e.pk.clear();
        e.pv.clear();
        e.e0.clear();
        e.e1.clear();
        for (int k = 0; k < 4; ++k) {
          if (L.amp(k) == 0.0) continue;
          const cplx kap = L.lambda(k);
          double dist = 1e300;
          for (Eigen::Index j = 0; j < e.modes.lu.size(); ++j)
            dist = std::min({dist, std::abs(kap - e.modes.lu(j)), std::abs(kap - e.modes.ld(j))});
          FG_REQUIRE(dist >= 1e-8 * lmax, "gradient: scattered field: background exponent resonates with a slice mode");
          const VecC pk = (kap * MatC::Identity(N, N) - M).partialPivLu().solve(sourceOf(M, MB, e.bgl, L.Psi.col(k))) * L.amp(k);
          const double yr = k < 2 ? L.y1 : L.y0;
          e.pk.push_back(k);
          e.pv.push_back(pk);
          e.e0.push_back(std::exp(kap * (e.y0n - yr)));
          e.e1.push_back(std::exp(kap * (e.y1n - yr)));
          e.F0 += pk * e.e0.back();
          e.F1 += pk * e.e1.back();
        }
        const VecC g0 = ref.PsiInv * e.F0, g1 = ref.PsiInv * e.F1;
        e.su = g1.tail(nh) - e.S.S11 * g0.tail(nh) - e.S.S12 * g1.head(nh);
        e.sv = g0.head(nh) - e.S.S21 * g0.tail(nh) - e.S.S22 * g1.head(nh);
      }
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
  std::vector<VecC> Pu;    // scattered field: source vectors u of the prefix products
  SMatrix Plast;
  VecC PuLast, PvLast;
  auto push = [&](Elem&& e) {
    const std::size_t i = el.size();
    if (scat) {
      const VecC z = VecC::Zero(nh);
      const SourcedSMatrix E{e.S, e.su.size() ? e.su : z, e.sv.size() ? e.sv : z};
      const SourcedSMatrix C = i == 0 ? E : star(SourcedSMatrix{Plast, PuLast, PvLast}, E);
      Plast = C.S;
      PuLast = C.u;
      PvLast = C.v;
      Pu.push_back(PuLast);
    } else {
      Plast = i == 0 ? e.S : star(Plast, e.S);
    }
    P.push_back(i % ckpt == 0 ? Plast : SMatrix{});
    el.push_back(std::move(e));
    if (ckpt > 1 && i > 0) drop(el.back());
  };
  {
    Elem e0;
    e0.S = fromTransfer(ref.PsiInv * sm.Psi);
    push(std::move(e0));
  }
  double yLayer = 0;  // normalized height of the layer bottom
  std::vector<double> layerBottom(models.size(), 0.0);
  for (std::size_t k = 0; k < models.size(); ++k) {
    const auto& m = *models[k];
    const double H = m.thickness_um() * k0;
    const double yBottom = yLayer;
    layerBottom[k] = yBottom;
    yLayer += H;
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
      auto Mfun = [&m, &Mof, k0](double yl) { return Mof(m.Q(yl / k0)); };
      double rho = 0;
      // transformed coordinates: |A| up to |alpha| / min f' = |alpha| / (1 - eta) (as in the solver)
      const double sx = tr.x ? 1.0 / (1.0 - tr.x->eta()) : 1.0, sz = tr.z ? 1.0 / (1.0 - tr.z->eta()) : 1.0;
      for (int r = 0; r < nf; ++r)
        for (const cplx eps : m.epsList) {
          const cplx a = sx * al(r), g = sz * ga(r);
          rho = std::max(rho, std::abs(std::sqrt(a * a + g * g - eps)));
        }
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
        e.y0n = yBottom;  // RK: the layer bottom
        if (scat) e.bgl = bgOf[k];
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
      e.y0n = yBottom + H * s / K;
      e.y1n = yBottom + H * (s + 1) / K;
      if (scat) e.bgl = bgOf[k];
      rebuild(e);
      push(std::move(e));
    }
  }
  const VecC b = scat ? PuLast : VecC(Plast.S12 * a_in), t = scat ? PvLast : VecC(Plast.S22 * a_in);
  VecC bP = b, tP = t;  // physical Rayleigh amplitudes (Plan eq:ct_rayleigh)
  if (trans) {
    bP = (refPhys.PsiInv * detail::componentwise(tc.toPhys, ref.Psi.rightCols(nh) * b)).tail(nh);
    tP = (subPhys.PsiInv * detail::componentwise(tc.toPhys, sm.Psi.leftCols(nh) * t)).head(nh);
  }

  // ---------------------------------------------------------------- efficiencies and objective
  GradientResult out;
  Result& res = out.result;
  const auto incOrder = planeWaveOrder(sup.eps, al(c), ga(c), rp.beta(c), -1, a_ph(c), a_ph(nf + c));
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
  for (int r = 0; r < nf && !ampObj; ++r) {
    const auto [n, m] = map.orders(r);
    const auto ro = planeWaveOrder(sup.eps, al(r), ga(r), rp.beta(r), +1, bP(r), bP(nf + r));
    if (ro.propagating) {
      const double e = ro.Sy / Sinc;
      res.reflected.push_back({n, m, e, ro.amp_s, ro.amp_p});
      sumR += e;
      const double w = weight('R', n, m);
      out.value += w * e;
      if (w != 0) {  // e = B^H G B / Sinc -> Bbar = 2 w G B / Sinc
        const Eigen::Matrix2cd G = fluxForm(sup.eps, al(r), ga(r), rp.beta(r), +1);
        const Eigen::Vector2cd B(bP(r), bP(nf + r)), Bb = 2 * w / Sinc * (G * B);
        bbar(r) += Bb(0);
        bbar(nf + r) += Bb(1);
      }
    }
    if (!subLossless) continue;
    const auto to = planeWaveOrder(sub.eps, al(r), ga(r), sp.beta(r), -1, tP(r), tP(nf + r));
    if (to.propagating) {
      const double e = -to.Sy / Sinc;
      res.transmitted.push_back({n, m, e, to.amp_s, to.amp_p});
      sumT += e;
      const double w = weight('T', n, m);
      out.value += w * e;
      if (w != 0) {
        const Eigen::Matrix2cd G = fluxForm(sub.eps, al(r), ga(r), sp.beta(r), -1);
        const Eigen::Vector2cd B(tP(r), tP(nf + r)), Bb = -2 * w / Sinc * (G * B);
        tbar(r) += Bb(0);
        tbar(nf + r) += Bb(1);
      }
    }
  }
  res.absorption = 1.0 - sumR - sumT;
  res.diag.energyBalance = sumR + sumT;
  AmpBar ampBar;
  if (ampObj) {  // aperiodic: objective on the amplitudes in X
    out.value = (*ampObj)(AmpContext{ref, sm, a_in}, b, t, bbar, tbar, ampBar, res);
    res.absorption = 0;
    res.diag.energyBalance = 1.0;
  }
  // physical amplitudes (global parameters with a transformation); aperiodic: none
  const VecC bbarP = ampObj ? VecC(VecC::Zero(nh)) : bbar, tbarP = ampObj ? VecC(VecC::Zero(nh)) : tbar;
  if (trans && !ampObj) {  // back to the amplitudes in X
    bbar = Tb.adjoint() * bbar;
    tbar = Tt.adjoint() * tbar;
  }

  // ---------------------------------------------------------------- reverse sweep
  out.grad.assign(ps.size(), 0.0);
  // M19 Phase 2: global parameters (wavelength, angles, superstrate/substrate) act on every element
  bool global = false;
  for (const auto& p : ps)
    global = global || p.kind == Param::Kind::Lambda || p.kind == Param::Kind::Phi || p.kind == Param::Kind::Psi ||
             p.kind == Param::Kind::Delta ||
             (p.kind == Param::Kind::Material && (p.material == st.superstrate || p.material == st.substrate));
  global = global || halfSup || halfSub;
  VecC alb = VecC::Zero(nf), gab = VecC::Zero(nf);
  MatC Abar, Gbar;  // transformed coordinates: cotangents of A = [[1/f']] alpha, G = [[1/g']] gamma
  MatC refPsib, refPsiInvb, smPsib;
  double dLdlnk0 = 0;  // all normalized thicknesses scale with k0
  // material cotangents (complex, Plan convention) per material name
  std::map<std::string, cplx> epsBar;
  std::map<std::string, Mat3C> tensorBar;  // anisotropic layers: cotangents of the region tensors (M19 Phase 2)
  /// d eps / d theta of a material parameter as a tensor (isotropic: d eps 1)
  auto materialDirection = [&](const Param& p) -> Mat3C {
    if (p.tensor) return p.dE;
    cplx deps = p.comp == 0 ? cplx(1.0) : I;
    if (p.viaN) deps *= 2.0 * std::sqrt(st.material(p.material).eps);  // eps = n^2
    return deps * Mat3C::Identity();
  };
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
  /// d g / d shift for a translation along x (axis 0) or z (axis 1): g_pq -> -i (p Kx | q Kz) g_pq
  auto translated = [&](Grid g, int axis) {
    for (int pp = -g.Px(); pp <= g.Px(); ++pp)
      for (int qq = -g.Pz(); qq <= g.Pz(); ++qq)
        g(pp, qq) *= -I * (axis == 0 ? 2 * PI * pp / st.dx_um : 2 * PI * qq / st.dz_um);
    return g;
  };
  auto gridInner = [&](const Grid& gb, const Grid& chi) {  // sum gbar conj(chi)
    cplx s = 0;
    if (gb.Px() < Px || gb.Pz() < Pz) return s;
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q) s += gb(p, q) * std::conj(chi(p, q));
    return s;
  };
  auto epsOf = [&](const std::string& name) { return st.material(name).eps; };
  /// a jump at the physical x (um): its position X, dX/dx and the intervals left and right of it in breaks
  struct JumpAt {
    double X, dXdx;
    std::size_t left, right;
  };
  auto jumpAt = [&](double x, const std::vector<double>& breaks) {
    JumpAt j{x, 1.0, 0, 0};
    if (tr.x) {
      j.X = tr.x->inverse(x);
      j.X -= st.dx_um * std::floor(j.X / st.dx_um);
      j.dXdx = 1.0 / tr.x->df(j.X);
    }
    auto at = [&](double X) {
      X -= st.dx_um * std::floor(X / st.dx_um);
      std::size_t k = 0;
      while (k + 2 < breaks.size() && X >= breaks[k + 1]) ++k;
      return k;
    };
    const double eps = 1e-9 * st.dx_um;
    j.left = at(j.X - eps);
    j.right = at(j.X + eps);
    return j;
  };
  /// d chi / d theta of a moving region (transformed coordinates: the mapped contour, central differences)
  auto regionDerivative = [&](const geometry::Region& r0, const geometry::Region& rPlus, const geometry::Region& rMinus,
                              double delta, int px, int pz) {
    if (!tr.x && !tr.z) return adjoint::shapeDerivative(r0, rPlus, rMinus, delta, st.dx_um, st.dz_um, px, pz);
    const Grid gp = fourier::regionCoefficients(rPlus, st.dx_um, st.dz_um, px, pz, tr.x.get(), tr.z.get());
    const Grid gm = fourier::regionCoefficients(rMinus, st.dx_um, st.dz_um, px, pz, tr.x.get(), tr.z.get());
    Grid d(px, pz);
    for (int p = -px; p <= px; ++p)
      for (int q = -pz; q <= pz; ++q) d(p, q) = (gp(p, q) - gm(p, q)) / (2 * delta);
    return d;
  };
  /// Shape parameters with a raster normal field (Algorithm NF): the layer at theta +- H, H = raster cell (M19:
  /// the normal field moves with the geometry; NF is not smooth below its raster resolution)
  auto nfModels = [&](const Param& p, std::size_t layer, const detail::LayerModel& model, double& H) {
    H = model.normalFieldStep();
    std::pair<std::unique_ptr<detail::LayerModel>, std::unique_ptr<detail::LayerModel>> m;
    if (H <= 0 || !p.stackAt) return m;
    m.first = detail::buildLayerModel(p.stackAt(H), layer, map, opt.rcwaSlices);
    m.second = detail::buildLayerModel(p.stackAt(-H), layer, map, opt.rcwaSlices);
    return m;
  };
  /// dL/dQ of the slice of layer k at height y_um -> parameters (materials, geometry, profiles)
  auto fromQbar = [&](std::size_t layer, bool homogeneous, double y_um, const QBlocks& Qb,
                      const detail::LayerModel* other = nullptr) {
    const auto& model = other ? *other : *models[layer];
    if (homogeneous) {
      const Material& mat = model.material();
      const std::string& name = std::get<HomogeneousLayer>(st.layers[layer].spec).material;
      if (mat.isotropic) {
        epsBar[name] += Qb(0, 0).trace() + Qb(1, 1).trace() + Qb(2, 2).trace();
      } else {  // Q_ij = eps_ij 1
        Mat3C& T = tensorBar.try_emplace(name, Mat3C::Zero()).first->second;
        for (int i = 0; i < 3; ++i)
          for (int jj = 0; jj < 3; ++jj) T(i, jj) += Qb(i, jj).size() ? Qb(i, jj).trace() : cplx(0);
      }
      return;
    }
    detail::AdjointSlice a;
    FG_REQUIRE(model.adjointSlice(y_um, a), "gradient: layer " + std::to_string(layer) +
                                                " has no adjoint");
    if (a.aniso) {
      std::vector<Grid> gb;
      adjoint::qAnisotropic(a.fac, a.grids, map, Qb, gb);
      auto dotAll = [&](const std::vector<Grid>& d) {  // Re sum conj(gbar) d
        double v = 0;
        for (std::size_t m = 0; m < gb.size(); ++m) v += gridDot(gb[m], d[m]);
        return v;
      };
      // material parameters: directional derivative of the grids along d eps / d theta of their regions
      // (holomorphic in the tensors: central step along the complex direction; two grid evaluations each)
      for (std::size_t q = 0; q < ps.size(); ++q) {
        const Param& p = ps[q];
        if (p.kind != Param::Kind::Material || std::find(a.materials.begin(), a.materials.end(), p.material) == a.materials.end())
          continue;
        const Mat3C D = materialDirection(p);
        constexpr double h = 1e-6;
        auto tp = a.tensors, tm = a.tensors;
        for (std::size_t r = 0; r < tp.size(); ++r)
          if (a.materials[r] == p.material) {
            tp[r] += h * D;
            tm[r] -= h * D;
          }
        const auto gp = a.gridsFor(tp), gm = a.gridsFor(tm);
        std::vector<Grid> d = gp;
        for (std::size_t m = 0; m < d.size(); ++m)
          for (int pp = -d[m].Px(); pp <= d[m].Px(); ++pp)
            for (int qq = -d[m].Pz(); qq <= d[m].Pz(); ++qq) d[m](pp, qq) = (gp[m](pp, qq) - gm[m](pp, qq)) / (2 * h);
        out.grad[q] += dotAll(d);
      }
      for (std::size_t q = 0; q < ps.size(); ++q) {
        const Param& p = ps[q];
        if (!onLayer(p, layer)) continue;
        if (p.kind == Param::Kind::Breakpoint) {  // jump between regions j-1 and j moves
          const double xp = std::get<Lamellar1DLayer>(st.layers[layer].spec).breakpoints_um[static_cast<std::size_t>(p.index)];
          const JumpAt ja = jumpAt(xp, a.breaks);
          const double x = ja.X;
          std::vector<cplx> fl(gb.size()), fr(gb.size());
          a.evalAt(ja.left, x, fl.data());
          a.evalAt(ja.right, x, fr.data());
          std::vector<Grid> d(gb.size(), Grid(Px, Pz));
          for (std::size_t m = 0; m < gb.size(); ++m)
            for (int pp = -Px; pp <= Px; ++pp)
              d[m](pp, 0) = ja.dXdx * std::exp(-I * (2 * PI * pp * x / st.dx_um)) / st.dx_um * (fl[m] - fr[m]);
          out.grad[q] += dotAll(d);
        } else if (p.kind == Param::Kind::ProfileShift && a.gridsShift) {  // transformed: no translation symmetry
          const auto gp = a.gridsShift(p.delta), gm = a.gridsShift(-p.delta);
          std::vector<Grid> d = gp;
          for (std::size_t m = 0; m < d.size(); ++m)
            for (int pp = -d[m].Px(); pp <= d[m].Px(); ++pp)
              for (int qq = -d[m].Pz(); qq <= d[m].Pz(); ++qq) d[m](pp, qq) = (gp[m](pp, qq) - gm[m](pp, qq)) / (2 * p.delta);
          out.grad[q] += dotAll(d);
        } else if (p.kind == Param::Kind::ProfileShift) {  // translation: g_p -> -i p K g_p
          std::vector<Grid> d = a.grids;
          for (auto& g : d) g = translated(g, p.index);
          out.grad[q] += dotAll(d);
        } else if (p.kind == Param::Kind::ProfileDepth) {  // same crossings, the normal field moves
          const auto gp = a.gridsDepth(p.delta), gm = a.gridsDepth(-p.delta);
          std::vector<Grid> d = gp;
          for (std::size_t m = 0; m < d.size(); ++m)
            for (int pp = -d[m].Px(); pp <= d[m].Px(); ++pp)
              for (int qq = -d[m].Pz(); qq <= d[m].Pz(); ++qq) d[m](pp, qq) = (gp[m](pp, qq) - gm[m](pp, qq)) / (2 * p.delta);
          out.grad[q] += dotAll(d);
        } else if (p.kind == Param::Kind::Shape) {  // 2D: Laurent eps_ij = sum_k t_k,ij chi_k; FFF_B: gridsShape
          FG_REQUIRE(a.fac == Factorization::Laurent || a.gridsShape,
                     "gradient: shapes of anisotropic layers need the Laurent or FFF_B factorization");
          const auto R0 = regionsOf(st.layers[layer].spec), Rp = regionsOf(p.plus->layers[layer].spec),
                     Rm = regionsOf(p.minus->layers[layer].spec);
          for (std::size_t r = 0; r < R0.size(); ++r) {
            if (sameControlPoints(R0[r], Rp[r]) && sameControlPoints(R0[r], Rm[r])) continue;
            if (a.gridsShape) {  // indicator part (raster fixed); the raster part below
              out.grad[q] += dotAll(a.gridsShape(r + 1, regionDerivative(R0[r], Rp[r], Rm[r], p.delta, a.shapePx, a.shapePz)));
              continue;
            }
            const Grid dchi = regionDerivative(R0[r], Rp[r], Rm[r], p.delta, Px, Pz);
            const Mat3C dt = a.tensors[r + 1] - a.tensors[0];
            std::vector<Grid> d(9, Grid(Px, Pz));
            for (int e = 0; e < 9; ++e)
              for (int pp = -Px; pp <= Px; ++pp)
                for (int qq = -Pz; qq <= Pz; ++qq) d[static_cast<std::size_t>(e)](pp, qq) = dt(e / 3, e % 3) * dchi(pp, qq);
            out.grad[q] += dotAll(d);
          }
          if (a.gridsRasterOf) {  // the raster normal field moves with the geometry
            double H = 0;
            const auto nm = nfModels(p, layer, model, H);
            if (nm.first) {
              const auto gp = a.gridsRasterOf(*nm.first), gm = a.gridsRasterOf(*nm.second);
              std::vector<Grid> d = gp;
              for (std::size_t m = 0; m < d.size(); ++m)
                for (int pp = -d[m].Px(); pp <= d[m].Px(); ++pp)
                  for (int qq = -d[m].Pz(); qq <= d[m].Pz(); ++qq) d[m](pp, qq) = (gp[m](pp, qq) - gm[m](pp, qq)) / (2 * H);
              out.grad[q] += dotAll(d);
            }
          }
        }
      }
      return;
    }
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
        const double xp = std::get<Lamellar1DLayer>(st.layers[layer].spec).breakpoints_um[static_cast<std::size_t>(p.index)];
        const JumpAt ja = jumpAt(xp, a.breaks);
        out.grad[q] += ja.dXdx * jump(ja.X, a.materials[ja.left], a.materials[ja.right]);
      }
      if (p.hasDnf) {  // the profile normal field moves with the parameter
        for (std::size_t ii = 0; ii < 3; ++ii) {
          for (std::size_t kk = 0; kk < 3; ++kk) out.grad[q] += gridDot(g.NN[ii][kk], p.dnf.NN[ii][kk]);
          out.grad[q] += gridDot(g.N[ii], p.dnf.N[ii]);
        }
      }
      if (p.kind == Param::Kind::ProfileShift && std::holds_alternative<Profile2DLayer>(st.layers[layer].spec)) {
        // profile2d: translation of all grids of the slice (the analytic normal field moves along)
        auto dot = [&](const Grid& gb, const Grid& cg) {
          return cg.Px() < Px || cg.Pz() < Pz ? 0.0 : gridDot(gb, translated(cg, p.index));
        };
        out.grad[q] += dot(g.eps, a.sc.eps) + dot(g.invEps, a.sc.invEps);
        if (a.sc.normal)
          for (std::size_t ii = 0; ii < 3; ++ii) {
            for (std::size_t kk = 0; kk < 3; ++kk) out.grad[q] += dot(g.NN[ii][kk], a.sc.normal->NN[ii][kk]);
            out.grad[q] += dot(g.N[ii], a.sc.normal->N[ii]);
          }
      } else if (p.kind == Param::Kind::ProfileShift) {
        // the crossings move with the profile: dx_j / dx0 = 1, dX_j / dx0 = 1 / f'(X_j) (the period ends 0
        // and d are no crossings; the extra breaks of the map have the same material on both sides)
        for (std::size_t jj = 1; jj + 1 < a.breaks.size(); ++jj)
          out.grad[q] += (tr.x ? 1.0 / tr.x->df(a.breaks[jj]) : 1.0) * jump(a.breaks[jj], a.materials[jj - 1], a.materials[jj]);
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
          const Grid dchi = regionDerivative(R0[r], Rp[r], Rm[r], p.delta, Px, Pz);
          const cplx er = epsOf(a.materials[idx]);
          Grid de = dchi, di = dchi;
          for (int pp = -Px; pp <= Px; ++pp)
            for (int qq = -Pz; qq <= Pz; ++qq) {
              de(pp, qq) *= er - eb;
              di(pp, qq) *= 1.0 / er - 1.0 / eb;
            }
          out.grad[q] += gridDot(g.eps, de) + gridDot(g.invEps, di);
        }
        double H = 0;
        const auto nm = nfModels(p, layer, model, H);
        if (nm.first) {
          const NormalFieldCoefficients &np = nm.first->NN(y_um), &nmi = nm.second->NN(y_um);
          auto ddot = [&](const Grid& gb, const Grid& x, const Grid& y) {
            if (x.Px() < Px || y.Px() < Px || x.Pz() < Pz || y.Pz() < Pz) return 0.0;
            Grid d(Px, Pz);
            for (int pp = -Px; pp <= Px; ++pp)
              for (int qq = -Pz; qq <= Pz; ++qq) d(pp, qq) = (x(pp, qq) - y(pp, qq)) / (2 * H);
            return gridDot(gb, d);
          };
          for (std::size_t ii = 0; ii < 3; ++ii) {
            for (std::size_t kk = 0; kk < 3; ++kk) out.grad[q] += ddot(g.NN[ii][kk], np.NN[ii][kk], nmi.NN[ii][kk]);
            out.grad[q] += ddot(g.N[ii], np.N[ii], nmi.N[ii]);
          }
        }
      }
    }
  };
  // cotangents of the background data per background layer (amplitudes, exponents, mode vectors, boundaries), of the
  // background propagation matrices, of Tsrc and of the Bloch-basis wave numbers
  const std::size_t nbg = bgf.size();
  std::vector<VecC> bgAmpB(nbg, VecC::Zero(4)), bgKapB(nbg, VecC::Zero(4));
  std::vector<MatC> bgPsiB(nbg, MatC::Zero(4, 4)), bgMb(nbg), bgMBb(nbg);
  std::vector<double> bgY0B(nbg, 0.0), bgY1B(nbg, 0.0);
  MatC TsrcB = bloch ? MatC(MatC::Zero(nf, nf)) : MatC();
  VecC alBb = VecC::Zero(nf), gaBb = VecC::Zero(nf);
  /// backward of element i >= 1 from the cotangent Eb of its S-matrix
  auto backward = [&](std::size_t i, const SMatrix& Eb, const MatC* Msrc = nullptr, const MatC* MBsrc = nullptr,
                      const VecC* fbar = nullptr) {
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
    if (!own && !global && !fbar) return;
    // eigen-decompositions / stage matrices only if Q matters (materials, geometry, profile normal field)
    bool needQ = false;
    for (const auto& p : ps)
      needQ = needQ || (p.kind == Param::Kind::Material && std::find(mats.begin(), mats.end(), p.material) != mats.end()) ||
              (onLayer(p, e.layer) && (p.hasDnf || (p.kind != Param::Kind::Thickness && p.kind != Param::Kind::ProfileDepth)));
    needQ = needQ || global || fbar;
    // dL/dM of a node or slice -> Q (own parameters) and alpha, gamma (global parameters)
    auto fromMbar = [&](const QBlocks& Q, const MatC& Mb, double y_um, bool homogeneous) {
      QBlocks Qb;
      if (trans && global) adjoint::propagationMatrix(Q, tc.A, tc.G, Mb, Qb, Abar, Gbar);
      else if (trans) adjoint::propagationMatrix(Q, tc.A, tc.G, Mb, Qb);
      else if (global) adjoint::propagationMatrix(Q, al, ga, Mb, Qb, alb, gab);
      else adjoint::propagationMatrix(Q, al, ga, Mb, Qb);
      if (own) fromQbar(e.layer, homogeneous, y_um, Qb);
    };
    if (e.rk && fbar) {
      // scattered field: the augmented system (F; 1), T_aug = [[T, f], [0, 1]]; the nodes [[M, sigma], [0, 0]]
      const MatC Tref = ref.PsiInv * e.T * ref.Psi;
      MatC Trefb;
      adjoint::fromTransfer(Tref, Eb, Trefb);
      if (global) {
        adjoint::accumulate(refPsiInvb, MatC(Trefb * (e.T * ref.Psi).adjoint()));
        adjoint::accumulate(refPsib, MatC((ref.PsiInv * e.T).adjoint() * Trefb));
      }
      const Eigen::Index n = 4 * nf;
      auto augAt = [&](double yl) {
        const MatC Mm = e.constant ? e.Mc : MatC(Mof(model.Q(yl / k0)));
        const MatC MB = bloch ? buildPropagationMatrix(model.Q(yl / k0), alB, gaB) : MatC();
        MatC A = MatC::Zero(n + 1, n + 1);
        A.topLeftCorner(n, n) = Mm;
        A.topRightCorner(n, 1) = sourceOf(Mm, MB, e.bgl, bgFieldAt(e.bgl, e.y0n + yl));
        return A;
      };
      auto nodesAt = [&](std::size_t k) {
        std::vector<MatC> nodes;
        for (int s2 = 0; s2 < tab.stages(); ++s2) nodes.push_back(augAt(e.steps.y[k] + tab.c(s2) * e.steps.h[k]));
        return nodes;
      };
      const std::size_t ms = e.steps.h.size();
      std::vector<MatC> Pk(ms), Lk(ms);
      for (std::size_t k = 0; k < ms; ++k) {
        Pk[k] = rkStepPropagator(nodesAt(k), e.steps.h[k], tab);
        Lk[k] = k == 0 ? Pk[k] : MatC(Pk[k] * Lk[k - 1]);
      }
      MatC G = MatC::Zero(n + 1, n + 1);
      G.topLeftCorner(n, n) = ref.PsiInv.adjoint() * Trefb * ref.Psi.adjoint();
      G.topRightCorner(n, 1) = *fbar;
      const double tk = model.thickness_um();
      double dLdH = 0;
      MatC Msum = MatC::Zero(n, n);
      const BgLayerField& L = bgf[e.bgl];
      for (std::size_t k = ms; k-- > 0;) {
        const MatC Pb = k == 0 ? G : MatC(G * Lk[k - 1].adjoint());
        std::vector<MatC> nb;
        double hb = 0;
        adjoint::rkStep(nodesAt(k), e.steps.h[k], tab, Pb, nb, &hb);
        dLdH += hb * e.steps.h[k];  // h_k proportional to the thickness
        for (int s2 = 0; s2 < tab.stages(); ++s2) {
          const double yl = e.steps.y[k] + tab.c(s2) * e.steps.h[k], y_um = yl / k0, y = e.y0n + yl;
          const MatC& N2 = nb[static_cast<std::size_t>(s2)];
          MatC Mb = N2.topLeftCorner(n, n);
          const VecC sigb = N2.topRightCorner(n, 1);
          // sigma = T(s), s = sum_q w_q (M - M_bg)_{:, q nf + c}, w = sum_k Psi_k amp_k E_k(y)
          VecC Ek;
          const VecC w = bgFieldAt(e.bgl, y, &Ek);
          const QBlocks Qy = model.Q(e.constant ? 0.0 : y_um);
          const MatC Ms = bloch ? buildPropagationMatrix(Qy, alB, gaB) : MatC(e.constant ? e.Mc : MatC(Mof(Qy)));
          const MatC& Lm = bloch ? bgMB[e.bgl] : bgM[e.bgl];
          MatC& Lmb = bloch ? bgMBb[e.bgl] : bgMb[e.bgl];
          if (Lmb.size() == 0) Lmb = MatC::Zero(n, n);
          VecC spre = VecC::Zero(n);
          for (int q = 0; q < 4; ++q) spre += w(q) * VecC(Ms.col(q * nf + c) - Lm.col(q * nf + c));
          if (bloch)
            for (int q = 0; q < 4; ++q) TsrcB += sigb.segment(q * nf, nf) * spre.segment(q * nf, nf).adjoint();
          const VecC sb = bloch ? detail::componentwise(MatC(Tsrc.adjoint()), sigb) : sigb;
          VecC wb(4);
          MatC MBb = bloch ? MatC(MatC::Zero(n, n)) : MatC();
          for (int q = 0; q < 4; ++q) {
            wb(q) = VecC(Ms.col(q * nf + c) - Lm.col(q * nf + c)).dot(sb);
            (bloch ? MBb : Mb).col(q * nf + c) += std::conj(w(q)) * sb;
            Lmb.col(q * nf + c) -= std::conj(w(q)) * sb;
          }
          if (bloch) {
            QBlocks Qb;
            if (global) adjoint::propagationMatrix(Qy, alB, gaB, MBb, Qb, alBb, gaBb);
            else adjoint::propagationMatrix(Qy, alB, gaB, MBb, Qb);
            if (own) fromQbar(e.layer, false, y_um, Qb);
          }
          double yb = 0;  // the node height
          for (int k2 = 0; k2 < 4; ++k2) {
            cplx Eb2 = 0, ab = 0;
            for (int q = 0; q < 4; ++q) {
              Eb2 += std::conj(L.Psi(q, k2) * L.amp(k2)) * wb(q);
              ab += std::conj(L.Psi(q, k2) * Ek(k2)) * wb(q);
              bgPsiB[e.bgl](q, k2) += wb(q) * std::conj(L.amp(k2) * Ek(k2));
            }
            bgAmpB[e.bgl](k2) += ab;
            const double yr = k2 < 2 ? L.y1 : L.y0;
            bgKapB[e.bgl](k2) += Eb2 * std::conj(Ek(k2) * (y - yr));
            const double d = std::real(std::conj(Eb2) * Ek(k2) * L.lambda(k2));
            yb += d;
            (k2 < 2 ? bgY1B : bgY0B)[e.bgl] -= d;
          }
          dLdlnk0 += yb * y;
          for (std::size_t q = 0; q < ps.size(); ++q) {
            if (ps[q].kind != Param::Kind::Thickness) continue;
            if (ps[q].layer < e.layer) out.grad[q] += k0 * yb;
            else if (ps[q].layer == e.layer) out.grad[q] += yb * yl / tk;  // the steps scale with the thickness
          }
          if (e.constant) Msum += Mb;
          else fromMbar(Qy, Mb, y_um, false);
        }
        G = Pk[k].adjoint() * G;
      }
      thicknessParam(e.layer, dLdH / tk);
      dLdlnk0 += dLdH;
      if (e.constant) fromMbar(model.Q(0.0), Msum, 0.0, false);
      return;
    }
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
          nodes.push_back(e.constant ? e.Mc : MatC(Mof(model.Q((e.steps.y[k] + tab.c(s) * e.steps.h[k]) / k0))));
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
    MatC Mb = adjoint::eig(W, lam, Wb, lamb, cls);
    if (Msrc) Mb += *Msrc;
    fromMbar(e.Q, Mb, e.y_um, e.homogeneous);
    if (MBsrc) {  // the sources in the Bloch basis
      QBlocks Qb;
      if (global) adjoint::propagationMatrix(e.Q, alB, gaB, *MBsrc, Qb, alBb, gaBb);
      else adjoint::propagationMatrix(e.Q, alB, gaB, *MBsrc, Qb);
      if (own) fromQbar(e.layer, false, e.y_um, Qb);
    }
  };
  /// source amplitudes su, sv of element i (cotangents ub, vb): S-matrix cotangents into Sb, M (and the Bloch-basis
  /// M) cotangents of the particular solution into Msrc, MBsrc
  auto sourceBackward = [&](std::size_t i, const VecC& ub, const VecC& vb, SMatrix& Sb, MatC& Msrc, MatC& MBsrc) {
    const Elem& e = el[i];
    const VecC g0 = ref.PsiInv * e.F0, g1 = ref.PsiInv * e.F1;
    Sb.S11 -= ub * g0.tail(nh).adjoint();
    Sb.S12 -= ub * g1.head(nh).adjoint();
    Sb.S21 -= vb * g0.tail(nh).adjoint();
    Sb.S22 -= vb * g1.head(nh).adjoint();
    VecC g0b(2 * nh), g1b(2 * nh);
    g0b << vb, -(e.S.S11.adjoint() * ub) - e.S.S21.adjoint() * vb;
    g1b << -(e.S.S12.adjoint() * ub) - e.S.S22.adjoint() * vb, ub;
    if (global) adjoint::accumulate(refPsiInvb, MatC(g0b * e.F0.adjoint() + g1b * e.F1.adjoint()));
    const VecC F0b = ref.PsiInv.adjoint() * g0b, F1b = ref.PsiInv.adjoint() * g1b;
    const MatC M = Mof(e.Q);
    const Eigen::Index N = M.rows();
    const BgLayerField& L = bgf[e.bgl];
    Msrc = MatC::Zero(N, N);
    MBsrc = bloch ? MatC(MatC::Zero(N, N)) : MatC();
    const MatC MB = bloch ? buildPropagationMatrix(e.Q, alB, gaB) : MatC();
    const MatC& Ms = bloch ? MB : M;
    const MatC& Lm = bloch ? bgMB[e.bgl] : bgM[e.bgl];
    MatC& Lmb = bloch ? bgMBb[e.bgl] : bgMb[e.bgl];
    if (Lmb.size() == 0) Lmb = MatC::Zero(N, N);
    double y0b = 0, y1b = 0;  // slice boundaries
    for (std::size_t j = 0; j < e.pv.size(); ++j) {
      const int k = e.pk[j];
      const cplx kap = L.lambda(k);
      const VecC& pv = e.pv[j];
      const VecC pb = std::conj(e.e0[j]) * F0b + std::conj(e.e1[j]) * F1b;
      const VecC rb = (kap * MatC::Identity(N, N) - M).adjoint().partialPivLu().solve(pb);  // r = sigma amp
      Msrc += rb * pv.adjoint();
      // e = exp(kappa (y - y_ref)): exponent, slice boundaries and the background-layer boundary
      const double yr = k < 2 ? L.y1 : L.y0;
      const cplx e0b = pv.dot(F0b), e1b = pv.dot(F1b);
      cplx kapb = e0b * std::conj(e.e0[j] * (e.y0n - yr)) + e1b * std::conj(e.e1[j] * (e.y1n - yr)) - pv.dot(rb);
      const double d0 = std::real(std::conj(e0b) * e.e0[j] * kap), d1 = std::real(std::conj(e1b) * e.e1[j] * kap);
      y0b += d0;
      y1b += d1;
      (k < 2 ? bgY1B : bgY0B)[e.bgl] -= d0 + d1;
      bgKapB[e.bgl](k) += kapb;
      // r = sigma amp, sigma = T(s), s = sum_q w_q (M - M_bg)_{:, q nf + c}
      VecC spre = VecC::Zero(N);
      for (int q = 0; q < 4; ++q) spre += L.Psi(q, k) * VecC(Ms.col(q * nf + c) - Lm.col(q * nf + c));
      const VecC sig = bloch ? detail::componentwise(Tsrc, spre) : spre;
      bgAmpB[e.bgl](k) += sig.dot(rb);
      const VecC sigb = std::conj(L.amp(k)) * rb;
      if (bloch)
        for (int q = 0; q < 4; ++q) TsrcB += sigb.segment(q * nf, nf) * spre.segment(q * nf, nf).adjoint();
      const VecC sb = bloch ? detail::componentwise(MatC(Tsrc.adjoint()), sigb) : sigb;
      MatC& Mt = bloch ? MBsrc : Msrc;
      for (int q = 0; q < 4; ++q) {
        bgPsiB[e.bgl](q, k) += VecC(Ms.col(q * nf + c) - Lm.col(q * nf + c)).dot(sb);
        Mt.col(q * nf + c) += std::conj(L.Psi(q, k)) * sb;
        Lmb.col(q * nf + c) -= std::conj(L.Psi(q, k)) * sb;
      }
    }
    // slice boundaries y0n = sum_{l<k} k0 t_l + k0 t_k s/K (normalized: ~ k0)
    dLdlnk0 += y0b * e.y0n + y1b * e.y1n;
    for (std::size_t q = 0; q < ps.size(); ++q) {
      if (ps[q].kind != Param::Kind::Thickness) continue;
      if (ps[q].layer < e.layer) {
        out.grad[q] += k0 * (y0b + y1b);
      } else if (ps[q].layer == e.layer) {
        const double tk = models[e.layer]->thickness_um();
        out.grad[q] += (y0b * (e.y0n - layerBottom[e.layer]) + y1b * (e.y1n - layerBottom[e.layer])) / tk;
      }
    }
  };
  // reverse sweep over the star chain, per checkpoint segment (recompute the dropped slices and prefixes)
  SMatrix E0bar;
  {
    SMatrix Pbar;
    Pbar.S11 = MatC::Zero(nh, nh);
    Pbar.S21 = MatC::Zero(nh, nh);
    Pbar.S12 = bbar * a_in.adjoint();
    Pbar.S22 = tbar * a_in.adjoint();
    VecC ub = bbar, vb = tbar;  // scattered field: cotangents of the prefix source vectors
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
        if (scat) {
          VecC Aub, Avb, Bub, Bvb;
          const VecC z = VecC::Zero(nh);
          adjoint::starSources(P[j - 1], el[j].S, Pu[j - 1], el[j].sv.size() ? el[j].sv : z, ub, vb, Ab, Bb, Aub, Avb, Bub, Bvb);
          ub = std::move(Aub);
          vb = std::move(Avb);
          if (el[j].su.size() && el[j].rk) {  // affine sources u = g_b - S12 g_a, v = -S22 g_a, g = Psi^-1 f
            const Elem& e = el[j];
            const VecC g = ref.PsiInv * e.F0;
            Bb.S12 -= Bub * g.head(nh).adjoint();
            Bb.S22 -= Bvb * g.head(nh).adjoint();
            VecC gb(2 * nh);
            gb << -(e.S.S12.adjoint() * Bub) - e.S.S22.adjoint() * Bvb, Bub;
            if (global) adjoint::accumulate(refPsiInvb, MatC(gb * e.F0.adjoint()));
            const VecC fb = ref.PsiInv.adjoint() * gb;
            Pbar = std::move(Ab);
            backward(j, Bb, nullptr, nullptr, &fb);
            continue;
          }
          if (el[j].su.size()) {
            MatC Msrc, MBsrc;
            sourceBackward(j, Bub, Bvb, Bb, Msrc, MBsrc);
            Pbar = std::move(Ab);
            backward(j, Bb, &Msrc, bloch ? &MBsrc : nullptr);
            continue;
          }
        }
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
  if (scat) {
    // background propagation matrices M_bg (X) and the Bloch-basis ones -> background permittivities
    for (std::size_t l = 0; l < nbg; ++l) {
      const cplx eps = st.material(bgf[l].material).eps;
      const QBlocks Q = homogeneousIsotropicQ(eps, nf);
      QBlocks Qb;
      if (bgMb[l].size()) {
        if (trans && global) adjoint::propagationMatrix(Q, tc.A, tc.G, bgMb[l], Qb, Abar, Gbar);
        else if (trans) adjoint::propagationMatrix(Q, tc.A, tc.G, bgMb[l], Qb);
        else if (global) adjoint::propagationMatrix(Q, al, ga, bgMb[l], Qb, alb, gab);
        else adjoint::propagationMatrix(Q, al, ga, bgMb[l], Qb);
        epsBar[bgf[l].material] += Qb(0, 0).trace() + Qb(1, 1).trace() + Qb(2, 2).trace();
      }
      if (bgMBb[l].size()) {
        QBlocks Qb2;
        if (global) adjoint::propagationMatrix(Q, alB, gaB, bgMBb[l], Qb2, alBb, gaBb);
        else adjoint::propagationMatrix(Q, alB, gaB, bgMBb[l], Qb2);
        epsBar[bgf[l].material] += Qb2(0, 0).trace() + Qb2(1, 1).trace() + Qb2(2, 2).trace();
      }
    }
    // the background solution (N = 0) and Tsrc: central differences per parameter acting on them
    std::vector<std::string> bgMats{st.superstrate, st.substrate};
    for (const auto& L : bgf) bgMats.push_back(L.material);
    for (std::size_t q = 0; q < ps.size(); ++q) {
      const Param& p = ps[q];
      const bool acts = p.kind == Param::Kind::Thickness || p.kind == Param::Kind::Lambda || p.kind == Param::Kind::Phi ||
                        p.kind == Param::Kind::Psi || p.kind == Param::Kind::Delta ||
                        (p.kind == Param::Kind::Material && std::find(bgMats.begin(), bgMats.end(), p.material) != bgMats.end());
      if (!acts) continue;
      const double h = (p.kind == Param::Kind::Thickness || p.kind == Param::Kind::Lambda) ? 1e-6 : 1e-5;
      double val = 0;
      for (const double sg : {1.0, -1.0}) {
        Stack s2;
        Wavelength w2;
        perturbedSetup(st, wl, p, sg * h, s2, w2);
        const Material& sp2 = s2.material(s2.superstrate);
        const IncidentWave i2 = makeIncidence(std::sqrt(sp2.eps.real()), w2.phi, w2.psi, w2.delta);
        const auto b2 = backgroundField(s2, opt, w2, i2.alpha0, i2.gamma0);
        FG_REQUIRE(b2.size() == nbg, "gradient: the parameter changes the background layers");
        double v = 0;
        for (std::size_t l = 0; l < nbg; ++l) {
          v += std::real(bgAmpB[l].dot(b2[l].amp) + bgKapB[l].dot(b2[l].lambda)) + bgY0B[l] * b2[l].y0 + bgY1B[l] * b2[l].y1;
          v += std::real((bgPsiB[l].adjoint() * b2[l].Psi).trace());
        }
        if (bloch) {
          const double k2 = 2 * PI / w2.lambda_um;
          v += std::real((TsrcB.adjoint() * makeTsrc(k2 * i2.alpha0, k2 * i2.gamma0)).trace());
        }
        val += sg * v;
      }
      out.grad[q] += val / (2 * h);
    }
  }
  if (global) {
    // first element S0 = fromTransfer(PsiInv_sup Psi_sub)
    MatC Xb;
    adjoint::fromTransfer(MatC(ref.PsiInv * sm.Psi), E0bar, Xb);
    adjoint::accumulate(refPsiInvb, MatC(Xb * sm.Psi.adjoint()));
    adjoint::accumulate(smPsib, MatC(ref.PsiInv.adjoint() * Xb));
    cplx epsSupb = 0, epsSubb = 0;
    // amplitudes b = S12 a_in, t = S22 a_in
    VecC ainb = Plast.S12.adjoint() * bbar + Plast.S22.adjoint() * tbar;
    if (scat) ainb.setZero();  // the sources drive the scattered field
    dLdlnk0 += ampBar.dlnk0;
    double dLdkx0 = 0, dLdkz0 = 0;  // transformed coordinates: conversions eq:ct_rayleigh / eq:ct_inc at k0 alpha0
    if (trans) {
      MatC refPhysPsib = MatC::Zero(4 * nf, 4 * nf), refPhysPsiInvb = MatC::Zero(4 * nf, 4 * nf),
           subPhysPsiInvb = MatC::Zero(4 * nf, 4 * nf), toPhysb = MatC::Zero(nf, nf), fromPhysb = MatC::Zero(nf, nf);
      adjoint::accumulate(refPsib, MatC(MatC::Zero(4 * nf, 4 * nf)));
      adjoint::accumulate(refPsiInvb, MatC(MatC::Zero(4 * nf, 4 * nf)));
      adjoint::accumulate(smPsib, MatC(MatC::Zero(4 * nf, 4 * nf)));
      if (ampBar.refPsi.size()) refPsib += ampBar.refPsi;
      if (ampBar.subPsi.size()) smPsib += ampBar.subPsi;
      if (ampBar.aIn.size()) ainb += ampBar.aIn;
      auto cwT = [&](const MatC& T, const VecC& v) { return detail::componentwise(T, v); };
      auto outer = [&](MatC& Tbr, const VecC& yb, const VecC& x) {  // T applied per component: Tbar += sum yb_k x_k^H
        for (int k = 0; k < 4; ++k) Tbr += yb.segment(k * nf, nf) * x.segment(k * nf, nf).adjoint();
      };
      // reflected: bP = (PsiInv_phys toPhys(Psi_up b))_up
      {
        const VecC x = ref.Psi.rightCols(nh) * b, y = cwT(tc.toPhys, x);
        VecC z = VecC::Zero(2 * nh);
        z.tail(nh) = bbarP;
        refPhysPsiInvb += z * y.adjoint();
        const VecC yb = refPhys.PsiInv.adjoint() * z;
        outer(toPhysb, yb, x);
        const VecC xb = cwT(MatC(tc.toPhys.adjoint()), yb);
        refPsib.rightCols(nh) += xb * b.adjoint();
      }
      // transmitted: tP = (PsiInv_phys,sub toPhys(Psi_sub,down t))_down
      {
        const VecC x = sm.Psi.leftCols(nh) * t, y = cwT(tc.toPhys, x);
        VecC z = VecC::Zero(2 * nh);
        z.head(nh) = tbarP;
        subPhysPsiInvb += z * y.adjoint();
        const VecC yb = subPhys.PsiInv.adjoint() * z;
        outer(toPhysb, yb, x);
        const VecC xb = cwT(MatC(tc.toPhys.adjoint()), yb);
        smPsib.leftCols(nh) += xb * t.adjoint();
      }
      // incidence: a_in = mask (PsiInv_sup fromPhys(Psi_phys,down a_ph))_down
      VecC aphb;
      {
        const VecC v = refPhys.Psi.leftCols(nh) * a_ph, w = cwT(tc.fromPhys, v);
        VecC z = VecC::Zero(2 * nh);
        for (int q = 0; q < nh; ++q)
          if (!(ref.lambda(q).real() > 1e-8)) z(q) = ainb(q);
        refPsiInvb += z * w.adjoint();
        const VecC wb = ref.PsiInv.adjoint() * z;
        outer(fromPhysb, wb, v);
        const VecC vb = cwT(MatC(tc.fromPhys.adjoint()), wb);
        refPhysPsib.leftCols(nh) += vb * a_ph.adjoint();
        aphb = refPhys.Psi.leftCols(nh).adjoint() * vb;
      }
      ainb = (incident || scat) ? VecC(VecC::Zero(nh)) : aphb;  // the physical incident amplitudes (modal: fixed)
      // numerical bases of the half spaces (Plan section "M17"): eigenvectors of M_sup, M_sub in X; structured
      // half spaces (M18): Q of their layer model, back to its materials
      const std::size_t noLayer = st.layers.size();
      auto halfSpace = [&](const ModeBasis& mb, MatC Psib, const MatC& PsiInvb, cplx eps, cplx& epsb,
                           const detail::LayerModel* hm, const QBlocks& Qh, const VecC& lamb) {
        Psib -= mb.PsiInv.adjoint() * PsiInvb * mb.PsiInv.adjoint();  // PsiInv = Psi^-1
        std::vector<int> cls(static_cast<std::size_t>(2 * nh), 0);
        for (int q = nh; q < 2 * nh; ++q) cls[static_cast<std::size_t>(q)] = 1;
        const MatC Mb = adjoint::eig(mb.Psi, mb.lambda, Psib, lamb, cls);
        QBlocks Qb;
        adjoint::propagationMatrix(hm ? Qh : homogeneousIsotropicQ(eps, nf), tc.A, tc.G, Mb, Qb, Abar, Gbar);
        if (hm) fromQbar(noLayer, false, 0.0, Qb, hm);
        else epsb += Qb(0, 0).trace() + Qb(1, 1).trace() + Qb(2, 2).trace();
      };
      halfSpace(ref, refPsib, refPsiInvb, sup.eps, epsSupb, halfSup.get(), Qsup, ampBar.refLam);
      halfSpace(sm, smPsib, MatC(MatC::Zero(4 * nf, 4 * nf)), sub.eps, epsSubb, halfSub.get(), Qsub, ampBar.subLam);
      adjoint::isotropicModes(sup.eps, al, ga, refPhys, refPhysPsib, refPhysPsiInvb, VecC(), epsSupb, alb, gab);
      adjoint::isotropicModes(sub.eps, al, ga, subPhys, MatC(), subPhysPsiInvb, VecC(), epsSubb, alb, gab);
      // A = C_x diag(alpha), G = C_z diag(gamma)
      VecC ones = VecC::Ones(nf);
      const auto c1 = detail::transformedCoordinates(map, tr, pmlX.get(), pmlZ.get(), ones, ones, inc.alpha0 * k0, inc.gamma0 * k0);
      if (Abar.size())
        for (int cc = 0; cc < nf; ++cc)
          for (int r = 0; r < nf; ++r) {
            alb(cc) += Abar(r, cc) * std::conj(c1.A(r, cc));
            gab(cc) += Gbar(r, cc) * std::conj(c1.G(r, cc));
          }
      // conversions at kx0 = k0 alpha0, kz0 = k0 gamma0 (central differences of the 1D matrices)
      auto convDot = [&](double hx, double hz) {
        const auto cp = detail::transformedCoordinates(map, tr, pmlX.get(), pmlZ.get(), ones, ones, inc.alpha0 * k0 + hx,
                                                       inc.gamma0 * k0 + hz);
        return std::real((toPhysb.adjoint() * cp.toPhys).trace() + (fromPhysb.adjoint() * cp.fromPhys).trace());
      };
      constexpr double hk = 1e-6;
      if (tr.x) dLdkx0 = (convDot(hk, 0) - convDot(-hk, 0)) / (2 * hk);
      if (tr.z) dLdkz0 = (convDot(0, hk) - convDot(0, -hk)) / (2 * hk);
    } else {
      adjoint::isotropicModes(sup.eps, al, ga, ref, refPsib, refPsiInvb, VecC(), epsSupb, alb, gab);
      adjoint::isotropicModes(sub.eps, al, ga, sm, smPsib, MatC(), VecC(), epsSubb, alb, gab);
    }
    // directional derivatives w.r.t. the scalars the bases, M and the amplitudes depend on
    double dAlpha0 = (oblX ? 0.0 : alb.real().sum()) + alBb.real().sum() + k0 * dLdkx0,
           dGamma0 = (oblZ ? 0.0 : gab.real().sum()) + gaBb.real().sum() + k0 * dLdkz0, dKx = 0, dKz = 0;
    dLdlnk0 += inc.alpha0 * k0 * dLdkx0 + inc.gamma0 * k0 * dLdkz0;
    for (int r = 0; r < nf; ++r) {
      const auto [n, m] = map.orders(r);
      dKx += n * (alb(r).real() + alBb(r).real());
      dKz += m * (gab(r).real() + gaBb(r).real());
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
          const auto ro = planeWaveOrder(g.epsSup, a(r), gg(r), betaBranch(g.epsSup, a(r), gg(r)), +1, bP(r), bP(nf + r));
          if (ro.propagating) v += w * ro.Sy / S;
        }
        if (!subLossless) continue;
        if (const double w = weight('T', n, m); w != 0) {
          const auto to = planeWaveOrder(g.epsSub, a(r), gg(r), betaBranch(g.epsSub, a(r), gg(r)), -1, tP(r), tP(nf + r));
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
        case Param::Kind::Lambda: out.grad[q] += (dKx * Kx + dKz * Kz - dLdlnk0) / wl.lambda_um; break;
        case Param::Kind::Phi: out.grad[q] += incidenceDot(0, 1, 0, 0) * deg; break;
        case Param::Kind::Psi: out.grad[q] += incidenceDot(0, 0, 1, 0) * deg; break;
        case Param::Kind::Delta: out.grad[q] += incidenceDot(0, 0, 0, 1) * deg; break;
        default: break;
      }
    }
  }
  for (std::size_t q = 0; q < ps.size(); ++q) {
    const Param& p = ps[q];
    if (p.kind != Param::Kind::Material) continue;
    const auto tb = tensorBar.find(p.material);
    if (p.tensor) {  // dL = Re sum conj(Tbar_ij) dE_ij
      if (tb != tensorBar.end()) out.grad[q] += std::real((tb->second.conjugate().cwiseProduct(p.dE)).sum());
      continue;
    }
    // isotropic material in anisotropic layers: eps 1 -> trace of the tensor cotangent
    const cplx eb = epsBar[p.material] + (tb != tensorBar.end() ? tb->second.trace() : cplx(0));
    if (!p.viaN) {
      out.grad[q] += p.comp == 0 ? eb.real() : eb.imag();
    } else {
      const cplx nidx = std::sqrt(st.material(p.material).eps);  // eps = n^2
      const cplx deps = p.comp == 0 ? 2.0 * nidx : 2.0 * nidx * I;
      out.grad[q] += std::real(std::conj(eb) * deps);
    }
  }
  out.notes = {"RCWA slices fixed (Plan M19 \"Methode\" 1.)",
               "normal field: profiles differentiated, Algorithm NF over one raster cell, CAD sections frozen (Plan M19 "
               "\"Methode\" 2., v2.12)"};
  return out;
}

/// Modal excitation (Plan section "Gradienten aperiodischer Strukturen", Stufe A): a guided mode of the
/// structured superstrate comes in; L = sum w |r_j|^2 + sum w |t_j|^2 with the powers of the solver
/// (|b_k|^2 P_k / P_in, flux over the physical region). The selection of the guided modes is fixed; the powers
/// do not depend on phase and scale of the eigenvectors, so the eigenvector reverse rule applies.
GradientResult modalGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                             const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps) {
  FG_PROF_SCOPE("gradient.modal");
  const AperiodicSpec& ap = st.aperiodic;
  FG_REQUIRE(objective.scattering.empty(), "gradient: modal excitation: use mode objectives (|r_j|^2, |t_j|^2)");
  for (const auto& p : ps)
    FG_REQUIRE(p.kind != Param::Kind::Phi && p.kind != Param::Kind::Psi && p.kind != Param::Kind::Delta,
               "gradient: modal excitation ignores the incidence angles");
  const IndexMap map(opt.Nx, st.is2D() ? opt.Nz : 0);
  const int nf = map.size(), nh = 2 * nf;
  const double x0 = ap.x.active() ? ap.x.width_um : 0.0, x1 = ap.x.active() ? st.dx_um - ap.x.width_um : st.dx_um;
  const double z0 = ap.z.active() ? ap.z.width_um : 0.0, z1 = ap.z.active() ? st.dz_um - ap.z.width_um : st.dz_um;
  auto axisInt = [](double kappa, double a, double b) -> cplx {
    if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
    return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
  };
  MatC W(nf, nf);
  for (int r = 0; r < nf; ++r)
    for (int cc = 0; cc < nf; ++cc) {
      const auto [nr, mr] = map.orders(r);
      const auto [nc, mc] = map.orders(cc);
      W(r, cc) = axisInt(2 * PI * (nc - nr) / st.dx_um, x0, x1);
      if (st.is2D()) W(r, cc) *= axisInt(2 * PI * (mc - mr) / st.dz_um, z0, z1);
    }
  auto flux = [&](const VecC& F) {
    return 0.5 * std::real(F.segment(2 * nf, nf).dot(W * F.segment(nf, nf)) - F.segment(3 * nf, nf).dot(W * F.head(nf)));
  };
  auto fluxBar = [&](const VecC& F) {
    VecC Fb(4 * nf);
    Fb.head(nf) = -0.5 * (W * F.segment(3 * nf, nf));
    Fb.segment(nf, nf) = 0.5 * (W * F.segment(2 * nf, nf));
    Fb.segment(2 * nf, nf) = 0.5 * (W * F.segment(nf, nf));
    Fb.segment(3 * nf, nf) = -0.5 * (W * F.head(nf));
    return Fb;
  };
  // guided modes of a basis among the columns [first, first + nh) (criteria of the solver, Plan section "M18")
  const double cell = st.dx_um * (st.is2D() ? st.dz_um : 1.0);
  auto guided = [&](const ModeBasis& b, int first, const std::string& clad) {
    const double nclad = std::sqrt(std::max(st.material(clad).tensor()(0, 0).real(), 0.0));
    std::vector<std::pair<double, int>> g;
    for (int k = first; k < first + nh; ++k) {
      const cplx l = b.lambda(k);
      const VecC& F = b.Psi.col(k);
      const VecC Ex = F.head(nf), Ez = F.segment(nf, nf);
      const double in = std::real(Ex.dot(W * Ex) + Ez.dot(W * Ez)) / (cell * (Ex.squaredNorm() + Ez.squaredNorm()));
      if (std::abs(l.real()) <= 1e-3 * std::abs(l.imag()) && std::abs(l.imag()) > nclad * (1 + 1e-9) && in >= 0.9)
        g.push_back({std::abs(l.imag()), k});
    }
    std::sort(g.begin(), g.end(), [](const auto& a, const auto& c) { return a.first > c.first; });
    std::vector<int> out;
    for (const auto& [n, k] : g) out.push_back(k);
    return out;
  };
  auto describe = [&](const ModeBasis& b, int col, int index) {
    ModeCoefficient m;
    m.index = index;
    m.nEff = cplx(std::abs(b.lambda(col).imag()), std::abs(b.lambda(col).real()));
    const VecC E = b.Psi.col(col);
    const double ex = E.segment(0, nf).squaredNorm(), ez = E.segment(nf, nf).squaredNorm();
    m.teFraction = ex + ez > 0 ? ez / (ex + ez) : 0.0;
    return m;
  };
  int cin = -1, jin = 0;
  const IncidentHook pick = [&](const ModeBasis& ref) {
    const std::vector<int> down = guided(ref, 0, st.superstrate);
    FG_REQUIRE(!down.empty(), "modal excitation: the structured superstrate has no guided mode");
    if (ap.nEff > 0) {
      for (int k = 1; k < static_cast<int>(down.size()); ++k)
        if (std::abs(std::abs(ref.lambda(down[static_cast<std::size_t>(k)]).imag()) - ap.nEff) <
            std::abs(std::abs(ref.lambda(down[static_cast<std::size_t>(jin)]).imag()) - ap.nEff))
          jin = k;
    } else {
      FG_REQUIRE(ap.mode >= 0 && ap.mode < static_cast<int>(down.size()),
                 "modal excitation: mode " + std::to_string(ap.mode) + " not guided");
      jin = ap.mode;
    }
    cin = down[static_cast<std::size_t>(jin)];
    VecC a = VecC::Zero(nh);
    a(cin) = 1.0;
    return a;
  };
  auto weight = [&](char side, int index) {
    double w = 0;
    for (const auto& [m, v] : objective.modes)
      if (m.side == side && m.index == index) w += v;
    return w;
  };
  const AmpObjective modes = [&](const AmpContext& cx, const VecC& b, const VecC& t, VecC& bbar, VecC& tbar, AmpBar& ab,
                                 Result& res) {
    const VecC& Fin = cx.ref.Psi.col(cin);
    const double Pin = -flux(Fin);
    FG_REQUIRE(Pin > 0, "modal excitation: the incident mode carries no power through the physical region");
    ModalResult& mr = res.modal;
    mr.present = true;
    mr.incident = describe(cx.ref, cin, jin);
    mr.incident.amplitude = 1.0;
    mr.incident.power = 1.0;
    ab.refPsi = MatC::Zero(4 * nf, 4 * nf);
    ab.subPsi = MatC::Zero(4 * nf, 4 * nf);
    double L = 0, sum = 0;
    int used = 0;
    auto side = [&](const ModeBasis& mb, const std::vector<int>& cols, int offset, double sign, const VecC& amp, VecC& ampBar,
                    MatC& PsiBar, char s, std::vector<ModeCoefficient>& list) {
      for (std::size_t k = 0; k < cols.size(); ++k) {
        const int col = cols[k], i = col - offset;
        const VecC& F = mb.Psi.col(col);
        const double Pk = sign * flux(F), power = std::norm(amp(i)) * Pk / Pin;
        ModeCoefficient m = describe(mb, col, static_cast<int>(k));
        m.amplitude = amp(i) * std::sqrt(std::abs(Pk) / Pin);
        m.power = power;
        sum += power;
        list.push_back(m);
        const double w = weight(s, static_cast<int>(k));
        if (w == 0) continue;
        ++used;
        L += w * power;
        ampBar(i) += 2 * w * Pk / Pin * amp(i);
        PsiBar.col(col) += (sign * w * std::norm(amp(i)) / Pin) * fluxBar(F);
      }
    };
    side(cx.ref, guided(cx.ref, nh, st.superstrate), nh, 1.0, b, bbar, ab.refPsi, 'R', mr.reflected);
    if (st.substrateStructure)
      side(cx.sub, guided(cx.sub, 0, st.substrate), 0, -1.0, t, tbar, ab.subPsi, 'T', mr.transmitted);
    FG_REQUIRE(used == static_cast<int>(objective.modes.size()) || objective.modes.empty(),
               "gradient: a mode of the objective is not guided");
    mr.radiated = 1.0 - sum;
    ab.refPsi.col(cin) += (L / Pin) * fluxBar(Fin);  // P_in = -q_W(psi_in)
    return L;
  };
  GradientResult out = gradientCore(st, opt, wl, prebuilt, objective, std::move(ps), &modes, &pick);
  out.notes.push_back("modal excitation: the selection of the guided modes is fixed");
  return out;
}

/// Scattered-field formulation (Plan section "Gradienten aperiodischer Strukturen", Stufe B): the scattered field
/// directly, sampled just above (y = top + 1 nm) and below (y = -1 nm) the stack as the solver;
/// P_up = q_W(F_top) / I_0, P_down = -q_W(F_bottom) / I_0 with the intensity I_0 of the incident plane wave. The
/// sources sigma = (M - M_bg) F_bg of every structured slice are differentiated with respect to M, the background
/// matrices and the background field (Stufe B2: its N = 0 solution and Tsrc by central differences).
GradientResult scatteredGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                                 const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps,
                                 double wUp, double wDown) {
  FG_PROF_SCOPE("gradient.scattered");
  const AperiodicSpec& ap = st.aperiodic;
  FG_REQUIRE(!st.superstrateStructure && !st.substrateStructure,
             "gradient: structured superstrate/substrate with the scattered-field formulation");
  const std::vector<Param> psIn = ps;
  const IndexMap map(opt.Nx, st.is2D() ? opt.Nz : 0);
  const int nf = map.size(), nh = 2 * nf;
  const double x0 = ap.x.active() ? ap.x.width_um : 0.0, x1 = ap.x.active() ? st.dx_um - ap.x.width_um : st.dx_um;
  const double z0 = ap.z.active() ? ap.z.width_um : 0.0, z1 = ap.z.active() ? st.dz_um - ap.z.width_um : st.dz_um;
  auto axisInt = [](double kappa, double a, double b) -> cplx {
    if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
    return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
  };
  MatC W(nf, nf);
  for (int r = 0; r < nf; ++r)
    for (int cc = 0; cc < nf; ++cc) {
      const auto [nr, mr] = map.orders(r);
      const auto [nc, mc] = map.orders(cc);
      W(r, cc) = axisInt(2 * PI * (nc - nr) / st.dx_um, x0, x1);
      if (st.is2D()) W(r, cc) *= axisInt(2 * PI * (mc - mr) / st.dz_um, z0, z1);
    }
  auto flux = [&](const VecC& F) {
    return 0.5 * std::real(F.segment(2 * nf, nf).dot(W * F.segment(nf, nf)) - F.segment(3 * nf, nf).dot(W * F.head(nf)));
  };
  auto fluxBar = [&](const VecC& F) {
    VecC Fb(4 * nf);
    Fb.head(nf) = -0.5 * (W * F.segment(3 * nf, nf));
    Fb.segment(nf, nf) = 0.5 * (W * F.segment(2 * nf, nf));
    Fb.segment(2 * nf, nf) = 0.5 * (W * F.segment(nf, nf));
    Fb.segment(3 * nf, nf) = -0.5 * (W * F.head(nf));
    return Fb;
  };
  // intensity of the incident plane wave (as the solver): |S_y| sqrt(eps_sup) / beta_0
  const Material& sup = st.material(st.superstrate);
  const IncidentWave inc = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
  double I0 = 0;
  {
    const ModeBasis s1 = isotropicModes(sup.eps, VecC::Constant(1, inc.alpha0), VecC::Constant(1, inc.gamma0));
    VecC ainc(2);
    ainc << inc.A(0), inc.A(2);
    const VecC Fi = s1.Psi.leftCols(2) * ainc;
    const double beta0 = std::sqrt(std::max(sup.eps.real() - inc.alpha0 * inc.alpha0 - inc.gamma0 * inc.gamma0, 1e-300));
    I0 = std::abs(0.5 * std::real(Fi(1) * std::conj(Fi(2)) - Fi(0) * std::conj(Fi(3)))) * std::sqrt(sup.eps.real()) / beta0;
  }
  const double k0 = 2 * PI / wl.lambda_um, dTop = 1e-9 * k0, dBot = -1e-9 * k0;
  const AmpObjective scatter = [&](const AmpContext& cx, const VecC& u, const VecC& v, VecC& ubar, VecC& vbar, AmpBar& ab,
                                   Result& res) {
    VecC phU(nh), phD(nh);
    for (int q = 0; q < nh; ++q) {
      phU(q) = std::exp(cx.ref.lambda(nh + q) * dTop);
      phD(q) = std::exp(cx.sub.lambda(q) * dBot);
    }
    const MatC up = cx.ref.Psi.rightCols(nh), down = cx.sub.Psi.leftCols(nh);
    const VecC Fu = up * phU.cwiseProduct(u), Fd = down * phD.cwiseProduct(v);
    ScatteringResult& sc = res.scattering;
    sc.present = true;
    sc.up = flux(Fu) / I0;
    sc.down = -flux(Fd) / I0;
    const VecC Fub = fluxBar(Fu) * (wUp / I0), Fdb = fluxBar(Fd) * (-wDown / I0);
    const VecC Gu = up.adjoint() * Fub, Gd = down.adjoint() * Fdb;
    ubar = phU.conjugate().cwiseProduct(Gu);
    vbar = phD.conjugate().cwiseProduct(Gd);
    // global parameters: the bases, their eigenvalues in the phases e^{lambda d} and d ~ k0
    ab.refPsi = MatC::Zero(4 * nf, 4 * nf);
    ab.subPsi = MatC::Zero(4 * nf, 4 * nf);
    ab.refPsi.rightCols(nh) = Fub * phU.cwiseProduct(u).adjoint();
    ab.subPsi.leftCols(nh) = Fdb * phD.cwiseProduct(v).adjoint();
    ab.refLam = VecC::Zero(2 * nh);
    ab.subLam = VecC::Zero(2 * nh);
    for (int q = 0; q < nh; ++q) {
      const cplx pu = std::conj(u(q)) * Gu(q), pd = std::conj(v(q)) * Gd(q);  // cotangents of the phases
      ab.refLam(nh + q) = std::conj(phU(q) * dTop) * pu;
      ab.subLam(q) = std::conj(phD(q) * dBot) * pd;
      ab.dlnk0 += std::real(std::conj(pu) * phU(q) * cx.ref.lambda(nh + q) * dTop) +
                  std::real(std::conj(pd) * phD(q) * cx.sub.lambda(q) * dBot);
    }
    return wUp * sc.up + wDown * sc.down;
  };
  GradientResult out = gradientCore(st, opt, wl, prebuilt, objective, std::move(ps), &scatter);
  out.notes.push_back("aperiodic: lateral scattering and absorption are not computed by the gradient");
  // I_0 depends on the superstrate and the incidence: dL/dI_0 = -L / I_0 (central differences)
  for (std::size_t q = 0; q < psIn.size(); ++q) {
    const Param& p = psIn[q];
    if (!(p.kind == Param::Kind::Phi || p.kind == Param::Kind::Psi || p.kind == Param::Kind::Delta ||
          (p.kind == Param::Kind::Material && p.material == st.superstrate)))
      continue;
    constexpr double h = 1e-5;
    double d = 0;
    for (const double sg : {1.0, -1.0}) {
      Stack s2;
      Wavelength w2;
      perturbedSetup(st, wl, p, sg * h, s2, w2);
      const Material& sp2 = s2.material(s2.superstrate);
      const IncidentWave i2 = makeIncidence(std::sqrt(sp2.eps.real()), w2.phi, w2.psi, w2.delta);
      const ModeBasis s1 = isotropicModes(sp2.eps, VecC::Constant(1, i2.alpha0), VecC::Constant(1, i2.gamma0));
      VecC ainc(2);
      ainc << i2.A(0), i2.A(2);
      const VecC Fi = s1.Psi.leftCols(2) * ainc;
      const double beta0 = std::sqrt(std::max(sp2.eps.real() - i2.alpha0 * i2.alpha0 - i2.gamma0 * i2.gamma0, 1e-300));
      d += sg * std::abs(0.5 * std::real(Fi(1) * std::conj(Fi(2)) - Fi(0) * std::conj(Fi(3)))) * std::sqrt(sp2.eps.real()) / beta0;
    }
    out.grad[q] += -out.value / I0 * d / (2 * h);
  }
  return out;
}

/// Aperiodic structures, total-field formulation at normal incidence (Plan section "Gradienten aperiodischer
/// Strukturen", Stufe A): P_up = q_W(Psi_up (b - b_bg)) / S_inc, P_down = -q_W(Psi_down (t - t_bg)) / S_inc with the
/// window flux q_W of the solver. Three passes with the same driver: the background stack forward (b_bg, t_bg), the
/// stack with bbar = dL/db, and the background stack with bbar_bg = -bbar (only material and thickness parameters
/// act on it; the structures do not reach into the zones).
GradientResult aperiodicGradient(const Stack& st, const SolverOptions& opt, const Wavelength& wl,
                                 const detail::LayerModels* prebuilt, const Objective& objective, std::vector<Param> ps) {
  FG_PROF_SCOPE("gradient.aperiodic");
  const AperiodicSpec& ap = st.aperiodic;
  FG_REQUIRE(objective.weights.empty(), "gradient: aperiodic structures have no diffraction orders (use scattering)");
  if (ap.modeExcitation) return modalGradient(st, opt, wl, prebuilt, objective, std::move(ps));
  FG_REQUIRE(objective.modes.empty(), "gradient: guided-mode objectives need a modal excitation");
  FG_REQUIRE(!st.superstrateStructure && !st.substrateStructure,
             "gradient: structured superstrate/substrate with plane-wave excitation are not available");
  const Material& sup = st.material(st.superstrate);
  const IncidentWave inc0 = makeIncidence(std::sqrt(sup.eps.real()), wl.phi, wl.psi, wl.delta);
  const bool scattered = opt.scatteredField || (ap.x.active() && std::abs(inc0.alpha0) >= 1e-12) ||
                         (ap.z.active() && std::abs(inc0.gamma0) >= 1e-12);
  double wUp = 0, wDown = 0;
  for (const auto& [r, w] : objective.scattering) {
    FG_REQUIRE(r.kind == ScatterRef::Kind::Up || r.kind == ScatterRef::Kind::Down,
               "gradient: lateral scattering and absorption are not differentiated (Plan v2.13, Stufe C)");
    (r.kind == ScatterRef::Kind::Up ? wUp : wDown) += w;
  }
  if (scattered) return scatteredGradient(st, opt, wl, prebuilt, objective, std::move(ps), wUp, wDown);
  const Stack bg = backgroundStack(st);
  // layer k of the stack -> layer of the background stack (layers of zero thickness are dropped)
  std::vector<std::size_t> bgLayer(st.layers.size(), kNone);
  {
    std::size_t n = 0;
    for (std::size_t k = 0; k < st.layers.size(); ++k) {
      Stack one = st;
      one.layers = {st.layers[k]};
      if (!std::holds_alternative<Profile1DLayer>(st.layers[k].spec)) bgLayer[k] = n;
      n += backgroundStack(one).layers.size();
    }
    FG_REQUIRE(n == bg.layers.size(), "gradient: background stack mismatch");
  }
  std::vector<Param> psBg;
  std::vector<std::size_t> toBg;  // index into ps of each background parameter
  for (std::size_t q = 0; q < ps.size(); ++q) {
    const Param& p = ps[q];
    switch (p.kind) {
      case Param::Kind::Material:
      case Param::Kind::Lambda:
      case Param::Kind::Psi:
      case Param::Kind::Delta:
        psBg.push_back(p);
        toBg.push_back(q);
        break;
      case Param::Kind::Phi:
        throw Error("gradient: phi of an aperiodic structure at normal incidence (phi != 0 switches to the "
                    "scattered-field formulation, Plan v2.13 Stufe B)");
      case Param::Kind::Thickness:
        FG_REQUIRE(bgLayer[p.layer] != kNone, "gradient: thickness of a profile layer with PML");
        psBg.push_back(p);
        psBg.back().layer = bgLayer[p.layer];
        toBg.push_back(q);
        break;
      case Param::Kind::Breakpoint:
      case Param::Kind::Shape:
        break;  // inside the physical region: the background does not move
      default:
        throw Error("gradient: this parameter is not available for aperiodic structures yet (Plan v2.13)");
    }
  }
  // window flux q_W(F) = 1/2 Re(Hx^H W Ez - Hz^H W Ex) over the physical region and its cotangent (as the solver)
  const IndexMap map(opt.Nx, st.is2D() ? opt.Nz : 0);
  const int nf = map.size(), nh = 2 * nf;
  const double x0 = ap.x.active() ? ap.x.width_um : 0.0, x1 = ap.x.active() ? st.dx_um - ap.x.width_um : st.dx_um;
  const double z0 = ap.z.active() ? ap.z.width_um : 0.0, z1 = ap.z.active() ? st.dz_um - ap.z.width_um : st.dz_um;
  auto axisInt = [](double kappa, double a, double b) -> cplx {
    if (std::abs(kappa) * (b - a) < 1e-12) return b - a;
    return (std::exp(I * (kappa * b)) - std::exp(I * (kappa * a))) / (I * kappa);
  };
  MatC W(nf, nf);
  for (int r = 0; r < nf; ++r)
    for (int cc = 0; cc < nf; ++cc) {
      const auto [nr, mr] = map.orders(r);
      const auto [nc, mc] = map.orders(cc);
      W(r, cc) = axisInt(2 * PI * (nc - nr) / st.dx_um, x0, x1);
      if (st.is2D()) W(r, cc) *= axisInt(2 * PI * (mc - mr) / st.dz_um, z0, z1);
    }
  const double area = (x1 - x0) * (st.is2D() ? z1 - z0 : 1.0);
  auto flux = [&](const VecC& F) {
    return 0.5 * std::real(F.segment(2 * nf, nf).dot(W * F.segment(nf, nf)) - F.segment(3 * nf, nf).dot(W * F.head(nf)));
  };
  auto fluxBar = [&](const VecC& F) {  // dq = Re(Fbar^H dF), W Hermitian
    VecC Fb(4 * nf);
    Fb.head(nf) = -0.5 * (W * F.segment(3 * nf, nf));
    Fb.segment(nf, nf) = 0.5 * (W * F.segment(2 * nf, nf));
    Fb.segment(2 * nf, nf) = 0.5 * (W * F.segment(nf, nf));
    Fb.segment(3 * nf, nf) = -0.5 * (W * F.head(nf));
    return Fb;
  };
  // pass 1: background amplitudes
  VecC bb, tb;
  const AmpObjective record = [&](const AmpContext&, const VecC& b, const VecC& t, VecC& bbar, VecC& tbar, AmpBar&,
                                  Result&) {
    bb = b;
    tb = t;
    bbar.setZero();
    tbar.setZero();
    return 0.0;
  };
  gradientCore(bg, opt, wl, nullptr, Objective{}, {}, &record);
  // pass 2: the stack
  VecC cb, ct;
  const AmpObjective scatter = [&](const AmpContext& cx, const VecC& b, const VecC& t, VecC& bbar, VecC& tbar, AmpBar& ab,
                                   Result& res) {
    const VecC Fi = cx.ref.Psi.leftCols(nh) * cx.aIn;
    const double qi = flux(Fi), Sinc = std::abs(qi) / area;
    const MatC up = cx.ref.Psi.rightCols(nh), down = cx.sub.Psi.leftCols(nh);
    const VecC Fu = up * (b - bb), Fd = down * (t - tb);
    ScatteringResult& sc = res.scattering;
    sc.present = true;
    sc.up = flux(Fu) / Sinc;
    sc.down = -flux(Fd) / Sinc;
    const double L = wUp * sc.up + wDown * sc.down;
    const VecC Fub = (wUp / Sinc) * fluxBar(Fu), Fdb = (-wDown / Sinc) * fluxBar(Fd);
    bbar = up.adjoint() * Fub;
    tbar = down.adjoint() * Fdb;
    cb = bbar;
    ct = tbar;
    // through the bases and S_inc = |q_W(Psi_down a_in)| / area (wavelength, polarization, half spaces)
    const VecC Fib = (-L / Sinc * (qi < 0 ? -1.0 : 1.0) / area) * fluxBar(Fi);
    ab.refPsi = MatC::Zero(4 * nf, 4 * nf);
    ab.subPsi = MatC::Zero(4 * nf, 4 * nf);
    ab.refPsi.rightCols(nh) += Fub * (b - bb).adjoint();
    ab.refPsi.leftCols(nh) += Fib * cx.aIn.adjoint();
    ab.subPsi.leftCols(nh) += Fdb * (t - tb).adjoint();
    ab.aIn = cx.ref.Psi.leftCols(nh).adjoint() * Fib;
    return L;
  };
  GradientResult out = gradientCore(st, opt, wl, prebuilt, objective, std::move(ps), &scatter);
  out.notes.push_back("aperiodic: lateral scattering and absorption are not computed by the gradient");
  // pass 3: the background stack with the negated cotangents
  if (!psBg.empty()) {
    const AmpObjective back = [&](const AmpContext&, const VecC&, const VecC&, VecC& bbar, VecC& tbar, AmpBar&, Result&) {
      bbar = -cb;
      tbar = -ct;
      return 0.0;
    };
    const GradientResult g = gradientCore(bg, opt, wl, nullptr, Objective{}, std::move(psBg), &back);
    for (std::size_t i = 0; i < toBg.size(); ++i) out.grad[toBg[i]] += g.grad[i];
  }
  return out;
}

}  // namespace

}  // namespace fg
