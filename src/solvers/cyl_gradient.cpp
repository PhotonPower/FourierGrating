// FourierGrating - gradients of cylindrical structures (M20, Plan sections "Differenzierbarkeit" and V21i): the
// discrete adjoint of the beam computation (total field, RCWA). Per harmonic: powers <- outgoing amplitudes <-
// total S-matrix <- star products <- layer S-matrices <- modes <- M <- Q <- Galerkin matrices <- ring
// permittivities and radii (Plan eq:cyl_shape); thicknesses through the layer S-matrices and the phase of the
// incident beam (referenced at the top of the stack).
#include "cyl_gradient.hpp"

#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <sstream>
#include <variant>

#include "cyl_basis.hpp"
#include "fg/adjoint/primitives.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/cyl/excitation.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::detail {

namespace {

using json = nlohmann::json;

struct CylParam {
  enum class Kind { Thickness, Radius, Material, Shape };
  Kind kind = Kind::Thickness;
  std::size_t layer = 0;  ///< expanded stack
  std::size_t index = 0;  ///< radius index
  std::string material;
  bool viaN = false;
  int comp = 0;  ///< 0: real part, 1: imaginary part
  /// Shape (M21, polygon2d / contour2d in the class path): the boundary velocity of every region of the layer
  std::vector<cyl::BoundaryVelocity> velocity;
};

std::vector<std::string> pointerTokens(const std::string& ptr) {
  std::vector<std::string> t;
  std::stringstream ss(ptr);
  std::string s;
  while (std::getline(ss, s, '/'))
    if (!s.empty()) t.push_back(s);
  return t;
}

CylParam resolveCyl(const json& j, const io::RunConfig& cfg, const std::string& ptr) {
  const json::json_pointer jp(ptr);
  FG_REQUIRE(j.contains(jp) && j.at(jp).is_number(), "gradient: '" + ptr + "' is not a number of the configuration");
  const auto t = pointerTokens(ptr);
  const Stack& st = cfg.stack;
  CylParam p;
  if (t.size() >= 3 && t[0] == "materials") {
    p.kind = CylParam::Kind::Material;
    p.material = t[1];
    FG_REQUIRE(t[2] == "n" || t[2] == "eps", "gradient: material parameters are n or eps ('" + ptr + "')");
    FG_REQUIRE(p.material != st.superstrate && p.material != st.substrate,
               "gradient (cylindrical): materials of the half spaces are not parameters yet ('" + ptr + "')");
    p.viaN = t[2] == "n";
    p.comp = t.size() >= 4 ? std::stoi(t[3]) : 0;
    return p;
  }
  FG_REQUIRE(t.size() >= 3 && t[0] == "layers", "gradient (cylindrical): unsupported parameter '" + ptr + "'");
  for (const auto& l : j["layers"])  // a cad block expands into several layers: the indices would not match
    FG_REQUIRE(l.value("type", std::string()) != "cad", "gradient: layer parameters with cad layers (M19 Phase 2)");
  const std::size_t k = std::stoul(t[1]);
  FG_REQUIRE(k < st.layers.size(), "gradient: layer index out of range in '" + ptr + "'");
  FG_REQUIRE(!std::holds_alternative<RingsProfileLayer>(st.layers[k].spec),
             "gradient (cylindrical): parameters of rings_profile layers are not available ('" + ptr + "')");
  // index in the expanded stack (rings_profile -> staircase, as the solver)
  std::size_t e = 0;
  for (std::size_t i = 0; i < k; ++i)
    if (const auto* rp = std::get_if<RingsProfileLayer>(&st.layers[i].spec))
      e += ringsStaircase(*rp, cfg.solver.rcwaSlices).size();
    else
      ++e;
  p.layer = e;
  if (t[2] == "thickness_um" && t.size() == 3) {
    p.kind = CylParam::Kind::Thickness;
  } else if (t[2] == "radii_um" && t.size() == 4) {
    FG_REQUIRE(std::holds_alternative<RingsLayer>(st.layers[k].spec), "gradient: radii_um needs a rings layer");
    p.kind = CylParam::Kind::Radius;
    p.index = std::stoul(t[3]);
  } else if (t[2] == "regions" || t[2] == "inclusions") {
    FG_REQUIRE(std::holds_alternative<Polygon2DLayer>(st.layers[k].spec) ||
                   std::holds_alternative<Contour2DLayer>(st.layers[k].spec),
               "gradient: shape parameters need a polygon2d or contour2d layer ('" + ptr + "')");
    p.kind = CylParam::Kind::Shape;
  } else {
    throw Error("gradient (cylindrical): unsupported layer parameter '" + ptr + "'");
  }
  return p;
}

/// Background quantities of the exponentials of ring layer l (plane wave): kappa, the field products amplitude x
/// (c_+, c_-, e_y) and the end factors e^{kappa (y_a - y_r)}, e^{kappa (y_b - y_r)}; independent of the
/// normalization of the background modes (Stufe 2b differentiates them by central differences).
struct BgValues {
  cplx kappa, cp, cm, ey, e0, e1;
};
std::vector<BgValues> bgValues(const CylProblem& P, std::size_t l) {
  std::vector<BgValues> out;
  const double ya = P.y0[l], yb = ya + P.layers[l]->thickness_um() * P.k0;
  for (const CylProblem::BgTerm& t : P.bg[l])
    out.push_back({t.kappa, t.amp * t.cp, t.amp * t.cm, t.amp * t.ey, std::exp(t.kappa * (ya - t.yr)),
                   std::exp(t.kappa * (yb - t.yr))});
  return out;
}

/// Boundary velocities of the regions of angular layer l for the shape parameter ptr (M21, V22g): the layers at
/// theta +- delta (delta = 1e-6 max(1, |theta|), as the Cartesian shape parameters), V = (C(+) - C(-)) / (2 delta) per
/// segment and curve parameter. A parameter that breaks the P-fold symmetry fails in buildCylLayers.
std::vector<cyl::BoundaryVelocity> shapeVelocity(const json& j, const std::string& ptr, const io::RunConfig& cfg,
                                                 const CylProblem& P, std::size_t l) {
  const json::json_pointer jp(ptr);
  const double v = j.at(jp).get<double>(), delta = 1e-6 * std::max(1.0, std::abs(v));
  auto at = [&](double d) {
    json a = j;
    a[jp] = v + d;
    const io::RunConfig c = io::parseConfig(a.dump());
    auto layers = buildCylLayers(expandRingsProfiles(c.stack, cfg.solver.rcwaSlices), P.k0);
    FG_REQUIRE(l < layers.size() && layers[l]->angular, "gradient: '" + ptr + "' changes the layers");
    return std::make_shared<const cyl::AngularSection>(layers[l]->section);
  };
  const auto plus = at(delta), minus = at(-delta);
  const cyl::AngularSection& s0 = P.layers[l]->section;
  auto same = [](const geometry::Contour2D& a, const geometry::Contour2D& b) {
    if (a.segments.size() != b.segments.size()) return false;
    for (std::size_t s = 0; s < a.segments.size(); ++s)
      if (a.segments[s].P.size() != b.segments[s].P.size() || a.segments[s].w.size() != b.segments[s].w.size())
        return false;
    return true;
  };
  std::vector<cyl::BoundaryVelocity> out;
  for (std::size_t r = 0; r < s0.regions.size(); ++r) {
    for (const auto* s : {plus.get(), minus.get()}) {
      bool ok = s->regions.size() == s0.regions.size() && same(s->regions[r].outer, s0.regions[r].outer) &&
                s->regions[r].holes.size() == s0.regions[r].holes.size();
      for (std::size_t h = 0; ok && h < s0.regions[r].holes.size(); ++h)
        ok = same(s->regions[r].holes[h], s0.regions[r].holes[h]);
      FG_REQUIRE(ok, "gradient: '" + ptr + "' changes the segment structure of the cross section");
    }
    out.push_back([plus, minus, r, delta](std::size_t c, std::size_t s, double t) {
      const auto& a = plus->regions[r], &b = minus->regions[r];
      const geometry::Point2 p = (c ? a.holes[c - 1] : a.outer).segments[s].point(t),
                             q = (c ? b.holes[c - 1] : b.outer).segments[s].point(t);
      return geometry::Point2{(p[0] - q[0]) / (2 * delta), (p[1] - q[1]) / (2 * delta)};
    });
  }
  return out;
}

/// State of harmonic b (4 N_k) of a class state (4 n_f), and its transpose (scatter into a class state).
VecC blockState(const VecC& F, int b, Eigen::Index n) {
  const Eigen::Index nf = F.size() / 4;
  VecC out(4 * n);
  for (Eigen::Index u = 0; u < 4; ++u) out.segment(u * n, n) = F.segment(u * nf + b * n, n);
  return out;
}
void addBlockState(VecC& F, const VecC& f, int b, Eigen::Index n) {
  const Eigen::Index nf = F.size() / 4;
  for (Eigen::Index u = 0; u < 4; ++u) F.segment(u * nf + b * n, n) += f.segment(u * n, n);
}

}  // namespace

GradientResult cylindricalGradient(const std::string& configJson, const io::RunConfig& cfg, const Wavelength& wl,
                                   const Objective& objective, const std::vector<std::string>& params) {
  FG_PROF_SCOPE("gradient.cylindrical");
  FG_REQUIRE(objective.weights.empty() && objective.scattering.empty() && objective.modes.empty() &&
                 !objective.cylindrical.empty(),
             "gradient (cylindrical): the objective takes the powers R and T (side \"C:R\", \"C:T\", per harmonic "
             "\"C:RM\", \"C:TM\")");
  FG_REQUIRE(cfg.solver.kind == SolverKind::RCWA, "gradient (cylindrical): RCWA only (V21i)");
  const json j = json::parse(configJson);
  std::vector<CylParam> ps;
  for (const auto& ptr : params) ps.push_back(resolveCyl(j, cfg, ptr));

  CylProblem P = cylProblem(cfg.solver, cfg.stack, wl);
  FG_REQUIRE(!P.beamScattered, "gradient (cylindrical): beams in the total-field formulation only (scattered_field false)");
  const bool pw = P.planeWave;
  const Eigen::Index n = P.grid.k.size();
  const double k0 = P.k0;
  // material names of the (expanded) layers: homogeneous -> one, rings -> inner to outer, polygon2d / contour2d ->
  // the regions, the background last
  std::vector<std::vector<std::string>> names;
  for (const LayerEntry& le : P.stack.layers) {
    if (const auto* h = std::get_if<HomogeneousLayer>(&le.spec)) {
      names.push_back({h->material});
    } else if (const auto* r = std::get_if<RingsLayer>(&le.spec)) {
      names.push_back(r->materials);
    } else if (const auto* pg = std::get_if<Polygon2DLayer>(&le.spec)) {
      names.emplace_back();
      for (const auto& inc : pg->inclusions) names.back().push_back(inc.material);
      names.back().push_back(pg->background);
    } else {
      const auto& c = std::get<Contour2DLayer>(le.spec);
      names.emplace_back();
      for (const auto& reg : c.regions) names.back().push_back(reg.material);
      names.back().push_back(c.background);
    }
  }
  for (std::size_t q = 0; q < ps.size(); ++q)
    if (ps[q].kind == CylParam::Kind::Shape) ps[q].velocity = shapeVelocity(j, params[q], cfg, P, ps[q].layer);
  // plane wave: thicknesses and background materials (homogeneous layers, outer ring materials) change the
  // background field (Stufe 2b: its quantities by central differences, the rest analytic)
  std::vector<bool> moveBg(ps.size(), false);
  for (std::size_t q = 0; q < ps.size(); ++q) {
    const CylParam& p = ps[q];
    if (p.kind == CylParam::Kind::Radius)
      FG_REQUIRE(p.index < P.layers[p.layer]->rings.radii.size(), "gradient: radius index out of range");
    if (!pw) continue;
    if (p.kind == CylParam::Kind::Thickness) moveBg[q] = true;
    if (p.kind == CylParam::Kind::Material)
      for (const auto& nm : names) moveBg[q] = moveBg[q] || nm.back() == p.material;
  }
  // cotangents of the background quantities per ring layer and exponential (summed over the harmonics)
  std::vector<std::vector<BgValues>> bgBar(P.layers.size());
  for (std::size_t l = 0; l < P.layers.size(); ++l) bgBar[l].assign(P.bg[l].size(), BgValues{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

  // normalization: beam -> 1 / P_inc (the phase of the beam does not change P_inc), plane wave -> um^2
  double norm = P.toArea;
  if (!pw) {
    double pInc = 0;
    const ModeBasis ref = cyl::homogeneousModes(P.epsSup, P.grid.k);
    for (std::size_t u = 0; u < P.classes.size(); ++u) {
      if (P.aIn[u].isZero(0.0)) continue;
      const int nb = P.classes[u].blocks();
      for (int b = 0; b < nb; ++b) {  // per harmonic (the class state is block-wise that of its harmonics)
        VecC V = VecC::Zero(4 * n);
        V.segment(0, n) = P.aIn[u].segment(b * n, n);
        V.segment(n, n) = P.aIn[u].segment((nb + b) * n, n);
        pInc -= cyl::flux(P.grid, ref.Psi * V);
      }
    }
    FG_REQUIRE(pInc > 0, "cylindrical: the incident beam carries no power (harmonics?)");
    norm = 1.0 / pInc;
  }
  auto weightOf = [&](char side, int m) {
    double w = 0;
    for (const auto& [r, v] : objective.cylindrical)
      if (r.side == side && (r.total || r.m == m)) w += v;
    return w * norm;
  };

  GradientResult out;
  out.grad.assign(ps.size(), 0.0);
  std::map<std::string, cplx> epsBar;
  std::vector<double> hBar(P.layers.size(), 0.0);                  // d L / d (normalized thickness)
  std::vector<std::vector<double>> rBar(P.layers.size());         // d L / d (normalized radius)
  for (std::size_t l = 0; l < P.layers.size(); ++l) rBar[l].assign(P.layers[l]->rings.radii.size(), 0.0);
  double yTopBar = 0;
  CylindricalResult& cr = out.result.cylindrical;
  cr.present = true;
  cr.planeWave = pw;
  cr.nodes = static_cast<int>(n);
  const bool rotation = P.stack.cylindrical.rotation > 0;
  if (rotation) {
    cr.rotation = P.stack.cylindrical.rotation;
    cr.classL = P.classes.front().L;
  }
  double value = 0;

  // one unit per harmonic (M20) or per symmetry class (M21): state index b N_k + j, n_f = (2L + 1) N_k
  for (std::size_t u = 0; u < P.classes.size(); ++u) {
    const cyl::HarmonicMap& map = P.classes[u];
    const int nb = map.blocks(), m = map.m0;
    const Eigen::Index nf = nb * n;
    FG_REQUIRE(!pw || nb == 1, "gradient (cylindrical): plane waves per harmonic only");
    const VecC zero = VecC::Zero(2 * nf);
    const VecC aIn = pw ? zero : P.aIn[u];
    if (!pw && aIn.isZero(0.0)) {
      if (!rotation) {  // M20 lists every harmonic of the incidence (as the solver)
        cr.harmonics.push_back(m);
        cr.reflectedPerHarmonic.push_back(0);
        cr.transmittedPerHarmonic.push_back(0);
      }
      continue;
    }
    const VecC kt = cyl::tiledK(P.grid.k, nb);
    // ---- forward: the chain of the solver (substrate -> reference, then the layers bottom to top)
    CylBasis basis(P.grid, map, P.epsSup, P.epsSub, P.layers, P.fac);
    if (pw) basis.setPlaneWave(P.kpar, P.Psi);
    const ModeBasis& ref = basis.ref();
    const std::size_t L = P.layers.size();
    std::vector<SourcedSMatrix> E(L), pre(L + 1);
    std::vector<LayerModes> modes(L);
    std::vector<cyl::CylQ> Q(L);
    std::vector<MatC> Ml(L);
    // plane wave: per ring layer the particular solution p_k = (kappa_k - M)^-1 sigma_k a_k (as modalSources)
    struct Src {
      cyl::PlaneWaveSource B;
      std::unique_ptr<Eigen::PartialPivLU<MatC>> Qyy;
      std::vector<VecC> p, sigma;
      std::vector<BgValues> bv;
      VecC g0, g1;
    };
    std::vector<Src> src(L);
    pre[0] = {fromTransfer(ref.PsiInv * basis.sub().Psi), zero, zero};
    for (std::size_t l = 0; l < L; ++l) {
      const CylLayer& lay = *P.layers[l];
      const double H = lay.thickness_um() * k0;
      LayerModes md;
      if (lay.homogeneous()) {
        md = basis.homogeneousModes(lay);
        const MatC Id = MatC::Identity(nf, nf);
        Q[l] = {lay.eps * Id, MatC::Zero(nf, nf), MatC::Zero(nf, nf), lay.eps * Id, lay.eps * Id, MatC()};
      } else {
        Q[l] = basis.Q(l);  // cylQ (block-diagonal over a class) or coupledQ (angular layers)
        Ml[l] = cyl::propagationMatrix(kt, Q[l]);
        md = basis.modes(Ml[l]);
      }
      RCWASliceRecord rec;
      E[l].S = stableLayerSMatrix(md, H, ref, &rec);
      modes[l] = std::move(rec.modes);
      E[l].u = zero;
      E[l].v = zero;
      if (pw && !lay.homogeneous()) {
        Src& sr = src[l];
        sr.B = cyl::planeWaveSource(P.grid, m, lay.rings, P.kpar, P.Psi, P.fac, Q[l]);
        sr.Qyy = std::make_unique<Eigen::PartialPivLU<MatC>>(Q[l].Qyy);
        const double lmax = std::max(modes[l].lu.cwiseAbs().maxCoeff(), modes[l].ld.cwiseAbs().maxCoeff());
        sr.bv = bgValues(P, l);
        VecC F0 = VecC::Zero(4 * n), F1 = VecC::Zero(4 * n);
        for (const BgValues& t : sr.bv) {
          double dist = 1e300;
          for (Eigen::Index q = 0; q < modes[l].lu.size(); ++q)
            dist = std::min({dist, std::abs(t.kappa - modes[l].lu(q)), std::abs(t.kappa - modes[l].ld(q))});
          FG_REQUIRE(dist >= 1e-8 * lmax, "gradient (cylindrical plane wave): a background exponent resonates with "
                                          "a slice mode (the solver integrates that slice with RK)");
          const VecC sig = cyl::scatteredSource(P.grid.k, *sr.Qyy, sr.B, t.cp, t.cm, t.ey);  // linear in (c, e_y)
          const VecC pk = (t.kappa * MatC::Identity(4 * n, 4 * n) - Ml[l]).partialPivLu().solve(sig);
          sr.sigma.push_back(sig);
          sr.p.push_back(pk);
          F0 += pk * t.e0;
          F1 += pk * t.e1;
        }
        sr.g0 = ref.PsiInv * F0;
        sr.g1 = ref.PsiInv * F1;
        const SMatrix& S = E[l].S;
        E[l].u = sr.g1.tail(2 * n) - S.S11 * sr.g0.tail(2 * n) - S.S12 * sr.g1.head(2 * n);
        E[l].v = sr.g0.head(2 * n) - S.S21 * sr.g0.tail(2 * n) - S.S22 * sr.g1.head(2 * n);
      }
      pre[l + 1] = star(pre[l], E[l]);
    }
    const SourcedSMatrix& T = pre[L];
    const VecC b = T.S.S12 * aIn + T.u, aSub = T.S.S22 * aIn + T.v;
    VecC V = VecC::Zero(4 * nf);
    V.tail(2 * nf) = b;
    const VecC Fup = ref.Psi * V;
    V.setZero();
    V.head(2 * nf) = aSub;
    const VecC Fdown = basis.sub().Psi * V;
    // powers per harmonic (orthogonal in the half spaces) and their cotangents
    VecC FupBar = VecC::Zero(4 * nf), FdownBar = VecC::Zero(4 * nf);
    for (int bk = 0; bk < nb; ++bk) {
      const int mb = map.m(bk);
      const VecC fu = blockState(Fup, bk, n), fd = blockState(Fdown, bk, n);
      const double up = cyl::flux(P.grid, fu), down = -cyl::flux(P.grid, fd);
      cr.harmonics.push_back(mb);
      cr.reflectedPerHarmonic.push_back(up * norm);
      cr.transmittedPerHarmonic.push_back(down * norm);
      const double wR = weightOf('R', mb), wT = weightOf('T', mb);
      value += wR * up + wT * down;
      addBlockState(FupBar, cyl::fluxBar(P.grid, fu, wR), bk, n);
      addBlockState(FdownBar, cyl::fluxBar(P.grid, fd, -wT), bk, n);
    }

    // ---- reverse
    const VecC bBar = (ref.Psi.adjoint() * FupBar).tail(2 * nf);
    const VecC aSubBar = (basis.sub().Psi.adjoint() * FdownBar).head(2 * nf);
    SourcedSMatrix cur;
    cur.S.S11 = MatC::Zero(2 * nf, 2 * nf);
    cur.S.S21 = MatC::Zero(2 * nf, 2 * nf);
    cur.S.S12 = bBar * aIn.adjoint();
    cur.S.S22 = aSubBar * aIn.adjoint();
    cur.u = bBar;
    cur.v = aSubBar;
    if (!pw) {  // the incident beam: a_j ~ exp(i beta_j (y_f - y_top)) (Plan sec:cylexc)
      const VecC aInBar = T.S.S12.adjoint() * bBar + T.S.S22.adjoint() * aSubBar;
      for (Eigen::Index q = 0; q < 2 * nf; ++q) {
        const cplx kq = P.grid.k(q % n), beta = std::sqrt(P.epsSup - kq * kq);
        yTopBar += std::real(std::conj(aInBar(q)) * (-I * beta * aIn(q)));
      }
    }
    for (std::size_t l = L; l-- > 0;) {
      SourcedSMatrix Ab, Eb;
      adjoint::star(pre[l].S, E[l].S, cur.S, Ab.S, Eb.S);
      if (pw) adjoint::starSources(pre[l].S, E[l].S, pre[l].u, E[l].v, cur.u, cur.v, Ab.S, Eb.S, Ab.u, Ab.v, Eb.u, Eb.v);
      cur = std::move(Ab);
      if (!pw) cur.u = cur.v = zero;
      const CylLayer& lay = *P.layers[l];
      const double H = lay.thickness_um() * k0;
      // plane wave: the source amplitudes u, v of the ring layer (cotangents Eb.u, Eb.v)
      MatC Msrc, Bpmbar;
      VecC Bybar;
      MatC Qyysrc;
      if (pw && !lay.homogeneous()) {
        const Src& sr = src[l];
        const SMatrix& S = E[l].S;
        const VecC& ub = Eb.u;
        const VecC& vb = Eb.v;
        const VecC g0t = sr.g0.tail(2 * n), g1h = sr.g1.head(2 * n);
        Eb.S.S11 -= ub * g0t.adjoint();
        Eb.S.S12 -= ub * g1h.adjoint();
        Eb.S.S21 -= vb * g0t.adjoint();
        Eb.S.S22 -= vb * g1h.adjoint();
        VecC g0b(4 * n), g1b(4 * n);
        g0b << vb, -(S.S11.adjoint() * ub + S.S21.adjoint() * vb);
        g1b << -(S.S12.adjoint() * ub + S.S22.adjoint() * vb), ub;
        const VecC F0b = ref.PsiInv.adjoint() * g0b, F1b = ref.PsiInv.adjoint() * g1b;
        Msrc = MatC::Zero(4 * n, 4 * n);
        Bpmbar = MatC::Zero(2 * n, 2);
        Bybar = VecC::Zero(n);
        Qyysrc = MatC::Zero(n, n);
        for (std::size_t t = 0; t < sr.p.size(); ++t) {
          const BgValues& term = sr.bv[t];
          BgValues& tb = bgBar[l][t];
          tb.e0 += sr.p[t].dot(F0b);  // F = sum p e: ebar = p^H Fbar
          tb.e1 += sr.p[t].dot(F1b);
          const VecC pb = std::conj(term.e0) * F0b + std::conj(term.e1) * F1b;
          // p = A^-1 sigma, A = kappa - M: rhs bar = A^-H pbar, Mbar += rhsbar p^H, kappabar = -p^H rhsbar
          const MatC A = term.kappa * MatC::Identity(4 * n, 4 * n) - Ml[l];
          const VecC rb = A.adjoint().partialPivLu().solve(pb);
          Msrc += rb * sr.p[t].adjoint();
          tb.kappa -= sr.p[t].dot(rb);
          const VecC& sb = rb;
          // sigma = (k o x, -k o x, Bpm (c+, c-) head, -tail), x = Q_yy^-1 (e_y B_y)
          const VecC x = sr.Qyy->solve(VecC(term.ey * sr.B.By));
          const VecC xb = P.grid.k.conjugate().cwiseProduct(VecC(sb.head(n) - sb.segment(n, n)));
          const VecC z = Q[l].Qyy.adjoint().partialPivLu().solve(xb);
          Bybar += std::conj(term.ey) * z;
          Qyysrc -= z * x.adjoint();
          tb.ey += VecC(sr.Qyy->solve(sr.B.By)).dot(xb);  // x = e_y Q_yy^-1 B_y
          VecC ppmb(2 * n);
          ppmb << sb.segment(2 * n, n), -sb.segment(3 * n, n);
          Bpmbar.col(0) += ppmb * std::conj(term.cp);
          Bpmbar.col(1) += ppmb * std::conj(term.cm);
          tb.cp += VecC(sr.B.Bpm.col(0)).dot(ppmb);
          tb.cm += VecC(sr.B.Bpm.col(1)).dot(ppmb);
        }
      }
      MatC Wub, Wdb;
      VecC lub, ldb;
      double hb = 0;
      adjoint::layerSMatrix(modes[l], H, ref, Eb.S, Wub, Wdb, lub, ldb, hb);
      hBar[l] += hb;
      const Eigen::Index nm = modes[l].lu.size();
      MatC W(4 * nf, 2 * nm), Wb(4 * nf, 2 * nm);
      W << modes[l].Wu, modes[l].Wd;
      Wb << Wub, Wdb;
      VecC lam(2 * nm), lamb(2 * nm);
      lam << modes[l].lu, modes[l].ld;
      lamb << lub, ldb;
      std::vector<int> cls(static_cast<std::size_t>(2 * nm), 0);
      for (Eigen::Index q = nm; q < 2 * nm; ++q) cls[static_cast<std::size_t>(q)] = 1;
      MatC Mb = adjoint::eig(W, lam, Wb, lamb, cls);
      if (Msrc.size()) Mb += Msrc;
      cyl::CylQ Qb = cyl::propagationMatrixBar(kt, Q[l], Mb);
      if (lay.homogeneous()) {
        epsBar[names[l][0]] += Qb.Qpp.trace() + Qb.Qmm.trace() + Qb.Qyy.trace();
        continue;
      }
      if (lay.angular) {  // M21: coupledQ (Laurent), permittivities and the shape of the regions
        const cyl::CoupledQBar cb = cyl::coupledQBar(P.grid, map, lay.section, Qb);
        for (std::size_t i = 0; i < cb.eps.size(); ++i) epsBar[names[l][i]] += cb.eps[i];
        epsBar[names[l].back()] += cb.epsBg;
        for (std::size_t q = 0; q < ps.size(); ++q)
          if (ps[q].kind == CylParam::Kind::Shape && ps[q].layer == l)
            out.grad[q] += cyl::coupledQShape(cb, map, lay.section, ps[q].velocity);
        continue;
      }
      if (nb > 1) {  // ring layer of a class: block-diagonal, cylQBar per harmonic
        for (int bk = 0; bk < nb; ++bk) {
          const cyl::CylQBar qb = cyl::cylQBar(P.grid, map.m(bk), lay.rings, P.fac, diagonalBlock(Q[l], bk, n),
                                               diagonalBlock(Qb, bk, n));
          for (std::size_t i = 0; i < qb.eps.size(); ++i) {
            epsBar[names[l][i]] += qb.eps[i];
            rBar[l][i] += qb.radii[i];
          }
          epsBar[names[l].back()] += qb.epsInf;
        }
        continue;
      }
      cyl::CylQBar qs;
      if (Msrc.size()) {
        Qb.Qyy += Qyysrc;
        cyl::planeWaveSourceBar(P.grid, m, lay.rings, P.kpar, P.Psi, P.fac, Q[l], src[l].B, Bpmbar, Bybar, Qb, qs);
      }
      const cyl::CylQBar qb = cyl::cylQBar(P.grid, m, lay.rings, P.fac, Q[l], Qb);
      for (std::size_t i = 0; i < qb.eps.size(); ++i) {
        epsBar[names[l][i]] += qb.eps[i] + (qs.eps.empty() ? cplx(0.0) : qs.eps[i]);
        rBar[l][i] += qb.radii[i] + (qs.radii.empty() ? 0.0 : qs.radii[i]);
      }
      epsBar[names[l].back()] += qb.epsInf + qs.epsInf;
    }
  }
  out.value = value;
  double sumR = 0, sumT = 0;
  for (double v : cr.reflectedPerHarmonic) sumR += v;
  for (double v : cr.transmittedPerHarmonic) sumT += v;
  if (pw) {  // cross sections in um^2 (volume integrals are not part of the gradient result)
    cr.scatteredUp_um2 = sumR;
    cr.scatteredDown_um2 = sumT;
    cr.scattering_um2 = sumR + sumT;
  } else {
    cr.reflected = sumR;
    cr.transmitted = sumT;
    cr.absorbed = 1.0 - sumR - sumT;
    out.result.absorption = cr.absorbed;
    out.result.diag.energyBalance = sumR + sumT;
  }

  // Stufe 2b: d(background quantities)/d theta by central differences of the background problem (as M19 B2)
  for (std::size_t q = 0; q < ps.size(); ++q) {
    if (!moveBg[q]) continue;
    const json::json_pointer jp(params[q]);
    const double v = j.at(jp).get<double>(), hstep = 1e-6 * std::max(1.0, std::abs(v));
    auto at = [&](double sgn) {  // two variables, no container (CylProblem is move-only)
      json a = j;
      a[jp] = v + sgn * hstep;
      const io::RunConfig ca = io::parseConfig(a.dump());
      return cylProblem(ca.solver, ca.stack, wl);
    };
    const CylProblem Pplus = at(1.0), Pminus = at(-1.0);
    double dL = 0;
    for (std::size_t l = 0; l < P.layers.size(); ++l) {
      if (bgBar[l].empty()) continue;
      const auto a = bgValues(Pplus, l), b = bgValues(Pminus, l);
      FG_REQUIRE(a.size() == bgBar[l].size() && b.size() == bgBar[l].size(),
                 "gradient (cylindrical plane wave): the background exponentials change with '" + params[q] + "'");
      for (std::size_t t = 0; t < a.size(); ++t) {
        const BgValues& tb = bgBar[l][t];
        auto d = [&](cplx x, cplx y) { return (x - y) / (2 * hstep); };
        dL += std::real(std::conj(tb.kappa) * d(a[t].kappa, b[t].kappa) + std::conj(tb.cp) * d(a[t].cp, b[t].cp) +
                        std::conj(tb.cm) * d(a[t].cm, b[t].cm) + std::conj(tb.ey) * d(a[t].ey, b[t].ey) +
                        std::conj(tb.e0) * d(a[t].e0, b[t].e0) + std::conj(tb.e1) * d(a[t].e1, b[t].e1));
      }
    }
    out.grad[q] += dL;
  }
  for (std::size_t q = 0; q < ps.size(); ++q) {
    const CylParam& p = ps[q];
    switch (p.kind) {
      case CylParam::Kind::Thickness:  // the layer and the top of the stack (phase of the beam) move
        out.grad[q] += (hBar[p.layer] + yTopBar) * k0;
        break;
      case CylParam::Kind::Radius:
        out.grad[q] += rBar[p.layer][p.index] * k0;
        break;
      case CylParam::Kind::Shape:  // added per class above (normalized velocities per unit of theta)
        break;
      case CylParam::Kind::Material: {
        const cplx eb = epsBar.count(p.material) ? epsBar[p.material] : cplx(0.0);
        cplx deps = p.comp == 0 ? cplx(1.0) : I;
        if (p.viaN) deps *= 2.0 * std::sqrt(cfg.stack.material(p.material).eps);  // eps = n^2
        out.grad[q] += std::real(std::conj(eb) * deps);
        break;
      }
    }
  }
  out.notes = {"cylindrical: spectral grid (k_max, panels, deformed path, graded panels) fixed",
               pw ? "cylindrical: plane wave (scattered field); the background solution of the incidence order by "
                    "central differences for thicknesses and background materials (V21i Stufe 2b)"
                  : "cylindrical: Gaussian beam, total field, RCWA (V21i)"};
  return out;
}

}  // namespace fg::detail
