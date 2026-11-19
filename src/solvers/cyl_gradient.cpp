// FourierGrating - gradients of cylindrical structures (M20, Plan sections "Differenzierbarkeit" and V21i): the
// discrete adjoint of the beam computation (total field, RCWA). Per harmonic: powers <- outgoing amplitudes <-
// total S-matrix <- star products <- layer S-matrices <- modes <- M <- Q <- Galerkin matrices <- ring
// permittivities and radii (Plan eq:cyl_shape); thicknesses through the layer S-matrices and the phase of the
// incident beam (referenced at the top of the stack).
#include "cyl_gradient.hpp"

#include <map>
#include <nlohmann/json.hpp>
#include <sstream>
#include <variant>

#include "cyl_basis.hpp"
#include "fg/adjoint/primitives.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::detail {

namespace {

using json = nlohmann::json;

struct CylParam {
  enum class Kind { Thickness, Radius, Material };
  Kind kind = Kind::Thickness;
  std::size_t layer = 0;  ///< expanded stack
  std::size_t index = 0;  ///< radius index
  std::string material;
  bool viaN = false;
  int comp = 0;  ///< 0: real part, 1: imaginary part
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
  } else {
    throw Error("gradient (cylindrical): unsupported layer parameter '" + ptr + "'");
  }
  return p;
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
  FG_REQUIRE(!P.planeWave, "gradient (cylindrical): plane waves (scattered field) follow; Gaussian beams only");
  FG_REQUIRE(!P.beamScattered, "gradient (cylindrical): the total-field formulation only (scattered_field false)");
  const Eigen::Index n = P.grid.k.size();
  const double k0 = P.k0;
  // material names of the (expanded) layers: homogeneous -> one, rings -> inner to outer
  std::vector<std::vector<std::string>> names;
  for (const LayerEntry& le : P.stack.layers) {
    if (const auto* h = std::get_if<HomogeneousLayer>(&le.spec)) names.push_back({h->material});
    else names.push_back(std::get<RingsLayer>(le.spec).materials);
  }
  for (const CylParam& p : ps)
    if (p.kind == CylParam::Kind::Radius)
      FG_REQUIRE(p.index < P.layers[p.layer]->rings.radii.size(), "gradient: radius index out of range");

  // incident power (independent of the parameters: the phase of the beam does not change it)
  double pInc = 0;
  for (std::size_t h = 0; h < P.ms.size(); ++h) {
    if (P.aIn[h].isZero(0.0)) continue;
    const ModeBasis ref = cyl::homogeneousModes(P.epsSup, P.grid.k);
    VecC V = VecC::Zero(4 * n);
    V.head(2 * n) = P.aIn[h];
    pInc -= cyl::flux(P.grid, ref.Psi * V);
  }
  FG_REQUIRE(pInc > 0, "cylindrical: the incident beam carries no power (harmonics?)");
  // weights of the raw powers per harmonic
  auto weightOf = [&](char side, int m) {
    double w = 0;
    for (const auto& [r, v] : objective.cylindrical)
      if (r.side == side && (r.total || r.m == m)) w += v;
    return w / pInc;
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
  cr.harmonics = P.ms;
  cr.nodes = static_cast<int>(n);
  double value = 0;

  for (std::size_t h = 0; h < P.ms.size(); ++h) {
    const int m = P.ms[h];
    const VecC& aIn = P.aIn[h];
    if (aIn.isZero(0.0)) {
      cr.reflectedPerHarmonic.push_back(0);
      cr.transmittedPerHarmonic.push_back(0);
      continue;
    }
    // ---- forward: the chain of the solver (substrate -> reference, then the layers bottom to top)
    const CylBasis basis(P.grid, m, P.epsSup, P.epsSub, P.layers, P.fac);
    const ModeBasis& ref = basis.ref();
    const std::size_t L = P.layers.size();
    std::vector<SMatrix> S(L), prefix(L + 1);
    std::vector<LayerModes> modes(L);
    std::vector<cyl::CylQ> Q(L);
    prefix[0] = fromTransfer(ref.PsiInv * basis.sub().Psi);
    for (std::size_t l = 0; l < L; ++l) {
      const CylLayer& lay = *P.layers[l];
      LayerModes md;
      if (lay.homogeneous()) {
        md = basis.homogeneousModes(lay);
        const MatC Id = MatC::Identity(n, n);
        Q[l] = {lay.eps * Id, MatC::Zero(n, n), MatC::Zero(n, n), lay.eps * Id, lay.eps * Id, MatC()};
      } else {
        Q[l] = cyl::cylQ(P.grid, m, lay.rings, P.fac);
        md = basis.modes(cyl::propagationMatrix(P.grid.k, Q[l]));
      }
      RCWASliceRecord rec;
      S[l] = stableLayerSMatrix(md, lay.thickness_um() * k0, ref, &rec);
      modes[l] = std::move(rec.modes);
      prefix[l + 1] = star(prefix[l], S[l]);
    }
    const SMatrix& T = prefix[L];
    const VecC b = T.S12 * aIn, aSub = T.S22 * aIn;
    VecC V = VecC::Zero(4 * n);
    V.tail(2 * n) = b;
    const VecC Fup = ref.Psi * V;
    V.setZero();
    V.head(2 * n) = aSub;
    const VecC Fdown = basis.sub().Psi * V;
    const double up = cyl::flux(P.grid, Fup), down = -cyl::flux(P.grid, Fdown);
    cr.reflectedPerHarmonic.push_back(up / pInc);
    cr.transmittedPerHarmonic.push_back(down / pInc);
    const double wR = weightOf('R', m), wT = weightOf('T', m);
    value += wR * up + wT * down;

    // ---- reverse
    const VecC bBar = (ref.Psi.adjoint() * cyl::fluxBar(P.grid, Fup, wR)).tail(2 * n);
    const VecC aSubBar = (basis.sub().Psi.adjoint() * cyl::fluxBar(P.grid, Fdown, -wT)).head(2 * n);
    SMatrix Tb;
    Tb.S11 = MatC::Zero(2 * n, 2 * n);
    Tb.S21 = MatC::Zero(2 * n, 2 * n);
    Tb.S12 = bBar * aIn.adjoint();
    Tb.S22 = aSubBar * aIn.adjoint();
    // the incident beam: a_j ~ exp(i beta_j (y_f - y_top)) (Plan sec:cylexc)
    const VecC aInBar = T.S12.adjoint() * bBar + T.S22.adjoint() * aSubBar;
    for (Eigen::Index j = 0; j < 2 * n; ++j) {
      const cplx kj = P.grid.k(j % n), beta = std::sqrt(P.epsSup - kj * kj);
      yTopBar += std::real(std::conj(aInBar(j)) * (-I * beta * aIn(j)));
    }
    SMatrix cur = Tb;
    for (std::size_t l = L; l-- > 0;) {
      SMatrix Pb, Sb;
      adjoint::star(prefix[l], S[l], cur, Pb, Sb);
      cur = std::move(Pb);
      const CylLayer& lay = *P.layers[l];
      MatC Wub, Wdb;
      VecC lub, ldb;
      double hb = 0;
      adjoint::layerSMatrix(modes[l], lay.thickness_um() * k0, ref, Sb, Wub, Wdb, lub, ldb, hb);
      hBar[l] += hb;
      const Eigen::Index nm = modes[l].lu.size();
      MatC W(4 * n, 2 * nm), Wb(4 * n, 2 * nm);
      W << modes[l].Wu, modes[l].Wd;
      Wb << Wub, Wdb;
      VecC lam(2 * nm), lamb(2 * nm);
      lam << modes[l].lu, modes[l].ld;
      lamb << lub, ldb;
      std::vector<int> cls(static_cast<std::size_t>(2 * nm), 0);
      for (Eigen::Index q = nm; q < 2 * nm; ++q) cls[static_cast<std::size_t>(q)] = 1;
      const MatC Mb = adjoint::eig(W, lam, Wb, lamb, cls);
      const cyl::CylQ Qb = cyl::propagationMatrixBar(P.grid.k, Q[l], Mb);
      if (lay.homogeneous()) {
        epsBar[names[l][0]] += Qb.Qpp.trace() + Qb.Qmm.trace() + Qb.Qyy.trace();
        continue;
      }
      const cyl::CylQBar qb = cyl::cylQBar(P.grid, m, lay.rings, P.fac, Q[l], Qb);
      for (std::size_t i = 0; i < qb.eps.size(); ++i) {
        epsBar[names[l][i]] += qb.eps[i];
        rBar[l][i] += qb.radii[i];
      }
      epsBar[names[l].back()] += qb.epsInf;
    }
  }
  out.value = value;
  cr.reflected = cr.transmitted = 0;
  for (double v : cr.reflectedPerHarmonic) cr.reflected += v;
  for (double v : cr.transmittedPerHarmonic) cr.transmitted += v;
  cr.absorbed = 1.0 - cr.reflected - cr.transmitted;
  out.result.absorption = cr.absorbed;
  out.result.diag.energyBalance = cr.reflected + cr.transmitted;

  for (std::size_t q = 0; q < ps.size(); ++q) {
    const CylParam& p = ps[q];
    switch (p.kind) {
      case CylParam::Kind::Thickness:  // the layer and the top of the stack (phase of the beam) move
        out.grad[q] = (hBar[p.layer] + yTopBar) * k0;
        break;
      case CylParam::Kind::Radius:
        out.grad[q] = rBar[p.layer][p.index] * k0;
        break;
      case CylParam::Kind::Material: {
        const cplx eb = epsBar.count(p.material) ? epsBar[p.material] : cplx(0.0);
        cplx deps = p.comp == 0 ? cplx(1.0) : I;
        if (p.viaN) deps *= 2.0 * std::sqrt(cfg.stack.material(p.material).eps);  // eps = n^2
        out.grad[q] = std::real(std::conj(eb) * deps);
        break;
      }
    }
  }
  out.notes = {"cylindrical: spectral grid (k_max, panels, deformed path, graded panels) fixed",
               "cylindrical: Gaussian beam, total field, RCWA (V21i)"};
  return out;
}

}  // namespace fg::detail
