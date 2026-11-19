// FourierGrating - Simulation: configuration and layer models kept between solves (Plan section
// "M12: API").
#include "fg/solvers/simulation.hpp"

#include <nlohmann/json.hpp>
#include <optional>
#include <type_traits>
#include <variant>

#include "fg/core/error.hpp"
#include "layer_model.hpp"

namespace fg {

using Json = nlohmann::json;

namespace {
/// Materials of a layer; nullopt if they are only known per section height (sectioned/CAD layers).
std::optional<std::vector<std::string>> materialsOf(const LayerSpec& spec) {
  return std::visit(
      [](const auto& l) -> std::optional<std::vector<std::string>> {
        using T = std::decay_t<decltype(l)>;
        if constexpr (std::is_same_v<T, HomogeneousLayer>) return std::vector<std::string>{l.material};
        else if constexpr (std::is_same_v<T, Lamellar1DLayer> || std::is_same_v<T, RingsLayer> ||
                           std::is_same_v<T, RingsProfileLayer>)
          return l.materials;
        else if constexpr (std::is_same_v<T, Profile1DLayer> || std::is_same_v<T, Profile2DLayer>)
          return std::vector<std::string>{l.below, l.above};
        else if constexpr (std::is_same_v<T, Polygon2DLayer>) {
          std::vector<std::string> v{l.background};
          for (const auto& i : l.inclusions) v.push_back(i.material);
          return v;
        } else if constexpr (std::is_same_v<T, Contour2DLayer>) {
          std::vector<std::string> v{l.background};
          for (const auto& r : l.regions) v.push_back(r.material);
          return v;
        } else if constexpr (std::is_same_v<T, Sectioned1DLayer>) {
          if (l.materials.empty()) return std::nullopt;
          return l.materials;
        } else {
          return std::nullopt;
        }
      },
      spec);
}
}  // namespace

struct Simulation::Impl {
  Json j;
  io::RunConfig cfg;
  detail::LayerModels models;
  long builds = 0;
  // M22 (Plan sec:disp_arch): with dispersive materials the models belong to one wavelength
  bool dispersive = false;
  double lambda = 0;  ///< wavelength of the models and of `resolved`
  Stack resolved;     ///< cfg.stack at lambda

  IndexMap map() const { return IndexMap(cfg.solver.Nx, cfg.stack.is2D() ? cfg.solver.Nz : 0); }
  const Stack& stack() const { return dispersive ? resolved : cfg.stack; }

  void buildAll() {
    models.clear();
    // the cylindrical path (M20) builds its layers per solve and resolves the materials itself
    dispersive = cfg.stack.dispersive() && !cfg.stack.cylindrical.enabled;
    if (cfg.stack.cylindrical.enabled) return;
    if (dispersive) {
      FG_REQUIRE(!cfg.scan.empty(), "configuration has no incidence");
      lambda = cfg.scan.front().lambda_um;
      resolved = cfg.stack.at(lambda);
    }
    detail::applyThreadSettings(cfg.solver);  // e.g. FFT sampling of the normal fields
    models = detail::buildLayerModels(stack(), map(), cfg.solver.rcwaSlices);
    builds += static_cast<long>(models.size());
  }

  void build(std::size_t k) {
    if (cfg.stack.cylindrical.enabled) return;
    detail::applyThreadSettings(cfg.solver);
    models[k] = detail::buildLayerModel(stack(), k, map(), cfg.solver.rcwaSlices);
    ++builds;
  }

  bool layerDispersive(std::size_t k) const {
    const auto mats = materialsOf(cfg.stack.layers[k].spec);
    if (!mats) return true;
    for (const auto& n : *mats)
      if (cfg.stack.material(n).dispersive()) return true;
    return false;
  }

  /// Models at lambda: rebuilds only the layers with a dispersive material (Plan sec:disp_arch).
  void at(double lam) {
    if (!dispersive || lam == lambda) return;
    Stack r = cfg.stack.at(lam);  // throws (range, pole) before anything changes
    resolved = std::move(r);
    lambda = lam;
    for (std::size_t k = 0; k < models.size(); ++k)
      if (layerDispersive(k)) build(k);
  }

  /// Everything that the layer models depend on besides their own layer entry.
  static Json modelContext(Json c) {
    for (const char* k : {"layers", "incidence", "output", "nearfield", "sweep", "reference"})
      c.erase(k);
    return c;
  }
};

Simulation::Simulation(const std::string& text) : impl_(std::make_unique<Impl>()) {
  impl_->j = Json::parse(text);
  impl_->cfg = io::parseConfig(text);
  impl_->buildAll();
}

Simulation::~Simulation() = default;
Simulation::Simulation(Simulation&&) noexcept = default;
Simulation& Simulation::operator=(Simulation&&) noexcept = default;

Result Simulation::solve(const Wavelength& wl) {
  Impl& m = *impl_;
  if (m.cfg.stack.cylindrical.enabled) return makeSolver(m.cfg.solver)->solve(m.cfg.stack, wl);
  m.at(wl.lambda_um);
  Result r = detail::solveWithModels(m.cfg.solver, m.stack(), wl, m.models);
  if (m.dispersive) detail::materialDiagnostics(m.cfg.stack, wl.lambda_um, r.diag);
  return r;
}

OrderSMatrix Simulation::orderSMatrix(const Wavelength& wl) {
  impl_->at(wl.lambda_um);
  return detail::orderSMatrixWithModels(impl_->cfg.solver, impl_->stack(), wl, impl_->models);
}

Result Simulation::solve() {
  FG_REQUIRE(!impl_->cfg.scan.empty(), "configuration has no incidence");
  return solve(impl_->cfg.scan.front());
}

void Simulation::setParameter(const std::string& pointer, double value) {
  const Json::json_pointer ptr(pointer);
  Impl& m = *impl_;
  FG_REQUIRE(m.j.contains(ptr), "set_parameter: '" + pointer + "' not found in the configuration");
  FG_REQUIRE(m.j.at(ptr).is_number(), "set_parameter: '" + pointer + "' is not a number");
  Json jn = m.j;
  jn[ptr] = value;
  io::RunConfig cn = io::parseConfig(jn.dump());  // validates before anything is changed
  const bool sameContext =
      Impl::modelContext(m.j) == Impl::modelContext(jn) &&
      m.j.value("layers", Json::array()).size() == jn.value("layers", Json::array()).size();
  std::vector<std::size_t> changed;
  if (sameContext && jn.contains("layers"))
    for (std::size_t k = 0; k < jn["layers"].size(); ++k)
      if (m.j["layers"][k] != jn["layers"][k]) changed.push_back(k);
  m.j = std::move(jn);
  m.cfg = std::move(cn);
  if (!sameContext) {
    m.buildAll();
    return;
  }
  if (m.dispersive) m.resolved = m.cfg.stack.at(m.lambda);  // same materials, new layer entries
  for (std::size_t k : changed) m.build(k);
}

const io::RunConfig& Simulation::config() const {
  return impl_->cfg;
}
std::string Simulation::json() const {
  return impl_->j.dump();
}
long Simulation::modelBuilds() const {
  return impl_->builds;
}

GradientResult Simulation::valueAndGrad(const Objective& objective, const std::vector<std::string>& params) {
  FG_REQUIRE(!impl_->cfg.scan.empty(), "configuration has no incidence");
  return valueAndGrad(impl_->cfg.scan.front(), objective, params);
}
GradientResult Simulation::valueAndGrad(const Wavelength& wl, const Objective& objective,
                                        const std::vector<std::string>& params) {
  if (impl_->cfg.stack.cylindrical.enabled)  // V21i: own forward chain (no layer models)
    return detail::gradientWithModels(impl_->j.dump(), impl_->cfg, wl, nullptr, objective, params);
  impl_->at(wl.lambda_um);  // M22: models at this wavelength; gradientWithModels resolves the configuration
  return detail::gradientWithModels(impl_->j.dump(), impl_->cfg, wl, &impl_->models, objective, params);
}
std::vector<double> Simulation::vjp(const Objective& cotangents, const std::vector<std::string>& params) {
  return valueAndGrad(cotangents, params).grad;
}

}  // namespace fg
