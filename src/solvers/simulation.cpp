// FourierGrating - Simulation: configuration and layer models kept between solves (Plan section
// "M12: API").
#include "fg/solvers/simulation.hpp"

#include <nlohmann/json.hpp>

#include "fg/core/error.hpp"
#include "layer_model.hpp"

namespace fg {

using Json = nlohmann::json;

struct Simulation::Impl {
  Json j;
  io::RunConfig cfg;
  detail::LayerModels models;
  long builds = 0;

  IndexMap map() const { return IndexMap(cfg.solver.Nx, cfg.stack.is2D() ? cfg.solver.Nz : 0); }

  void buildAll() {
    models.clear();
    if (cfg.stack.cylindrical.enabled) return;  // M20: the cylindrical path builds its layers per solve
    detail::applyThreadSettings(cfg.solver);  // e.g. FFT sampling of the normal fields
    models = detail::buildLayerModels(cfg.stack, map(), cfg.solver.rcwaSlices);
    builds += static_cast<long>(models.size());
  }

  void build(std::size_t k) {
    if (cfg.stack.cylindrical.enabled) return;
    detail::applyThreadSettings(cfg.solver);
    models[k] = detail::buildLayerModel(cfg.stack, k, map(), cfg.solver.rcwaSlices);
    ++builds;
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
  if (impl_->cfg.stack.cylindrical.enabled) return makeSolver(impl_->cfg.solver)->solve(impl_->cfg.stack, wl);
  return detail::solveWithModels(impl_->cfg.solver, impl_->cfg.stack, wl, impl_->models);
}

OrderSMatrix Simulation::orderSMatrix(const Wavelength& wl) {
  return detail::orderSMatrixWithModels(impl_->cfg.solver, impl_->cfg.stack, wl, impl_->models);
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
  FG_REQUIRE(!impl_->cfg.stack.cylindrical.enabled, "gradients of cylindrical structures follow in M20 (V21i)");
  return detail::gradientWithModels(impl_->j.dump(), impl_->cfg, wl, &impl_->models, objective, params);
}
std::vector<double> Simulation::vjp(const Objective& cotangents, const std::vector<std::string>& params) {
  return valueAndGrad(cotangents, params).grad;
}

}  // namespace fg
