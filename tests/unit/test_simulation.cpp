// Simulation (Plan section "M12: API"): kept layer models, set_parameter = rebuild (bit-identical).
#include <gtest/gtest.h>

#include <nlohmann/json.hpp>
#include <string>

#include "fg/core/error.hpp"
#include "fg/solvers/simulation.hpp"

namespace {

const char* kConfig = R"({
  "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 8},
  "period_um": {"x": 0.55},
  "incidence": {"lambda_um": 0.63, "phi_deg": 17, "psi_deg": 10, "delta_deg": 35},
  "rcwa": {"slices": 10},
  "materials": {"air": {"n": 1.0}, "glass": {"n": [1.5, 0.01]}, "al": {"n": [1.2, 7.0]}},
  "superstrate": "air", "substrate": "glass",
  "layers": [
    {"type": "homogeneous", "thickness_um": 0.07, "material": "glass"},
    {"type": "profile1d", "shape": "sinus", "depth_um": 0.12, "below": "al", "above": "air"},
    {"type": "lamellar1d", "thickness_um": 0.15, "breakpoints_um": [0.0, 0.2, 0.55], "materials": ["glass", "air"]}
  ]
})";

/// Bit-identical comparison of two results.
void expectIdentical(const fg::Result& a, const fg::Result& b) {
  ASSERT_EQ(a.reflected.size(), b.reflected.size());
  ASSERT_EQ(a.transmitted.size(), b.transmitted.size());
  for (std::size_t k = 0; k < a.reflected.size(); ++k) {
    EXPECT_EQ(a.reflected[k].efficiency, b.reflected[k].efficiency);
    EXPECT_EQ(a.reflected[k].amp_s, b.reflected[k].amp_s);
    EXPECT_EQ(a.reflected[k].amp_p, b.reflected[k].amp_p);
  }
  for (std::size_t k = 0; k < a.transmitted.size(); ++k)
    EXPECT_EQ(a.transmitted[k].efficiency, b.transmitted[k].efficiency);
  EXPECT_EQ(a.absorption, b.absorption);
}

std::string changed(const std::string& pointer, double v) {
  auto j = nlohmann::json::parse(kConfig);
  j[nlohmann::json::json_pointer(pointer)] = v;
  return j.dump();
}

}  // namespace

TEST(Simulation, EqualsSolver) {
  fg::Simulation sim(kConfig);
  const auto& c = sim.config();
  expectIdentical(sim.solve(), fg::makeSolver(c.solver)->solve(c.stack, c.scan.front()));
  // repeated solves with the kept models, other incidence
  fg::Wavelength w = c.scan.front();
  w.lambda_um = 0.7;
  expectIdentical(sim.solve(w), fg::makeSolver(c.solver)->solve(c.stack, w));
  EXPECT_EQ(sim.modelBuilds(), 3);
}

TEST(Simulation, SetParameterRebuildsOnlyTheLayer) {
  fg::Simulation sim(kConfig);
  sim.setParameter("/layers/1/depth_um", 0.09);
  EXPECT_EQ(sim.modelBuilds(), 4);
  expectIdentical(sim.solve(), fg::Simulation(changed("/layers/1/depth_um", 0.09)).solve());
  // incidence: no rebuild
  sim.setParameter("/incidence/lambda_um", 0.66);
  EXPECT_EQ(sim.modelBuilds(), 4);
  auto j = nlohmann::json::parse(changed("/layers/1/depth_um", 0.09));
  j["incidence"]["lambda_um"] = 0.66;
  expectIdentical(sim.solve(), fg::Simulation(j.dump()).solve());
}

TEST(Simulation, SetParameterOutsideLayersRebuildsAll) {
  fg::Simulation sim(kConfig);
  sim.setParameter("/period_um/x", 0.6);
  EXPECT_EQ(sim.modelBuilds(), 6);
  expectIdentical(sim.solve(), fg::Simulation(changed("/period_um/x", 0.6)).solve());
  sim.setParameter("/materials/al/n/0", 1.3);
  EXPECT_EQ(sim.modelBuilds(), 9);
}

TEST(Simulation, SetParameterErrors) {
  fg::Simulation sim(kConfig);
  EXPECT_THROW(sim.setParameter("/layers/7/thickness_um", 1.0), fg::Error);
  EXPECT_THROW(sim.setParameter("/superstrate", 1.0), fg::Error);
  expectIdentical(sim.solve(), fg::Simulation(kConfig).solve());
}
