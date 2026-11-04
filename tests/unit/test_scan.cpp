// Parameter scans and sweeps (M11, docs/config_format.md).
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "fg/core/error.hpp"
#include "fg/io/config.hpp"

namespace {
const char* kBase = R"({
  "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 6},
  "period_um": {"x": 0.5},
  "incidence": {"lambda_um": {"from": 0.6, "to": 0.7, "steps": 3}, "phi_deg": [0, 10],
                "psi_deg": 0, "delta_deg": [0, 90]},
  "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
  "superstrate": "air", "substrate": "glass",
  "layers": [{"type": "lamellar1d", "thickness_um": 0.2, "breakpoints_um": [0.0, 0.2, 0.5],
              "materials": ["glass", "air"]}]
  SWEEP
})";

std::string config(const std::string& sweep) {
  std::string s = kBase;
  s.replace(s.find("SWEEP"), 5, sweep);
  return s;
}
}  // namespace

TEST(Scan, IncidenceIsCartesianProduct) {
  const auto c = fg::io::parseConfig(config(""));
  ASSERT_EQ(c.scan.size(), 3u * 2u * 2u);
  EXPECT_NEAR(c.scan[0].lambda_um, 0.6, 1e-15);
  EXPECT_NEAR(c.scan.back().lambda_um, 0.7, 1e-15);
  EXPECT_NEAR(c.scan[1].delta, fg::PI / 2, 1e-15);  // innermost loop: delta
  EXPECT_NEAR(c.scan[2].phi, 10 * fg::PI / 180, 1e-15);
}

TEST(Scan, SweepOverLayerThickness) {
  const auto runs = fg::io::parseSweep(
      config(R"(, "sweep": {"parameter": "/layers/0/thickness_um", "values": {"from": 0.1, "to": 0.3, "steps": 5}})"));
  ASSERT_EQ(runs.size(), 5u);
  for (std::size_t k = 0; k < runs.size(); ++k) {
    const double t = 0.1 + 0.05 * static_cast<double>(k);
    EXPECT_NEAR(runs[k].sweepValue, t, 1e-15);
    EXPECT_EQ(runs[k].sweepParameter, "/layers/0/thickness_um");
    EXPECT_NEAR(std::get<fg::Lamellar1DLayer>(runs[k].stack.layers[0].spec).thickness_um, t, 1e-15);
    EXPECT_EQ(runs[k].scan.size(), 12u);
  }
  EXPECT_EQ(fg::io::parseSweep(config("")).size(), 1u);
  EXPECT_TRUE(fg::io::parseSweep(config("")).front().sweepParameter.empty());
  EXPECT_THROW(fg::io::parseSweep(config(R"(, "sweep": {"parameter": "/layers/5/x", "values": [1]})")), fg::Error);
  EXPECT_THROW(fg::io::parseSweep(config(R"(, "sweep": {"parameter": "/superstrate", "values": [1]})")), fg::Error);
}

// Swept results equal individually configured runs; CSV has the sweep columns.
TEST(Scan, SweepResultsAndCsv) {
  const auto runs =
      fg::io::parseSweep(config(R"(, "sweep": {"parameter": "/layers/0/thickness_um", "values": [0.15, 0.25]})"));
  std::vector<std::vector<fg::Result>> res;
  for (const auto& r : runs) res.push_back(fg::solveScan(r.solver, r.stack, r.scan));
  std::string single = config("");
  single.replace(single.find("\"thickness_um\": 0.2"), 19, "\"thickness_um\": 0.25");
  const auto c = fg::io::parseConfig(single);
  const auto ref = fg::makeSolver(c.solver)->solve(c.stack, c.scan[5]);
  ASSERT_EQ(ref.reflected.size(), res[1][5].reflected.size());
  for (std::size_t k = 0; k < ref.reflected.size(); ++k)
    EXPECT_EQ(ref.reflected[k].efficiency, res[1][5].reflected[k].efficiency);
  const auto file = (std::filesystem::temp_directory_path() / "fg_sweep.csv").string();
  fg::io::writeEfficiencies(file, runs, res);
  std::ifstream in(file);
  std::string header, row;
  std::getline(in, header);
  std::getline(in, row);
  EXPECT_EQ(header.rfind("sweep_parameter,sweep_value,lambda_um,", 0), 0u);
  EXPECT_EQ(row.rfind("/layers/0/thickness_um,0.15,0.6,", 0), 0u);
}
