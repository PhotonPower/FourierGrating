// FourierGrating - JSON configuration (see docs/config_format.md).
#pragma once

#include <string>
#include <vector>

#include "fg/field/nearfield.hpp"
#include "fg/geometry/stack.hpp"
#include "fg/solvers/solver.hpp"

namespace fg::io {

struct RunConfig {
  Stack stack;
  SolverOptions solver;
  std::vector<Wavelength> scan;            ///< one or more incidence conditions (angles in rad)
  bool nearField = false;
  NearFieldRequest field;
  std::string outputDir = "out";
  std::string efficiencyFile = "efficiencies.csv";
};

/// Parses a configuration file (units: micrometres, degrees).
RunConfig loadConfig(const std::string& path);
/// Parses a configuration from a JSON string (used by tests).
RunConfig parseConfig(const std::string& json);

/// Writes efficiencies of a scan as CSV (docs/config_format.md, "Ausgabe").
void writeEfficiencies(const std::string& file, const std::vector<Wavelength>& scan,
                       const std::vector<Result>& results);

}  // namespace fg::io
