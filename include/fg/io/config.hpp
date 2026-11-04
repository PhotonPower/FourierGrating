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
  std::vector<Wavelength> scan;            ///< one or more incidence conditions
  bool nearField = false;
  NearFieldRequest field;
  std::string outputDir = "out";
};

/// Parses a configuration file (units: micrometres, degrees). TODO(M2): full schema.
RunConfig loadConfig(const std::string& path);

}  // namespace fg::io
