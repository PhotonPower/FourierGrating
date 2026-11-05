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
  std::string sweepParameter;              ///< JSON pointer of the swept value (empty: no sweep)
  double sweepValue = 0;
};

/// Parses a configuration file (units: micrometres, degrees).
RunConfig loadConfig(const std::string& path);
/// Parses a configuration from a JSON string; relative file names (cad layers) against baseDir.
RunConfig parseConfig(const std::string& json, const std::string& baseDir = "");

/// Parameter sweep (M11): "sweep": {"parameter": "<JSON pointer>", "values": <number | array | range>}.
/// Returns one configuration per value (the value written at the pointer); without "sweep" exactly one.
std::vector<RunConfig> parseSweep(const std::string& json, const std::string& baseDir = "");
std::vector<RunConfig> loadSweep(const std::string& path);

/// Writes efficiencies of a scan as CSV (docs/config_format.md, "Ausgabe").
void writeEfficiencies(const std::string& file, const std::vector<Wavelength>& scan,
                       const std::vector<Result>& results);
/// Sweep: one block per configuration, with the leading columns sweep_parameter, sweep_value.
void writeEfficiencies(const std::string& file, const std::vector<RunConfig>& runs,
                       const std::vector<std::vector<Result>>& results);

/// Geometry check and band decomposition of the cad layer(s) of a configuration (M14, fgsolve
/// --cad-report, Python fg.cad.check) as JSON; stepText (if not empty) replaces the STEP file.
std::string cadReport(const std::string& json, const std::string& baseDir = "", const std::string& stepText = "");
/// SVG of the section at y (FG micrometres) of the first cad layer (Python fg.cad.preview).
std::string cadPreview(const std::string& json, const std::string& baseDir, const std::string& stepText, double y_um);

}  // namespace fg::io
