// FourierGrating - prepared computation for repeated solves (Plan section "M12: API").
#pragma once

#include <memory>
#include <string>

#include "fg/io/config.hpp"

namespace fg {

/// Configuration parsed once; the layer models (wavelength independent: Fourier coefficients,
/// normal fields, Q blocks) are kept between solve() calls. setParameter() changes one number of
/// the configuration (JSON pointer as in sweeps) and rebuilds only the affected layers; results
/// are bit-identical to a Simulation built from the changed configuration. Results equal
/// makeSolver(config().solver)->solve(...). Not thread-safe: one Simulation per thread.
class Simulation {
 public:
  /// Configuration in the JSON format of docs/config_format.md ("sweep" is ignored).
  explicit Simulation(const std::string& json);
  ~Simulation();
  Simulation(Simulation&&) noexcept;
  Simulation& operator=(Simulation&&) noexcept;

  /// Diffraction efficiencies for one incidence (angles in radians, as in Wavelength).
  Result solve(const Wavelength& wl);
  /// First scan point of the configuration.
  Result solve();
  /// Sets the number at a JSON pointer (e.g. "/layers/0/thickness_um") and rebuilds what depends
  /// on it: only the changed layers if nothing outside "/layers" (except the incidence) changed,
  /// otherwise all layer models.
  void setParameter(const std::string& pointer, double value);

  const io::RunConfig& config() const;
  /// Current configuration as JSON text.
  std::string json() const;
  /// Number of layer models built so far (rebuild accounting).
  long modelBuilds() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace fg
