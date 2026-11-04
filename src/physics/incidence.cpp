#include "fg/physics/incidence.hpp"

#include <cmath>

namespace fg {

IncidentWave makeIncidence(double nu, double phi, double psi, double delta) {
  IncidentWave w;
  w.alpha0 = nu * std::sin(phi) * std::cos(psi);
  w.beta0 = nu * std::cos(phi);
  w.gamma0 = nu * std::sin(phi) * std::sin(psi);
  w.A << std::sin(psi) * std::cos(delta) + std::sin(delta) * std::cos(phi) * std::cos(psi),
      std::sin(delta) * std::sin(phi),
      std::cos(phi) * std::sin(psi) * std::sin(delta) - std::cos(delta) * std::cos(psi);
  return w;
}

}  // namespace fg
