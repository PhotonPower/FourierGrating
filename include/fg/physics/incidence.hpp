// FourierGrating - incident plane wave (Plan eq:incidence, eq:polarization; corrects book IX.1).
#pragma once

#include "fg/core/types.hpp"

namespace fg {

struct IncidentWave {
  double alpha0 = 0, beta0 = 0, gamma0 = 0;  ///< normalized: k = (alpha0, -beta0, gamma0)
  Vec3D A;                                   ///< unit polarization vector, A . k = 0
};

/// nu_sup: real refractive index of the superstrate; angles in radians.
IncidentWave makeIncidence(double nu_sup, double phi, double psi, double delta);

}  // namespace fg
