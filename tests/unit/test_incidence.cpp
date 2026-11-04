#include <gtest/gtest.h>

#include <cmath>

#include "fg/physics/incidence.hpp"

TEST(Incidence, WaveVectorNormAndTransversality) {
  const double nu = 1.5, phi = 0.61, psi = 0.35;
  for (double delta : {0.0, 0.4, 1.5707963267948966}) {
    const auto w = fg::makeIncidence(nu, phi, psi, delta);
    const fg::Vec3D k(w.alpha0, -w.beta0, w.gamma0);
    EXPECT_NEAR(k.norm(), nu, 1e-14);   // book IX.1 with cos(Psi) would fail this
    EXPECT_NEAR(w.A.norm(), 1.0, 1e-14);
    EXPECT_NEAR(w.A.dot(k), 0.0, 1e-14);
  }
}

TEST(Incidence, NormalIncidenceDeltaIsAngleToZ) {
  const auto w = fg::makeIncidence(1.0, 0.0, fg::PI, 0.0);
  EXPECT_NEAR(w.A.z(), 1.0, 1e-14);
}
