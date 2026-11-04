// FourierGrating - power flux and diffraction efficiencies in homogeneous isotropic media
// (Plan eq:Sy, eq:eff, s/p decomposition).
#pragma once

#include "fg/core/types.hpp"

namespace fg {

struct PlaneWaveOrder {
  double Sy = 0;        ///< time-averaged flux 1/2 Re(Ez Hx* - Ex Hz*) (sign = direction)
  Vec3C E;              ///< complex electric field vector of the order
  cplx amp_s = 0, amp_p = 0;
  bool propagating = false;
};

/// Order with wave vector (alpha, dir*beta, gamma) in a medium eps and tangential E (Ex, Ez).
/// dir = +1 (up-going) or -1 (down-going). Normalized units, mu = 1.
PlaneWaveOrder planeWaveOrder(cplx eps, cplx alpha, cplx gamma, cplx beta, int dir, cplx Ex,
                              cplx Ez);

}  // namespace fg
