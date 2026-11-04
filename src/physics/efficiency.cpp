#include "fg/physics/efficiency.hpp"

#include <cmath>

namespace fg {

PlaneWaveOrder planeWaveOrder(cplx eps, cplx alpha, cplx gamma, cplx beta, int dir, cplx Ex,
                              cplx Ez) {
  PlaneWaveOrder o;
  const cplx ky = static_cast<double>(dir) * beta;
  // H_t = (1/k_y) K E_t (Plan eq:Kmat)
  const cplx Hx = (alpha * gamma * Ex + (eps - alpha * alpha) * Ez) / ky;
  const cplx Hz = (-(eps - gamma * gamma) * Ex - alpha * gamma * Ez) / ky;
  o.Sy = 0.5 * std::real(Ez * std::conj(Hx) - Ex * std::conj(Hz));
  const cplx Ey = -(alpha * Ex + gamma * Ez) / ky;
  o.E << Ex, Ey, Ez;
  o.propagating = std::abs(beta.imag()) < 1e-12 && beta.real() > 0 && std::abs(eps.imag()) < 1e-14;
  // s/p basis: s = (gamma, 0, -alpha)/|kappa| (z for kappa = 0), p = k_hat x s
  const double kx = alpha.real(), kz = gamma.real();
  const double kap = std::hypot(kx, kz);
  Vec3D s = kap < 1e-12 ? Vec3D(0, 0, 1) : Vec3D(kz / kap, 0, -kx / kap);
  Vec3D k(kx, ky.real(), kz);
  if (k.norm() > 0) k.normalize();
  const Vec3D p = k.cross(s);
  o.amp_s = (o.E.array() * s.cast<cplx>().array()).sum();  // E . s (Eigen dot() conjugates)
  o.amp_p = (o.E.array() * p.cast<cplx>().array()).sum();
  return o;
}

}  // namespace fg
