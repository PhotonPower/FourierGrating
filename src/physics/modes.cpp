#include "fg/physics/modes.hpp"

#include <cmath>

#include "fg/core/error.hpp"

namespace fg {

cplx betaBranch(cplx eps, cplx alpha, cplx gamma, double eta) {
  cplx b = std::sqrt(eps - alpha * alpha - gamma * gamma);
  if (b.imag() < 0.0 || (b.imag() == 0.0 && b.real() < 0.0)) b = -b;
  if (std::abs(b) < eta) b += cplx(0.0, eta);
  return b;
}

ModeBasis isotropicModes(cplx eps, const VecC& alpha, const VecC& gamma, double eta) {
  const Eigen::Index n = alpha.size();
  FG_REQUIRE(gamma.size() == n, "alpha/gamma size mismatch");
  ModeBasis mb;
  mb.beta.resize(n);
  // V maps (Ex;Ez) -> (Hx;Hz) for up-going waves: V = K / beta (Plan eq:Kmat)
  MatC V = MatC::Zero(2 * n, 2 * n), Vi = MatC::Zero(2 * n, 2 * n);
  for (Eigen::Index r = 0; r < n; ++r) {
    const cplx a = alpha(r), g = gamma(r);
    const cplx b = betaBranch(eps, a, g, eta);
    mb.beta(r) = b;
    V(r, r) = a * g / b;
    V(r, n + r) = (eps - a * a) / b;
    V(n + r, r) = -(eps - g * g) / b;
    V(n + r, n + r) = -a * g / b;
    const cplx f = 1.0 / (eps * b);
    Vi(r, r) = -a * g * f;
    Vi(r, n + r) = -(eps - a * a) * f;
    Vi(n + r, r) = (eps - g * g) * f;
    Vi(n + r, n + r) = a * g * f;
  }
  const MatC Id = MatC::Identity(2 * n, 2 * n);
  mb.Psi.resize(4 * n, 4 * n);
  mb.Psi << Id, Id, -V, V;
  mb.PsiInv.resize(4 * n, 4 * n);
  mb.PsiInv << 0.5 * Id, -0.5 * Vi, 0.5 * Id, 0.5 * Vi;
  mb.lambda.resize(4 * n);
  for (Eigen::Index r = 0; r < n; ++r) {
    mb.lambda(r) = mb.lambda(n + r) = -I * mb.beta(r);              // down (a_x, a_z)
    mb.lambda(2 * n + r) = mb.lambda(3 * n + r) = I * mb.beta(r);   // up   (b_x, b_z)
  }
  return mb;
}

ModeBasis anisotropicModes(const Mat3C&, const VecC&, const VecC&) {
  throw NotImplemented("anisotropicModes (milestone M5, Plan eq:classify)");
}

double poyntingY(const VecC& F) {
  const Eigen::Index n = F.size() / 4;
  double s = 0.0;
  for (Eigen::Index r = 0; r < n; ++r) {
    const cplx Ex = F(r), Ez = F(n + r), Hx = F(2 * n + r), Hz = F(3 * n + r);
    s += 0.5 * std::real(Ez * std::conj(Hx) - Ex * std::conj(Hz));
  }
  return s;
}

}  // namespace fg
