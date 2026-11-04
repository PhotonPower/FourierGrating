#include "fg/physics/anisotropic.hpp"

#include <cmath>

#include "fg/core/error.hpp"

namespace fg {

namespace {
// Plain (bilinear) cross product. Eigen's cross() conjugates its result for complex vectors!
Vec3C crossBilinear(const Vec3C& a, const Vec3C& b) {
  return Vec3C(a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0));
}
}  // namespace

Mat3C cMatrix(const Mat3C& eps, const Vec3D& N, const Vec3D& T1, const Vec3D& T2) {
  const Vec3C Nc = N.cast<cplx>();
  const Vec3C v = (Nc.transpose() * eps).transpose();  // v_j = sum_i N_i eps_ij
  const cplx xi = (v.array() * Nc.array()).sum();      // xi = sum_i v_i N_i (no conjugation!)
  FG_REQUIRE(std::abs(xi) > 1e-14 * eps.cwiseAbs().maxCoeff(), "xi = N eps N vanishes (Plan eq:Ceps)");
  const Vec3C c1 = crossBilinear(v, T2.cast<cplx>()) / xi;
  const Vec3C c3 = -crossBilinear(v, T1.cast<cplx>()) / xi;
  Mat3C C;
  C.col(0) = c1;
  C.col(1) = Nc / xi;
  C.col(2) = c3;
  return C;
}

ProjectorSet projectorMatrices(const Mat3C& eps, const Vec3D& N) {
  const Mat3C PN = (N * N.transpose()).cast<cplx>();
  const Vec3C Nc = N.cast<cplx>();
  const cplx xi = (Nc.transpose() * eps * Nc)(0, 0);
  FG_REQUIRE(std::abs(xi) > 1e-14 * eps.cwiseAbs().maxCoeff(), "xi = N eps N vanishes (Plan eq:projector)");
  ProjectorSet p;
  p.invXi = 1.0 / xi;
  p.H1 = Mat3C::Identity() - p.invXi * PN * eps;
  p.G1 = eps * p.H1;
  p.G2 = p.invXi * eps;
  p.PT = Mat3C::Identity() - PN;
  return p;
}

void profileTangents(const Vec3D& N, Vec3D& T1, Vec3D& T2) {
  const double r = std::hypot(N.x(), N.y());
  FG_REQUIRE(r > 1e-12, "profile normal must have a component in the x-y plane");
  T1 = Vec3D(N.y() / r, -N.x() / r, 0.0);
  T2 = T1.cross(N);  // = z for N_z = 0;  T2 x T1 = N
}

Mat3C rotatedTensor(cplx e1, cplx e2, cplx e3, double phi, double theta, double psi) {
  using std::cos;
  using std::sin;
  Eigen::Matrix3d Rz1, Rx, Rz2;
  Rz1 << cos(phi), -sin(phi), 0, sin(phi), cos(phi), 0, 0, 0, 1;
  Rx << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
  Rz2 << cos(psi), -sin(psi), 0, sin(psi), cos(psi), 0, 0, 0, 1;
  const Eigen::Matrix3d R = Rz1 * Rx * Rz2;
  Mat3C D = Mat3C::Zero();
  D(0, 0) = e1;
  D(1, 1) = e2;
  D(2, 2) = e3;
  return R.cast<cplx>() * D * R.transpose().cast<cplx>();
}

}  // namespace fg
