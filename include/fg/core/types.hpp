// FourierGrating - basic numeric types.
// All physics code works in normalized units (Plan, section "Konventionen"):
// lengths in units of 1/k0, H scaled by Z0, omega = mu0 = eps0 = 1.
#pragma once

#include <Eigen/Dense>
#include <complex>
#include <vector>

namespace fg {

using cplx = std::complex<double>;
using MatC = Eigen::MatrixXcd;
using VecC = Eigen::VectorXcd;
using MatD = Eigen::MatrixXd;
using VecD = Eigen::VectorXd;
using Mat3C = Eigen::Matrix3cd;
using Vec3C = Eigen::Vector3cd;
using Vec3D = Eigen::Vector3d;

inline constexpr cplx I{0.0, 1.0};
inline constexpr double PI = 3.14159265358979323846;

}  // namespace fg
