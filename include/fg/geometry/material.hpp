// FourierGrating - materials (relative permittivity; mu = 1, Plan section "Nicht-Ziele").
#pragma once

#include <array>
#include <memory>
#include <string>

#include "fg/core/types.hpp"

namespace fg {

class Dispersion;

struct Material {
  std::string name;
  bool isotropic = true;
  cplx eps = 1.0;                         ///< isotropic relative permittivity
  Mat3C epsTensor = Mat3C::Identity();    ///< anisotropic relative permittivity (x, y, z)

  /// Dispersive materials (Plan section "Dispersive Materialien (M22)"): the isotropic value (n or eps) or
  /// the principal values as functions of lambda; eps / epsTensor are then NaN until at() resolves them.
  std::shared_ptr<const Dispersion> model;
  std::array<std::shared_ptr<const Dispersion>, 3> principalModels;
  bool principalAxes = false;             ///< given as eps_principal (principal, euler valid)
  std::array<cplx, 3> principal{};        ///< principal values (dispersive entries: after at())
  std::array<double, 3> euler{};          ///< z-x-z Euler angles of the principal axes (radians)

  static Material fromIndex(const std::string& name, cplx n) {
    Material m;
    m.name = name;
    m.eps = n * n;
    m.epsTensor = m.eps * Mat3C::Identity();
    return m;
  }
  static Material fromTensor(const std::string& name, const Mat3C& eps) {
    Material m;
    m.name = name;
    m.isotropic = false;
    m.epsTensor = eps;
    m.eps = eps.trace() / 3.0;
    return m;
  }
  Mat3C tensor() const { return isotropic ? Mat3C(eps * Mat3C::Identity()) : epsTensor; }

  /// True if a value depends on the wavelength.
  bool dispersive() const;
  /// The material with fixed values at lambda (Plan sec:disp_arch); bit-identical to the configuration with
  /// these values written as numbers. Fixed materials are returned unchanged.
  Material at(double lambda_um) const;
};

}  // namespace fg
