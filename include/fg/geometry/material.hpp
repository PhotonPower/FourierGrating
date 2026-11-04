// FourierGrating - materials (relative permittivity; mu = 1, Plan section "Nicht-Ziele").
#pragma once

#include <string>

#include "fg/core/types.hpp"

namespace fg {

struct Material {
  std::string name;
  bool isotropic = true;
  cplx eps = 1.0;                         ///< isotropic relative permittivity
  Mat3C epsTensor = Mat3C::Identity();    ///< anisotropic relative permittivity (x, y, z)

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
};

}  // namespace fg
