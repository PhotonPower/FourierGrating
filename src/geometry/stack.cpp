#include "fg/geometry/stack.hpp"

#include <limits>

#include "fg/core/error.hpp"
#include "fg/geometry/dispersion.hpp"
#include "fg/physics/anisotropic.hpp"

namespace fg {

const Material& Stack::material(const std::string& name) const {
  const auto it = materials.find(name);
  FG_REQUIRE(it != materials.end(), "unknown material '" + name + "'");
  return it->second;
}

bool Material::dispersive() const {
  if (model) return true;
  for (const auto& p : principalModels)
    if (p) return true;
  return false;
}

Material Material::at(double lambda_um) const {
  if (!dispersive()) return *this;
  // the same constructions as the configuration parser for numbers (bit-identical, Plan V23a)
  if (model) {
    Material m;
    const cplx v = model->value(lambda_um);
    if (model->quantity() == Dispersion::Quantity::N) {
      m = fromIndex(name, v);
    } else {
      m = fromIndex(name, 1.0);
      m.eps = v;
      m.epsTensor = m.eps * Mat3C::Identity();
    }
    return m;
  }
  std::array<cplx, 3> p = principal;
  for (std::size_t i = 0; i < 3; ++i)
    if (principalModels[i]) p[i] = principalModels[i]->eps(lambda_um);
  Material m = fromTensor(name, rotatedTensor(p[0], p[1], p[2], euler[0], euler[1], euler[2]));
  m.principalAxes = true;
  m.principal = p;
  m.euler = euler;
  return m;
}

bool Stack::dispersive() const {
  for (const auto& [name, m] : materials)
    if (m.dispersive()) return true;
  return false;
}

Stack Stack::at(double lambda_um) const {
  Stack s = *this;
  for (auto& [name, m] : s.materials)
    if (m.dispersive()) m = m.at(lambda_um);
  return s;
}

}  // namespace fg
