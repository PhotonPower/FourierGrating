#include "fg/geometry/normal_field.hpp"

#include "fg/core/error.hpp"

namespace fg {

NormalFieldCoefficients normalField(const Profile1DLayer&, double, const IndexMap&) {
  throw NotImplemented("normalField(Profile1DLayer) - milestone M3, Plan eq:N1d");
}

NormalFieldCoefficients normalField(const Polygon2DLayer&, double, double, const IndexMap&) {
  throw NotImplemented("normalField(Polygon2DLayer) - milestone M7, Plan algorithm NF");
}

}  // namespace fg
