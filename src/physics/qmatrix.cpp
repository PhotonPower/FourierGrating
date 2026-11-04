#include "fg/physics/qmatrix.hpp"

#include "fg/core/error.hpp"

namespace fg {

QBlocks qLaurent(const SliceCoefficients&, const IndexMap&) {
  throw NotImplemented("qLaurent - milestone M2");
}
QBlocks qLiLamellar(const SliceCoefficients&, const IndexMap&) {
  throw NotImplemented("qLiLamellar - milestone M2");
}
QBlocks qFFFIsotropic(const SliceCoefficients&, const IndexMap&) {
  throw NotImplemented("qFFFIsotropic - milestone M3, Plan eq:Qiso");
}
QBlocks qFFFAnisotropicA(const SliceCoefficients&, const IndexMap&) {
  throw NotImplemented("qFFFAnisotropicA - milestone M5, Plan eq:QaniA");
}
QBlocks qFFFAnisotropicB(const SliceCoefficients&, const IndexMap&) {
  throw NotImplemented("qFFFAnisotropicB - milestone M5, Plan eq:QaniB");
}

}  // namespace fg
