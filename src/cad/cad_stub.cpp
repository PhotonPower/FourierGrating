// FourierGrating - CAD slicer without Nurbforge (FG_ENABLE_NURBFORGE=OFF).
#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"

namespace fg::cad {

namespace {
[[noreturn]] void unavailable() {
  throw NotImplemented("cad layers need a build with FG_ENABLE_NURBFORGE=ON (M14, docs/cad.md)");
}
}  // namespace

bool available() {
  return false;
}
Report check(const CadOptions&, const Stack&) {
  unavailable();
}
std::vector<LayerEntry> slice(const CadOptions&, const Stack&, Report*) {
  unavailable();
}
std::string previewSvg(const CadOptions&, const Stack&, double) {
  unavailable();
}

}  // namespace fg::cad
