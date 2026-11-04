#include "fg/geometry/stack.hpp"

#include "fg/core/error.hpp"

namespace fg {

const Material& Stack::material(const std::string& name) const {
  const auto it = materials.find(name);
  FG_REQUIRE(it != materials.end(), "unknown material '" + name + "'");
  return it->second;
}

}  // namespace fg
