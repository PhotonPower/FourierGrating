#include "fg/field/nearfield.hpp"

#include "fg/core/error.hpp"

namespace fg {

void NearFieldWriter::write(const FieldCut&, const std::string&, const std::string&) {
  throw NotImplemented("NearFieldWriter::write - milestone M6");
}

}  // namespace fg
