#include "fg/io/config.hpp"

#include <fstream>
#include <nlohmann/json.hpp>

#include "fg/core/error.hpp"

namespace fg::io {

RunConfig loadConfig(const std::string& path) {
  std::ifstream in(path);
  FG_REQUIRE(in.good(), "cannot open configuration file " + path);
  const nlohmann::json j = nlohmann::json::parse(in);
  FG_REQUIRE(j.contains("dimension"), "configuration needs 'dimension'");
  throw NotImplemented("loadConfig: schema parsing (milestone M2, docs/config_format.md)");
}

}  // namespace fg::io
