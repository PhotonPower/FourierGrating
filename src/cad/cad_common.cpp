// FourierGrating - CAD slicer: parts independent of Nurbforge (report JSON, LLM tool description).
#include <cmath>
#include <nlohmann/json.hpp>

#include "fg/cad/cad.hpp"

namespace fg::cad {

namespace {
nlohmann::json findings(const std::vector<Finding>& v) {
  nlohmann::json a = nlohmann::json::array();
  for (const auto& f : v) {
    nlohmann::json j{{"code", f.code}, {"message", f.message}};
    if (std::isfinite(f.where_um[0])) j["where_um"] = {f.where_um[0], f.where_um[1], f.where_um[2]};
    a.push_back(j);
  }
  return a;
}
}  // namespace

std::string Report::json() const {
  nlohmann::json j;
  j["ok"] = ok;
  j["errors"] = findings(errors);
  j["warnings"] = findings(warnings);
  j["cell_origin_um"] = {cellOrigin_um[0], cellOrigin_um[1]};
  j["bands"] = nlohmann::json::array();
  for (const auto& b : bands) {
    nlohmann::json bj{{"y0_um", b.y0_um}, {"y1_um", b.y1_um}, {"kind", b.kind}};
    if (!b.material.empty()) bj["material"] = b.material;
    nlohmann::json fr = nlohmann::json::object();
    for (const auto& [m, f] : b.areaFractions) fr[m] = f;
    bj["area_fractions"] = fr;
    j["bands"].push_back(bj);
  }
  return j.dump(2);
}

std::string llmToolSpec() {
  // JSON schema of the check as a tool (function calling); the model passes the STEP text it wrote
  // with Nurbforge (unit micrometre, named solids) and the FourierGrating configuration around the
  // cad block
  const nlohmann::json spec = {
      {"name", "fouriergrating_check_cad"},
      {"description",
       "Checks a STEP model (written with Nurbforge, unit micrometre, one named solid per body) "
       "for "
       "FourierGrating's CAD slicer and returns ok, errors and warnings with fixed codes "
       "(SHELL_OPEN, "
       "SHELL_ORIENTATION, BODY_OVERLAP, PERIODIC_OVERLAP, OUT_OF_RANGE, UNIT_SUSPECT, "
       "FEATURE_SMALL, "
       "BODY_UNASSIGNED, STEP_UNSUPPORTED), positions in micrometres (FG axes: y = stacking "
       "direction = CAD "
       "+Z, x = CAD +X, z = -CAD Y by default) and the band decomposition. Fix every error and "
       "check again."},
      {"parameters",
       {{"type", "object"},
        {"properties",
         {{"step", {{"type", "string"}, {"description", "STEP file contents (ISO 10303-21)"}}},
          {"config",
           {{"type", "object"},
            {"description",
             "FourierGrating configuration (docs/config_format.md): period_um, materials, "
             "superstrate, "
             "substrate, incidence and a layer {\"type\": \"cad\", \"bodies\": [{\"name\": ..., "
             "\"material\": ...}], "
             "\"background\": ...}"}}}}},
        {"required", {"step", "config"}}}}};
  return spec.dump(2);
}

}  // namespace fg::cad
