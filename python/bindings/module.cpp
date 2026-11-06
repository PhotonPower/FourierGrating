// fouriergrating._core - Python bindings (Plan section "M12: API"). Conversion and calls only
// (AGENTS.md section 7): JSON text in, Python objects / numpy arrays out; degrees <-> radians at
// this boundary; the GIL is released during every computation.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <chrono>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <vector>

#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/field/nearfield.hpp"
#include "fg/cad/cad.hpp"
#include "fg/io/config.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/solvers/simulation.hpp"
#include "fg/version.hpp"

namespace py = pybind11;

namespace {

constexpr double kDeg = fg::PI / 180.0;

py::list orders(const std::vector<fg::OrderEfficiency>& v) {
  py::list l;
  for (const auto& o : v)
    l.append(py::make_tuple(o.n, o.m, o.efficiency, o.eff_s, o.eff_p, o.amp_s, o.amp_p));
  return l;
}

/// Scope times of fg::prof (M15) as {name: {"total_s", "count", "max_s"}}.
py::dict scopes(const std::vector<fg::prof::Entry>& entries) {
  py::dict d;
  for (const auto& e : entries) {
    py::dict x;
    x["total_s"] = e.total_s;
    x["count"] = e.count;
    x["max_s"] = e.max_s;
    d[py::str(e.name)] = x;
  }
  return d;
}

/// Runs f with fg::prof enabled (if profile) and returns the collected entries.
template <typename F>
std::vector<fg::prof::Entry> profiled(bool profile, F&& f) {
  if (!profile) {
    f();
    return {};
  }
  fg::prof::reset();
  fg::prof::enable(true);
  try {
    f();
  } catch (...) {
    fg::prof::enable(false);
    throw;
  }
  fg::prof::enable(false);
  return fg::prof::snapshot();
}

py::dict toDict(const fg::Result& r, double seconds, const std::vector<fg::prof::Entry>& prof = {}) {
  py::dict d;
  d["reflected"] = orders(r.reflected);
  d["transmitted"] = orders(r.transmitted);
  d["absorption"] = r.absorption;
  py::dict diag;
  diag["energy_balance"] = r.diag.energyBalance;
  diag["max_condition_t11"] = r.diag.maxConditionT11;
  diag["rk_steps"] = r.diag.rkSteps;
  diag["warnings"] = r.diag.warnings;
  diag["precision"] = r.diag.precision;  // M16: working precision actually used
  diag["precision_fallback"] = r.diag.precisionFallback;
  diag["refinement_steps"] = r.diag.refinementSteps;
  diag["refinement_fallbacks"] = r.diag.refinementFallbacks;
  py::dict timings;
  timings["solve_s"] = seconds;
  if (!prof.empty()) timings["scopes"] = scopes(prof);
  diag["timings"] = timings;
  d["diagnostics"] = diag;
  return d;
}

py::dict wavelength(const fg::Wavelength& w) {
  py::dict d;
  d["lambda_um"] = w.lambda_um;
  d["phi_deg"] = w.phi / kDeg;
  d["psi_deg"] = w.psi / kDeg;
  d["delta_deg"] = w.delta / kDeg;
  return d;
}

template <typename T>
py::array_t<T> array(const std::vector<T>& v, const fg::FieldCut& c) {
  const std::size_t n1 = c.axis1_um.size(), n2 = c.axis2_um.size();
  std::vector<py::ssize_t> shape =
      n2 ? std::vector<py::ssize_t>{static_cast<py::ssize_t>(n2), static_cast<py::ssize_t>(n1)}
         : std::vector<py::ssize_t>{static_cast<py::ssize_t>(n1)};
  py::array_t<T> a(shape);
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

py::dict cutDict(const fg::FieldCut& c) {
  py::dict d;
  d["cut"] = c.cut;
  d["axis1_name"] = c.axis1Name;
  d["axis2_name"] = c.axis2Name;
  d["axis1_um"] =
      py::array_t<double>(static_cast<py::ssize_t>(c.axis1_um.size()), c.axis1_um.data());
  d["axis2_um"] =
      py::array_t<double>(static_cast<py::ssize_t>(c.axis2_um.size()), c.axis2_um.data());
  d["lambda_um"] = c.lambda_um;
  const char* names[3] = {"x", "y", "z"};
  for (int k = 0; k < 3; ++k) {
    d[py::str(std::string("E") + names[k])] = array(c.E[k], c);
    d[py::str(std::string("H") + names[k])] = array(c.H[k], c);
  }
  const fg::DerivedFields f = fg::derivedFields(c);
  d["absE2"] = array(f.absE2, c);
  d["absH2"] = array(f.absH2, c);
  d["Sy"] = array(f.Sy, c);
  if (!c.material.empty()) d["material"] = array(c.material, c);
  d["meta"] = fg::nearFieldMeta(c);
  return d;
}

template <typename F>
auto timed(F&& f, double& seconds) {
  const auto t0 = std::chrono::steady_clock::now();
  auto r = f();
  seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

/// Incidence of a Simulation call: given values in degrees / micrometres, the rest from the
/// first scan point of the configuration.
fg::Wavelength incidence(const fg::Simulation& s, std::optional<double> lambda,
                         std::optional<double> phi, std::optional<double> psi,
                         std::optional<double> delta) {
  const auto& scan = s.config().scan;
  FG_REQUIRE(!scan.empty(), "configuration has no incidence");
  fg::Wavelength w = scan.front();
  if (lambda) w.lambda_um = *lambda;
  if (phi) w.phi = *phi * kDeg;
  if (psi) w.psi = *psi * kDeg;
  if (delta) w.delta = *delta * kDeg;
  return w;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() =
      "FourierGrating core (C++); use the package " FG_PYTHON_PACKAGE ", not this module directly";

  // fg::Error -> FourierGratingError (RuntimeError); translators registered later are tried first
  py::register_exception<fg::Error>(m, "FourierGratingError", PyExc_RuntimeError);
  static PyObject* fgError = m.attr("FourierGratingError").ptr();  // owned by the module
  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p) std::rethrow_exception(p);
    } catch (const fg::NotImplemented& e) {
      PyErr_SetString(PyExc_NotImplementedError, e.what());
    } catch (const nlohmann::json::exception& e) {  // malformed or incomplete configuration
      PyErr_SetString(fgError, e.what());
    }
  });

  m.def("version", [] { return std::string(FG_VERSION); });
#ifdef FG_HAVE_CUDA
  m.def("variant", [] { return std::string("cuda12"); });
#else
  m.def("variant", [] { return std::string("cpu"); });
#endif
  m.def("cuda_available", [] { return fg::linalg::cudaAvailable(); });

  // CAD slicer (M14): check, preview, LLM tool description
  m.def("cad_available", [] { return fg::cad::available(); });
  m.def(
      "cad_report",
      [](const std::string& json, const std::string& base_dir, const std::string& step) {
        py::gil_scoped_release nogil;
        return fg::io::cadReport(json, base_dir, step);
      },
      py::arg("config_json"), py::arg("base_dir") = "", py::arg("step_text") = "");
  m.def(
      "cad_preview",
      [](const std::string& json, const std::string& base_dir, const std::string& step, double y_um) {
        py::gil_scoped_release nogil;
        return fg::io::cadPreview(json, base_dir, step, y_um);
      },
      py::arg("config_json"), py::arg("base_dir"), py::arg("step_text"), py::arg("y_um"));
  m.def("cad_llm_tool_spec", [] { return fg::cad::llmToolSpec(); });

  m.def(
      "solve",
      [](const std::string& json, bool profile) {
        double t = 0;
        fg::Result r;
        std::vector<fg::prof::Entry> prof;
        {
          py::gil_scoped_release nogil;
          prof = profiled(profile, [&] { r = timed([&] { return fg::Simulation(json).solve(); }, t); });
        }
        return toDict(r, t, prof);
      },
      py::arg("config_json"), py::arg("profile") = false, "First scan point of the configuration.");

  m.def(
      "solve_scan",
      [](const std::string& json, bool profile) {
        struct Run {
          fg::io::RunConfig cfg;
          std::vector<fg::Result> res;
          double seconds = 0;
          std::vector<fg::prof::Entry> prof;
        };
        std::vector<Run> runs;
        {
          py::gil_scoped_release nogil;
          for (auto& c : fg::io::parseSweep(json)) {
            Run r{std::move(c), {}, 0, {}};
            r.prof = profiled(profile, [&] {
              r.res = timed([&] { return fg::solveScan(r.cfg.solver, r.cfg.stack, r.cfg.scan); }, r.seconds);
            });
            runs.push_back(std::move(r));
          }
        }
        py::list out;
        for (const auto& r : runs) {
          py::dict d;
          d["sweep_parameter"] = r.cfg.sweepParameter;
          d["sweep_value"] = r.cfg.sweepValue;
          py::list pts;
          for (std::size_t i = 0; i < r.res.size(); ++i) {
            py::dict p = toDict(r.res[i], r.seconds / static_cast<double>(r.res.size()), r.prof);
            p["incidence"] = wavelength(r.cfg.scan[i]);
            pts.append(p);
          }
          d["points"] = pts;
          out.append(d);
        }
        return out;
      },
      py::arg("config_json"), py::arg("profile") = false,
      "All sweep values and scan points (as fgsolve); profile: scope times of the whole run in every point.");

  m.def(
      "near_field",
      [](const std::string& json) {
        fg::NearFieldResult nf;
        fg::Wavelength w;
        {
          py::gil_scoped_release nogil;
          const auto c = fg::io::parseConfig(json);
          FG_REQUIRE(!c.scan.empty(), "configuration has no incidence");
          w = c.scan.front();
          nf = fg::makeSolver(c.solver)->nearField(c.stack, w, c.field);
        }
        py::dict d;
        d["horizontal"] = cutDict(nf.horizontal);
        d["vertical"] = cutDict(nf.vertical);
        d["incidence"] = wavelength(w);
        return d;
      },
      py::arg("config_json"), "Near-field cuts for the first scan point (block \"nearfield\").");

  py::class_<fg::Simulation>(m, "Simulation")
      .def(py::init([](const std::string& json) {
             py::gil_scoped_release nogil;
             return std::make_unique<fg::Simulation>(json);
           }),
           py::arg("config_json"))
      .def(
          "solve",
          [](fg::Simulation& s, std::optional<double> lambda_um, std::optional<double> phi_deg,
             std::optional<double> psi_deg, std::optional<double> delta_deg) {
            const fg::Wavelength w = incidence(s, lambda_um, phi_deg, psi_deg, delta_deg);
            double t = 0;
            fg::Result r;
            {
              py::gil_scoped_release nogil;
              r = timed([&] { return s.solve(w); }, t);
            }
            py::dict d = toDict(r, t);
            d["incidence"] = wavelength(w);
            return d;
          },
          py::arg("lambda_um") = py::none(), py::arg("phi_deg") = py::none(),
          py::arg("psi_deg") = py::none(), py::arg("delta_deg") = py::none())
      .def(
          "set_parameter",
          [](fg::Simulation& s, const std::string& pointer, double value) {
            py::gil_scoped_release nogil;
            s.setParameter(pointer, value);
          },
          py::arg("pointer"), py::arg("value"))
      .def("config_json", &fg::Simulation::json)
      .def_property_readonly("model_builds", &fg::Simulation::modelBuilds);
}
