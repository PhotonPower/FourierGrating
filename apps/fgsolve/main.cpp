// fgsolve - command line front end of FourierGrating.
//   fgsolve --version
//   fgsolve --info                 (backends, build options)
//   fgsolve --cad-report <config>  (geometry check and band decomposition of the cad layers, JSON)
//   fgsolve <config.json>          (diffraction efficiencies, scans and sweeps, optional near fields;
//                                   near fields for the first scan point of the first sweep value)
//   fgsolve --profile <config.json> (as above, plus <output.dir>/profile.json: time per scope, M15)
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "fg/cad/cad.hpp"
#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/io/config.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/version.hpp"

namespace {

void printResult(const fg::Wavelength& w, const fg::Result& r) {
  const double deg = 180.0 / fg::PI;
  std::printf("lambda = %.6g um, phi = %.4g deg, psi = %.4g deg, delta = %.4g deg\n", w.lambda_um,
              w.phi * deg, w.psi * deg, w.delta * deg);
  std::printf("  side   n    m   efficiency\n");
  for (const auto& o : r.reflected) std::printf("  R   %4d %4d   %.10f\n", o.n, o.m, o.efficiency);
  for (const auto& o : r.transmitted) std::printf("  T   %4d %4d   %.10f\n", o.n, o.m, o.efficiency);
  std::printf("  sum R+T = %.12f   absorption = %.3e\n", r.diag.energyBalance, r.absorption);
  if (r.diag.precision != "fp64" || r.diag.precisionFallback)  // M16
    std::printf("  precision %s%s\n", r.diag.precision.c_str(), r.diag.precisionFallback ? " (repeated in fp64)" : "");
  for (const auto& msg : r.diag.warnings) std::printf("  warning: %s\n", msg.c_str());
  std::printf("\n");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || std::strcmp(argv[1], "--help") == 0) {
    std::cout << "usage: fgsolve [--version | --info | --cad-report <config.json> | [--profile] <config.json>]\n";
    return argc < 2 ? 1 : 0;
  }
  try {
    if (std::strcmp(argv[1], "--version") == 0) {
      std::cout << "FourierGrating " << FG_VERSION << "\n";
      return 0;
    }
    if (std::strcmp(argv[1], "--cad-report") == 0) {
      FG_REQUIRE(argc >= 3, "--cad-report needs a configuration file");
      std::ifstream in(argv[2]);
      FG_REQUIRE(in.good(), std::string("cannot open configuration file ") + argv[2]);
      std::stringstream ss;
      ss << in.rdbuf();
      const std::string report = fg::io::cadReport(ss.str(), std::filesystem::path(argv[2]).parent_path().string());
      std::cout << report << "\n";
      return report.find("\"ok\": true") != std::string::npos ? 0 : 4;
    }
    if (std::strcmp(argv[1], "--info") == 0) {
      std::cout << "FourierGrating " << FG_VERSION << "\n"
                << "cpu backend : " << fg::linalg::makeCpuBackend()->name() << "\n";
#ifdef FG_HAVE_CUDA
      if (fg::linalg::cudaAvailable())
        std::cout << "cuda backend: " << fg::linalg::makeCudaBackend(0)->name() << "\n";
      else
        std::cout << "cuda backend: built, but no CUDA device available (CPU is used)\n";
#else
      std::cout << "cuda backend: not built\n";
#endif
      std::cout << "cad slicer  : " << (fg::cad::available() ? "Nurbforge" : "not built (FG_ENABLE_NURBFORGE=OFF)")
                << "\n";
      return 0;
    }
    const bool profile = std::strcmp(argv[1], "--profile") == 0;
    FG_REQUIRE(!profile || argc >= 3, "--profile needs a configuration file");
    fg::prof::enable(profile);
    const auto t0 = std::chrono::steady_clock::now();
    // one configuration, or one per value of a parameter sweep (M11); scan points run in parallel
    const auto runs = fg::io::loadSweep(argv[profile ? 2 : 1]);
    const auto& cfg = runs.front();
    auto solver = fg::makeSolver(cfg.solver);
    std::vector<std::vector<fg::Result>> all;
    for (const auto& run : runs) {
      all.push_back(fg::solveScan(run.solver, run.stack, run.scan));
      if (!run.sweepParameter.empty()) std::printf("== sweep %s = %.10g\n", run.sweepParameter.c_str(), run.sweepValue);
      for (std::size_t i = 0; i < all.back().size(); ++i) printResult(run.scan[i], all.back()[i]);
    }
    const auto file = (std::filesystem::path(cfg.outputDir) / cfg.efficiencyFile).string();
    fg::io::writeEfficiencies(file, runs, all);
    std::cout << "written: " << file << "\n";
    if (cfg.nearField && !cfg.scan.empty()) {
      const auto nf = solver->nearField(cfg.stack, cfg.scan.front(), cfg.field);
      const auto base = std::filesystem::path(cfg.outputDir);
      fg::NearFieldWriter::write(nf.horizontal, (base / "nearfield_horizontal").string());
      fg::NearFieldWriter::write(nf.vertical, (base / "nearfield_vertical").string());
      std::cout << "written: " << (base / "nearfield_horizontal").string() << ", "
                << (base / "nearfield_vertical").string() << "\n";
    }
    if (profile) {
      // busy time per scope summed over all threads; "wall" is the elapsed time of the whole run
      fg::prof::record("wall", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
      const auto pf = std::filesystem::path(cfg.outputDir) / "profile.json";
      std::ofstream(pf) << fg::prof::json() << "\n";
      std::printf("profile (busy time summed over threads):\n");
      for (const auto& e : fg::prof::snapshot())
        std::printf("  %-26s %10.4f s  %8ld x  max %.4f s\n", e.name.c_str(), e.total_s, e.count, e.max_s);
      std::cout << "written: " << pf.string() << "\n";
    }
    return 0;
  } catch (const fg::NotImplemented& e) {
    std::cerr << "fgsolve: " << e.what() << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "fgsolve: error: " << e.what() << "\n";
    return 3;
  }
}
