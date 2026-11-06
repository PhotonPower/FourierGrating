// Precision study (Plan section "M16: Tests und Abnahme"): every example with a golden file plus
// the profile cases V9 and V13 (bench/) in fp64, mixed and fp32. Prints a Markdown table (maximal
// deviation of the efficiencies from fp64 in units of the incident power, energy balance, fallback)
// and checks the targets fixed after the study (Plan v2.8): mixed and fp32 <= 1e-4 for V1-V9.
// Backend as configured or from FG_BACKEND / FG_BACKEND_THRESHOLD (GPU: FG_BACKEND=cuda
// FG_BACKEND_THRESHOLD=0).
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "fg/cad/cad.hpp"
#include "fg/io/config.hpp"

namespace {

const std::string kSrc = FG_SOURCE_DIR;

std::vector<std::string> cases() {
  std::vector<std::string> v;
  std::ifstream in(kSrc + "/tests/regression/examples.txt");
  for (std::string s; std::getline(in, s);) {
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ')) s.pop_back();
    if (!s.empty() && s[0] != '#') v.push_back("examples/" + s);
  }
  v.push_back("bench/v9_sinus_x_2d_rk");
  v.push_back("bench/v13_gold_sinus_rk");
  return v;
}

struct Outcome {
  std::vector<double> eff;  ///< all orders and absorption of all scan points
  double energyDev = 0;  ///< max |sum R + T - 1| over lossless points (absorption < 1e-12 in fp64)
  bool fallback = false;
  std::string precision;
};

Outcome run(const std::string& name, const std::string& precision,
            const std::vector<bool>* lossless) {
  auto runs = fg::io::loadSweep(kSrc + "/" + name + ".json");
  Outcome r;
  std::size_t point = 0;
  for (auto& c : runs) {
    c.solver.precision = precision;
    c.solver.precisionFallback = precision == "fp64" ? 0 : -1;  // defaults (on for fp32)
    for (const auto& res : fg::solveScan(c.solver, c.stack, c.scan)) {
      for (const auto* side : {&res.reflected, &res.transmitted})
        for (const auto& o : *side) r.eff.push_back(o.efficiency);
      r.eff.push_back(res.absorption);
      r.fallback = r.fallback || res.diag.precisionFallback;
      r.precision = res.diag.precision;
      if (lossless && (*lossless)[point])
        r.energyDev = std::max(r.energyDev, std::abs(res.diag.energyBalance - 1.0));
      ++point;
    }
  }
  return r;
}

std::vector<bool> losslessPoints(const std::string& name) {
  std::vector<bool> v;
  for (auto& c : fg::io::loadSweep(kSrc + "/" + name + ".json")) {
    c.solver.precision = "fp64";
    for (const auto& res : fg::solveScan(c.solver, c.stack, c.scan))
      v.push_back(std::abs(res.absorption) < 1e-12);
  }
  return v;
}

}  // namespace

TEST(PrecisionStudy, AllCasesInAllModes) {
  std::printf(
      "\n| Fall | mixed: max. Abw. | fp32: max. Abw. | Energie fp32 (verlustfrei) | Rückfall fp32 "
      "|\n");
  std::printf("|---|---|---|---|---|\n");
  for (const std::string& name : cases()) {
    {
      std::ifstream in(kSrc + "/" + name + ".json");
      std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      if (s.find("\"cad\"") != std::string::npos && !fg::cad::available()) continue;
    }
    const auto lossless = losslessPoints(name);
    const Outcome ref = run(name, "fp64", nullptr);
    double dev[2] = {0, 0};
    Outcome r32;
    int k = 0;
    for (const char* p : {"mixed", "fp32"}) {
      const Outcome r = run(name, p, &lossless);
      ASSERT_EQ(r.eff.size(), ref.eff.size()) << name << " " << p;
      for (std::size_t i = 0; i < r.eff.size(); ++i)
        dev[k] = std::max(dev[k], std::abs(r.eff[i] - ref.eff[i]));
      if (k == 1) r32 = r;
      ++k;
    }
    bool anyLossless = false;
    for (bool b : lossless) anyLossless = anyLossless || b;
    std::printf("| %s | %.1e | %.1e | %s | %s |\n", name.c_str(), dev[0], dev[1],
                anyLossless ? (std::to_string(r32.energyDev).c_str()) : "–",
                r32.fallback ? "ja" : "nein");
    // targets fixed after the study (Plan section "M16: Tests und Abnahme", v2.8): V1-V9
    const bool v1to9 =
        name.find("/v1_") != std::string::npos || name.find("/v3_") != std::string::npos ||
        name.find("/v4_") != std::string::npos || name.find("/v5_") != std::string::npos ||
        name.find("/v7_") != std::string::npos || name.find("/v9_") != std::string::npos;
    if (v1to9) {
      EXPECT_LE(dev[0], 1e-4) << name << " mixed";
      EXPECT_LE(dev[1], 1e-4) << name << " fp32";
    }
  }
}
