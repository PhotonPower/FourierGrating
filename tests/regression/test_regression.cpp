// Regression protection for version 2 (Plan section "Querschnittsthemen der Version 2", AGENTS.md
// section 7): efficiencies of all examples with default options against the golden files of v0.10.
//   default                  1e-11 of the incident power, i.e. absolute on the efficiencies (other
//                            compiler / platform; measured: MSVC vs MinGW GCC <= 1.2e-14, Linux GCC 13
//                            vs MinGW GCC 16 <= 1.3e-12 for adaptive RK and the 2D hole array)
//   FG_REGRESSION_EXACT=1    bit-identical (same machine and compiler as the golden files)
//   FG_REGRESSION_UPDATE=1   rewrite the golden files (only with a reason in CHANGELOG.md)
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "fg/io/config.hpp"

namespace {

const std::string kSrc = FG_SOURCE_DIR;

bool envSet(const char* name) {
  const char* v = std::getenv(name);
  return v && std::string(v) == "1";
}

std::vector<std::string> exampleList() {
  std::ifstream in(kSrc + "/tests/regression/examples.txt");
  std::vector<std::string> names;
  for (std::string s; std::getline(in, s);) {
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ')) s.pop_back();
    if (!s.empty() && s[0] != '#') names.push_back(s);
  }
  return names;
}

std::string num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

/// One line per order (and one "A" line per scan point with the absorption), 17 significant digits.
std::vector<std::string> compute(const std::string& name) {
  const auto runs = fg::io::loadSweep(kSrc + "/examples/" + name + ".json");
  std::vector<std::string> lines;
  for (const auto& run : runs) {
    const auto res = fg::solveScan(run.solver, run.stack, run.scan);
    for (std::size_t i = 0; i < res.size(); ++i) {
      const auto& w = run.scan[i];
      const std::string pre = num(run.sweepValue) + ',' + num(w.lambda_um) + ',' + num(w.phi) + ',' + num(w.psi) +
                              ',' + num(w.delta) + ',';
      for (const char side : {'R', 'T'})
        for (const auto& o : side == 'R' ? res[i].reflected : res[i].transmitted)
          lines.push_back(pre + side + ',' + std::to_string(o.n) + ',' + std::to_string(o.m) + ',' +
                          num(o.efficiency) + ',' + num(o.eff_s) + ',' + num(o.eff_p));
      lines.push_back(pre + "A,0,0," + num(res[i].absorption) + ",0,0");
    }
  }
  return lines;
}

std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> f;
  std::stringstream ss(s);
  for (std::string t; std::getline(ss, t, ',');) f.push_back(t);
  return f;
}

class Regression : public ::testing::TestWithParam<std::string> {};

TEST_P(Regression, EfficienciesEqualGolden) {
  const std::string name = GetParam();
  const std::string file = kSrc + "/tests/regression/golden/" + name + ".csv";
  const auto lines = compute(name);
  if (envSet("FG_REGRESSION_UPDATE")) {
    std::ofstream out(file, std::ios::binary);
    out << "sweep_value,lambda_um,phi_rad,psi_rad,delta_rad,side,n,m,efficiency,eff_s,eff_p\n";
    for (const auto& l : lines) out << l << '\n';
    GTEST_SKIP() << "golden file written: " << file;
  }
  std::ifstream in(file);
  ASSERT_TRUE(in.good()) << "missing golden file " << file << " (FG_REGRESSION_UPDATE=1)";
  std::vector<std::string> golden;
  std::string header;
  std::getline(in, header);
  for (std::string s; std::getline(in, s);) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
    if (!s.empty()) golden.push_back(s);
  }
  ASSERT_EQ(lines.size(), golden.size()) << name;
  const bool exact = envSet("FG_REGRESSION_EXACT");
  for (std::size_t k = 0; k < lines.size(); ++k) {
    if (exact) {
      EXPECT_EQ(lines[k], golden[k]) << name << " line " << k + 2;
      continue;
    }
    const auto a = split(lines[k]), b = split(golden[k]);
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t c = 0; c < a.size(); ++c) {
      if (c >= 5 && c <= 7) {  // side, n, m
        EXPECT_EQ(a[c], b[c]) << name << " line " << k + 2;
        continue;
      }
      const double x = std::stod(a[c]), y = std::stod(b[c]);
      // angles and lambda are inputs (exact); efficiencies are fractions of the incident power
      EXPECT_LE(std::abs(x - y), c <= 4 ? 0.0 : 1e-11) << name << " line " << k + 2 << " column " << c;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Examples, Regression, ::testing::ValuesIn(exampleList()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

}  // namespace
