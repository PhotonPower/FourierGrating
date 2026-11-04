// fgbench - benchmarks of the linear algebra backends and of full 2D solves (M10, docs/benchmarks.md).
//   fgbench [--sizes 256,512,1024] [--orders 4,6,8] [--repeat 1]
// Prints Markdown tables: ZGEMM / LU solve / general eigenproblem per backend, and the time of one
// RCWA solve of the Ag hole array (V10) with backend cpu and cuda.
#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "fg/linalg/backend.hpp"
#include "fg/solvers/solver.hpp"
#include "fg/version.hpp"

namespace {

std::vector<int> parseList(const char* s) {
  std::vector<int> v;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ',')) v.push_back(std::stoi(t));
  return v;
}

double seconds(const std::function<void()>& f, int repeat) {
  double best = 1e300;
  for (int r = 0; r < repeat; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  return best;
}

fg::Stack holeArray() {
  fg::Stack st;
  st.dx_um = st.dz_um = 0.9;
  st.materials["air"] = fg::Material::fromIndex("air", 1.0);
  st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["ag"] = fg::Material::fromIndex("ag", fg::cplx(0.1, 8.94));
  st.superstrate = "air";
  st.substrate = "glass";
  fg::Polygon2DLayer g;
  g.thickness_um = 0.2;
  g.background = "ag";
  g.inclusions.push_back({{{0.325, 0.325}, {0.575, 0.325}, {0.575, 0.575}, {0.325, 0.575}}, "air"});
  fg::LayerEntry le;
  le.spec = g;
  st.layers = {le};
  return st;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<int> sizes{256, 512, 1024}, orders{4, 6, 8};
  int repeat = 1;
  for (int i = 1; i + 1 < argc; i += 2) {
    if (!std::strcmp(argv[i], "--sizes")) sizes = parseList(argv[i + 1]);
    else if (!std::strcmp(argv[i], "--orders")) orders = parseList(argv[i + 1]);
    else if (!std::strcmp(argv[i], "--repeat")) repeat = std::stoi(argv[i + 1]);
  }
  const bool gpu = fg::linalg::cudaAvailable();
  std::vector<std::pair<std::string, std::unique_ptr<fg::linalg::Backend>>> backends;
  backends.emplace_back("cpu", fg::linalg::makeCpuBackend());
  if (gpu) backends.emplace_back("cuda", fg::linalg::makeCudaBackend(0));
  std::printf("FourierGrating %s - fgbench\n\n", FG_VERSION);
  for (const auto& b : backends) std::printf("- %s: %s\n", b.first.c_str(), b.second->name().c_str());
  if (!gpu) std::printf("- cuda: not available\n");
  std::printf("\n| n | backend | ZGEMM n^3 (s) | GFLOP/s | LU solve, n rhs (s) | eig (s) |\n|---|---|---|---|---|---|\n");
  for (int n : sizes) {
    std::srand(1);
    const fg::MatC A = fg::MatC::Random(n, n), B = fg::MatC::Random(n, n);
    for (const auto& b : backends) {
      const double tg = seconds([&] { (void)b.second->gemm(A, B); }, repeat);
      const double ts = seconds([&] { (void)b.second->solve(A, B); }, repeat);
      const double te = seconds([&] { (void)b.second->eig(A); }, repeat);
      std::printf("| %d | %s | %.3f | %.1f | %.3f | %.3f |\n", n, b.first.c_str(), tg, 8.0 * n * n * double(n) / tg / 1e9,
                  ts, te);
      std::fflush(stdout);
    }
  }
  std::printf("\nRCWA, Ag hole array (V10), lambda = 1.388 um, one solve:\n\n| N | 4 n_f | backend | time (s) |\n|---|---|---|---|\n");
  for (int N : orders) {
    for (const auto& b : backends) {
      fg::SolverOptions o;
      o.Nx = o.Nz = N;
      o.backend = b.first;
      const double t = seconds([&] { (void)fg::makeSolver(o)->solve(holeArray(), {1.388, 0.0, fg::PI, 0.0}); }, repeat);
      std::printf("| %d | %d | %s | %.2f |\n", N, 4 * (2 * N + 1) * (2 * N + 1), b.first.c_str(), t);
      std::fflush(stdout);
    }
  }
  return 0;
}
