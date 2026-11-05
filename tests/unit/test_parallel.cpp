// Thread model and parallel assembly (Plan section "M15: Profiling und Parallelisierung"): results
// do not depend on the thread count (bit-identical), tree = sequential to 1e-12. Small sizes: this
// test also runs under ThreadSanitizer (CI job "tsan", std::thread instead of OpenMP).
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fg/core/parallel.hpp"
#include "fg/geometry/section.hpp"
#include "fg/io/config.hpp"

namespace {

/// Frustum pillar (tilted walls, 3D normals) as in V16a.
class Frustum final : public fg::SectionProvider {
 public:
  fg::CrossSection section(double y) const override {
    const double R = 0.15 - 0.08 * y / 0.2;
    fg::CrossSection cs;
    fg::geometry::Region r;
    r.material = "pillar";
    r.outer = fg::geometry::Contour2D::circle({0.25, 0.25}, R);
    cs.regions.push_back(std::move(r));
    fg::BoundaryChain c;
    for (int s = 0; s < 180; ++s) {
      const double t = 2 * fg::PI * s / 180;
      c.samples.push_back({0.25 + R * std::cos(t), 0.25 + R * std::sin(t),
                           fg::Vec3D(std::cos(t), 0.4, std::sin(t)).normalized(), false});
    }
    cs.chains.push_back(std::move(c));
    return cs;
  }
  double thickness_um() const override { return 0.2; }
};

fg::io::RunConfig config(const std::string& solver, const std::string& layers,
                         const std::string& dim = "2D") {
  return fg::io::parseConfig(R"({"dimension": ")" + dim + R"(", "solver": ")" + solver +
                             R"(", "truncation": {"Nx": 2, "Nz": )" + (dim == "2D" ? "2" : "0") +
                             R"(},
    "period_um": {"x": 0.5, "z": )" +
                             (dim == "2D" ? "0.5" : "0") + R"(},
    "incidence": {"lambda_um": 0.633, "phi_deg": 15, "psi_deg": 30, "delta_deg": 35},
    "rk": {"tol": 1e-6}, "rcwa": {"slices": 9},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "pillar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass", "layers": [)" +
                             layers + "]}");
}

void addSectioned(fg::io::RunConfig& c, int chebyshevNodes,
                  fg::Factorization fac = fg::Factorization::Auto) {
  fg::SectionedLayer l;
  l.provider = std::make_shared<Frustum>();
  l.background = "air";
  l.normalField.grid = 64;
  l.chebyshevNodes = chebyshevNodes;
  fg::LayerEntry le;
  le.spec = l;
  le.factorization = fac;
  c.stack.layers.insert(c.stack.layers.begin(), std::move(le));
}

std::vector<double> efficiencies(const fg::Result& r) {
  std::vector<double> v;
  for (const auto* side : {&r.reflected, &r.transmitted})
    for (const auto& o : *side) v.push_back(o.efficiency);
  v.push_back(r.absorption);
  return v;
}

/// Efficiencies with `threads` slice threads (fresh layer models each time).
std::vector<double> run(fg::io::RunConfig c, int threads, bool tree) {
  c.solver.threadsSlices = threads;
  c.solver.treeReduction = tree;
  return efficiencies(fg::makeSolver(c.solver)->solve(c.stack, c.scan.front()));
}

void expectThreadInvariant(const fg::io::RunConfig& c) {
  const auto ref = run(c, 1, false);
  for (int t : {2, 4}) EXPECT_EQ(run(c, t, false), ref) << "sequential chain, " << t << " threads";
  const auto tree = run(c, 1, true);
  for (int t : {2, 4}) EXPECT_EQ(run(c, t, true), tree) << "tree, " << t << " threads";
  ASSERT_EQ(tree.size(), ref.size());
  for (std::size_t k = 0; k < ref.size(); ++k) EXPECT_NEAR(tree[k], ref[k], 1e-12);
}

}  // namespace

TEST(ParallelFor, ResultsByIndexAndFirstExceptionRethrown) {
  std::vector<int> v(100, 0);
  fg::parallelFor(100, 4, [&](int i) { v[static_cast<std::size_t>(i)] = i * i; });
  for (int i = 0; i < 100; ++i) EXPECT_EQ(v[static_cast<std::size_t>(i)], i * i);
  std::atomic<int> done{0};
  try {
    fg::parallelFor(20, 4, [&](int i) {
      ++done;
      if (i == 7 || i == 13) throw std::runtime_error("task " + std::to_string(i));
    });
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
    EXPECT_STREQ(e.what(), "task 7");  // lowest index, after all tasks finished
  }
  EXPECT_EQ(done.load(), 20);
}

TEST(ParallelFor, InnerLoopsDoNotNest) {
  EXPECT_EQ(fg::innerThreads(), fg::availableThreads());
  std::vector<int> inner(8, -1);
  fg::parallelFor(8, 4, [&](int i) { inner[static_cast<std::size_t>(i)] = fg::innerThreads(); });
  for (int t : inner) EXPECT_EQ(t, 1);
}

TEST(ParallelAssembly, RcwaProfileAndSectionedLayer) {
  auto c = config("rcwa",
                  R"({"type": "profile2d", "shape": "egg_crate", "depth_um": 0.1, "below": "pillar",
                              "above": "air", "slices": 7},
                             {"type": "homogeneous", "thickness_um": 0.05, "material": "glass"})");
  addSectioned(c, 0);
  expectThreadInvariant(c);
}

TEST(ParallelAssembly, SectionedLayerWithChebyshevInterpolant) {
  // the interpolant is built once by whichever slice task needs it first (its nodes in parallel);
  // Laurent: no normal field per node (fast enough for ThreadSanitizer)
  auto c = config("rcwa", "");
  addSectioned(c, 4, fg::Factorization::Laurent);
  EXPECT_EQ(run(c, 4, false), run(c, 1, false));
}

TEST(ParallelAssembly, RungeKuttaSlices) {
  auto c =
      config("rk", R"({"type": "profile1d", "shape": "sinus", "depth_um": 0.3, "below": "pillar",
                            "above": "air"})",
             "1D");
  expectThreadInvariant(c);
}

TEST(ParallelAssembly, NearFieldRecordsAndScan) {
  auto c = config("rcwa",
                  R"({"type": "lamellar1d", "thickness_um": 0.2, "breakpoints_um": [0.0, 0.2, 0.5],
                              "materials": ["pillar", "air"]},
                             {"type": "profile2d", "shape": "egg_crate", "depth_um": 0.1, "below": "pillar",
                              "above": "air", "slices": 5})");
  std::vector<std::vector<fg::cplx>> fields;
  for (int t : {1, 3}) {
    c.solver.threadsSlices = t;
    const auto fs = fg::makeSolver(c.solver)->fieldCoefficients(c.stack, c.scan.front(),
                                                                {-0.05, 0.03, 0.12, 0.27, 0.4});
    std::vector<fg::cplx> v;
    for (const auto& s : fs.samples) v.insert(v.end(), s.F.data(), s.F.data() + s.F.size());
    fields.push_back(v);
  }
  EXPECT_EQ(fields[0], fields[1]);
  // scan: points in parallel (2 x 2 threads) = sequential
  std::vector<fg::Wavelength> scan;
  for (double l : {0.6, 0.633, 0.7}) scan.push_back({l, 0.2, 0.5, 0.6});
  c.solver.threadsScan = 1;
  c.solver.threadsSlices = 1;
  const auto seq = fg::solveScan(c.solver, c.stack, scan);
  c.solver.threadsScan = 2;
  c.solver.threadsSlices = 2;
  const auto par = fg::solveScan(c.solver, c.stack, scan);
  for (std::size_t i = 0; i < scan.size(); ++i)
    EXPECT_EQ(efficiencies(par[i]), efficiencies(seq[i]));
}
