// V15 (Plan section "M13: Tests und Abnahme"): pillars with circular cross-section. Contour2D with the
// exact circle against Polygon2D with L inscribed vertices: the efficiencies converge like O(L^-2)
// (expansion in even powers of 1/L), the Richardson extrapolation (two stages: L^-2, L^-4 from
// L = 128, 256, 512) agrees with Contour2D to 1e-8; RK = RCWA (y-invariant).
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <string>

#include "fg/io/config.hpp"

namespace {

const double kR = 0.14, kCx = 0.25, kCz = 0.25;

std::string base(const std::string& layer, const std::string& solver = "rcwa", const std::string& fac = "laurent") {
  return R"({"dimension": "2D", "solver": ")" + solver + R"(", "truncation": {"Nx": 4, "Nz": 4},
    "period_um": {"x": 0.5, "z": 0.5},
    "incidence": {"lambda_um": 0.633, "phi_deg": 12, "psi_deg": 25, "delta_deg": 40},
    "rk": {"tol": 1e-10},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "si": {"n": [3.0, 0.0]}},
    "superstrate": "air", "substrate": "glass",
    "layers": [)" + layer + R"(, "factorization": ")" + fac + R"("}]})";
}

std::string contourLayer() {
  char b[256];
  std::snprintf(b, sizeof b,
                R"({"type": "contour2d", "thickness_um": 0.2, "background": "air",
                    "regions": [{"material": "si", "circle": {"center_um": [%.17g, %.17g], "radius_um": %.17g}}])",
                kCx, kCz, kR);
  return b;
}

std::string polygonLayer(int L) {
  std::string v;
  for (int k = 0; k < L; ++k) {
    const double t = 2 * fg::PI * k / L;
    char b[96];
    std::snprintf(b, sizeof b, "%s[%.17g, %.17g]", k ? ", " : "", kCx + kR * std::cos(t), kCz + kR * std::sin(t));
    v += b;
  }
  return R"({"type": "polygon2d", "thickness_um": 0.2, "background": "air",
             "inclusions": [{"material": "si", "vertices_um": [)" + v + "]}]";
}

/// R_0, T_0, T_(-1,0) and the absorption-free energy balance.
std::vector<double> effs(const std::string& json) {
  const auto c = fg::io::parseConfig(json);
  const auto r = fg::makeSolver(c.solver)->solve(c.stack, c.scan.front());
  auto get = [](const std::vector<fg::OrderEfficiency>& v, int n, int m) {
    for (const auto& o : v)
      if (o.n == n && o.m == m) return o.efficiency;
    return 0.0;
  };
  return {get(r.reflected, 0, 0), get(r.transmitted, 0, 0), get(r.transmitted, -1, 0), r.diag.energyBalance};
}

}  // namespace

TEST(V15_CircularPillars, PolygonConvergesToContour) {
  const auto ref = effs(base(contourLayer()));
  EXPECT_NEAR(ref[3], 1.0, 1e-12);  // lossless
  std::vector<std::vector<double>> e;
  const int Ls[] = {16, 32, 64, 128, 256, 512};
  for (int L : Ls) e.push_back(effs(base(polygonLayer(L))));
  for (std::size_t q = 0; q < 3; ++q) {
    std::printf("order %zu: contour %.12f\n", q, ref[q]);
    for (std::size_t k = 0; k < e.size(); ++k) {
      const double d = e[k][q] - ref[q];
      std::printf("  L = %3d  diff %+.3e  L^2 diff %+.6f\n", Ls[k], d, d * Ls[k] * Ls[k]);
    }
    // O(L^-2): L^2 * diff approaches a constant (ratio of successive differences -> 4)
    const double r1 = (e[3][q] - ref[q]) / (e[4][q] - ref[q]), r2 = (e[4][q] - ref[q]) / (e[5][q] - ref[q]);
    EXPECT_NEAR(r1, 4.0, 0.05) << "order " << q;
    EXPECT_NEAR(r2, 4.0, 0.02) << "order " << q;
    // Richardson extrapolation: one stage (L^-2) leaves O(L^-4) ~ 1e-8 at L = 512, two stages remove it
    const double r256 = (4 * e[4][q] - e[3][q]) / 3, r512 = (4 * e[5][q] - e[4][q]) / 3;
    const double rich = (16 * r512 - r256) / 15;
    std::printf("  Richardson: one stage %+.3e, two stages %+.3e\n", r512 - ref[q], rich - ref[q]);
    EXPECT_NEAR(rich, ref[q], 1e-8) << "order " << q;
  }
}

TEST(V15_CircularPillars, RKEqualsRCWA) {
  const auto a = effs(base(contourLayer(), "rcwa", "FFF_sym"));
  const auto b = effs(base(contourLayer(), "rk", "FFF_sym"));
  for (std::size_t q = 0; q < 3; ++q) EXPECT_NEAR(a[q], b[q], 1e-8) << "order " << q;
  EXPECT_NEAR(a[3], 1.0, 1e-10);
}

TEST(V15_CircularPillars, FFFTrendAndPolygonRegion) {
  // FFF_sym: Polygon2D (NF from edges) approaches Contour2D (NF from curve samples)
  const auto ref = effs(base(contourLayer(), "rcwa", "FFF_sym"));
  const auto p64 = effs(base(polygonLayer(64), "rcwa", "FFF_sym"));
  const auto p256 = effs(base(polygonLayer(256), "rcwa", "FFF_sym"));
  for (std::size_t q = 0; q < 3; ++q) {
    std::printf("FFF_sym order %zu: contour %.10f, L=64 %+.2e, L=256 %+.2e\n", q, ref[q], p64[q] - ref[q], p256[q] - ref[q]);
    EXPECT_LT(std::abs(p256[q] - ref[q]), std::abs(p64[q] - ref[q]) + 1e-9);
    EXPECT_LT(std::abs(p256[q] - ref[q]), 1e-3);
  }
  // a contour2d polygon region is a Polygon2D (identical coefficients and normal field)
  const std::string poly = R"({"type": "contour2d", "thickness_um": 0.2, "background": "air",
      "regions": [{"material": "si", "polygon_um": [[0.1, 0.1], [0.4, 0.12], [0.35, 0.38], [0.12, 0.3]]}])";
  const std::string poly2 = R"({"type": "polygon2d", "thickness_um": 0.2, "background": "air",
      "inclusions": [{"material": "si", "vertices_um": [[0.1, 0.1], [0.4, 0.12], [0.35, 0.38], [0.12, 0.3]]}])";
  const auto c1 = effs(base(poly, "rcwa", "FFF_sym")), c2 = effs(base(poly2, "rcwa", "FFF_sym"));
  for (std::size_t q = 0; q < 3; ++q) EXPECT_NEAR(c1[q], c2[q], 1e-12) << "order " << q;
}
