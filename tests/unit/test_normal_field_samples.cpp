// Algorithm NF from boundary samples (Plan section "M13: Normalenfeld mit Randstichproben").
#include <gtest/gtest.h>

#include <chrono>
#include <cmath>

#include "fg/geometry/normal_field.hpp"

namespace {

fg::NormalFieldRaster sampled(const fg::Polygon2DLayer& l, double dx, double dz, int Gx, int Gz) {
  std::vector<fg::BoundaryChain> chains;
  for (const auto& inc : l.inclusions) {
    auto c = fg::polygonChains(inc.vertices);
    chains.insert(chains.end(), c.begin(), c.end());
  }
  return fg::normalFieldRaster(
      chains, [&](double x, double z) { return fg::polygonRegion(l, dx, dz, x, z); },
      [&](int k) { return k < 0 ? l.background : l.inclusions[static_cast<std::size_t>(k)].material; }, l.normalField,
      dx, dz, Gx, Gz);
}

void expectBitIdentical(const fg::Polygon2DLayer& l, double dx, double dz, int Gx, int Gz) {
  const auto a = fg::normalFieldRaster(l, dx, dz, Gx, Gz);
  const auto b = sampled(l, dx, dz, Gx, Gz);
  EXPECT_FALSE(b.has3D());
  EXPECT_EQ(a.region, b.region);
  EXPECT_EQ(a.nxx, b.nxx);
  EXPECT_EQ(a.nxz, b.nxz);
  EXPECT_EQ(a.nzz, b.nzz);
  EXPECT_EQ(a.nx, b.nx);
  EXPECT_EQ(a.nz, b.nz);
}

fg::Polygon2DLayer layout() {
  fg::Polygon2DLayer l;
  l.background = "air";
  // square on the grid (medial-axis ties), triangle across the cell corner, two touching boxes of the
  // same material (internal edge removed) and one of another material
  l.inclusions.push_back({{{0.1, 0.1}, {0.3, 0.1}, {0.3, 0.3}, {0.1, 0.3}}, "si"});
  l.inclusions.push_back({{{0.9, 0.85}, {1.15, 0.9}, {0.95, 1.1}}, "si"});
  l.inclusions.push_back({{{0.5, 0.5}, {0.7, 0.5}, {0.7, 0.6}, {0.5, 0.6}}, "au"});
  l.inclusions.push_back({{{0.7, 0.5}, {0.8, 0.5}, {0.8, 0.6}, {0.7, 0.6}}, "au"});
  l.inclusions.push_back({{{0.5, 0.6}, {0.7, 0.6}, {0.6, 0.75}}, "sio2"});
  return l;
}

}  // namespace

TEST(NormalFieldSamples, PolygonsBitIdenticalTensor) {
  fg::Polygon2DLayer l = layout();
  expectBitIdentical(l, 1.0, 0.9, 64, 64);
  expectBitIdentical(l, 1.0, 0.9, 128, 96);
}

TEST(NormalFieldSamples, PolygonsBitIdenticalVector) {
  fg::Polygon2DLayer l = layout();
  l.normalField.method = fg::NormalFieldOptions::Method::Vector;
  expectBitIdentical(l, 1.0, 0.9, 128, 128);
}

TEST(NormalFieldSamples, ThreeDimensionalProjector) {
  // circle of radius 0.2 with tilted normals (frustum wall): n = (cos t, -0.4, sin t) / |.|
  const double R = 0.2, cx = 0.5, cz = 0.5, ny = -0.4;
  fg::BoundaryChain c;
  const int S = 400;
  for (int s = 0; s < S; ++s) {
    const double t = 2 * fg::PI * s / S;
    const fg::Vec3D n = fg::Vec3D(std::cos(t), ny, std::sin(t)).normalized();
    c.samples.push_back({cx + R * std::cos(t), cz + R * std::sin(t), n, false});
  }
  fg::NormalFieldOptions opt;
  const auto inside = [&](double x, double z) { return std::hypot(x - cx, z - cz) < R ? 0 : -1; };
  const auto R3 = fg::normalFieldRaster({c}, inside, [](int k) { return std::string(k < 0 ? "air" : "si"); }, opt, 1.0,
                                        1.0, 64, 64);
  ASSERT_TRUE(R3.has3D());
  double maxTrace = 0, nyMin = 0;
  for (std::size_t j = 0; j < R3.nxx.size(); ++j) {
    maxTrace = std::max(maxTrace, std::abs(R3.nxx[j] + R3.nyy[j] + R3.nzz[j] - 1.0));
    nyMin = std::min(nyMin, R3.nyy[j]);
    EXPECT_NEAR(R3.nx[j] * R3.nx[j] + R3.ny[j] * R3.ny[j] + R3.nz[j] * R3.nz[j], 1.0, 1e-12);
  }
  EXPECT_LT(maxTrace, 1e-12);  // smoothing is linear: trace 1 everywhere
  EXPECT_GE(nyMin, 0.0);
  // close to the wall (radius 0.25, unsmoothed region) n_y^2 = 0.16 / 1.16
  const std::size_t j = static_cast<std::size_t>(48) * 64 + 32;  // (x, z) = (0.75, 0.5)
  EXPECT_NEAR(R3.nyy[j], 0.16 / 1.16, 2e-3);
  // coefficients: all six components
  const auto nf = fg::rasterCoefficients(R3, fg::IndexMap(3, 3));
  EXPECT_NEAR(nf.NN[1][1](0, 0).real(), 0.16 / 1.16, 2e-2);
  EXPECT_EQ(nf.NN[0][1](1, 1), nf.NN[1][0](1, 1));
}

TEST(NormalFieldSamples, Performance) {
  // Plan: 512^2 raster with 1e4 samples in under 1 s on one core (Release)
  fg::BoundaryChain c;
  const int S = 10000;
  for (int s = 0; s < S; ++s) {
    const double t = 2 * fg::PI * s / S;
    c.samples.push_back({0.5 + 0.3 * std::cos(t), 0.5 + 0.2 * std::sin(t), fg::Vec3D(std::cos(t), 0, std::sin(t)), false});
  }
  const auto inside = [](double x, double z) {
    return std::pow((x - 0.5) / 0.3, 2) + std::pow((z - 0.5) / 0.2, 2) < 1 ? 0 : -1;
  };
  const auto t0 = std::chrono::steady_clock::now();
  const auto R = fg::normalFieldRaster({c}, inside, [](int k) { return std::string(k < 0 ? "a" : "b"); },
                                       fg::NormalFieldOptions{}, 1.0, 1.0, 512, 512);
  const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  RecordProperty("seconds", std::to_string(s));
  std::printf("NF 512^2, 1e4 samples: %.3f s\n", s);
  EXPECT_EQ(R.nxx.size(), 512u * 512u);
#ifdef NDEBUG
  // measured (one core): GCC 0.47 s, MSVC 0.68 s; 2 s leaves room for slower CI machines
  EXPECT_LT(s, 2.0);
#endif
}
