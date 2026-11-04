// Algorithm NF for y-invariant 2D cross-sections (Plan alg:nfield) and Profile2D coefficients.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "../../src/solvers/layer_model.hpp"
#include "fg/geometry/normal_field.hpp"
#include "fg/geometry/stack.hpp"

using fg::cplx;
using Method = fg::NormalFieldOptions::Method;

namespace {

fg::Polygon2DLayer strip(double x0, double x1, double dz, Method m = Method::Tensor) {
  fg::Polygon2DLayer l;
  l.thickness_um = 0.1;
  l.background = "a";
  l.inclusions.push_back({{{x0, 0.0}, {x1, 0.0}, {x1, dz}, {x0, dz}}, "b"});
  l.normalField.method = m;
  return l;
}

fg::Polygon2DLayer square(double cx, double cz, double w, int grid = 256) {
  fg::Polygon2DLayer l;
  l.thickness_um = 0.1;
  l.background = "a";
  const double h = 0.5 * w;
  l.inclusions.push_back({{{cx - h, cz - h}, {cx + h, cz - h}, {cx + h, cz + h}, {cx - h, cz + h}}, "b"});
  l.normalField.grid = grid;
  return l;
}

double maxDev(const fg::fourier::CoefficientGrid& g, cplx c00) {
  double m = 0;
  for (int p = -g.Px(); p <= g.Px(); ++p)
    for (int q = -g.Pz(); q <= g.Pz(); ++q) m = std::max(m, std::abs(g(p, q) - (p == 0 && q == 0 ? c00 : 0.0)));
  return m;
}

}  // namespace

// z-invariant lamella: the tensor method gives exactly N N^T = x x^T (prerequisite of V9).
TEST(NormalField2D, LamellaTensorMethodIsExact) {
  const fg::IndexMap map(4, 3);
  const auto nf = fg::normalField(strip(0.2, 0.6, 0.8), 1.0, 0.8, map);
  EXPECT_LT(maxDev(nf.NN[0][0], 1.0), 1e-14);
  EXPECT_LT(maxDev(nf.NN[0][2], 0.0), 1e-14);
  EXPECT_LT(maxDev(nf.NN[2][2], 0.0), 1e-14);
  EXPECT_LT(maxDev(nf.NN[1][1], 0.0), 1e-15);
  EXPECT_FALSE(nf.hasSigned);
}

// Same lamella split into two touching rectangles: the internal edge is no material boundary.
TEST(NormalField2D, InternalEdgesAreRemoved) {
  fg::Polygon2DLayer l = strip(0.2, 0.6, 0.8);
  l.inclusions = {{{{0.2, 0.0}, {0.6, 0.0}, {0.6, 0.35}, {0.2, 0.35}}, "b"},
                  {{{0.2, 0.35}, {0.6, 0.35}, {0.6, 0.8}, {0.2, 0.8}}, "b"}};
  const auto nf = fg::normalField(l, 1.0, 0.8, fg::IndexMap(4, 3));
  EXPECT_LT(maxDev(nf.NN[0][0], 1.0), 1e-14);
  EXPECT_LT(maxDev(nf.NN[2][2], 0.0), 1e-14);
}

// Befund M7: the vector formulation of Plan v1.4 falls back to N_x^2 = 1/2 on medial lines that
// hit raster points (here x = 0.4 and x = 0.9 with G = 520), so a lamella is not reproduced exactly.
TEST(NormalField2D, VectorMethodFallbackOnMedialAxis) {
  fg::Polygon2DLayer l = strip(0.2, 0.6, 0.8, Method::Vector);
  l.normalField.grid = 520;
  const auto R = fg::normalFieldRaster(l, 1.0, 0.8, 520, 520);
  double mn = 1;
  for (double v : R.nxx) mn = std::min(mn, v);
  EXPECT_NEAR(mn, 0.5, 1e-15);
  const auto nf = fg::normalField(l, 1.0, 0.8, fg::IndexMap(4, 3));
  EXPECT_GT(maxDev(nf.NN[0][0], 1.0), 1e-4);
  // the tensor method on the same raster is exact
  l.normalField.method = Method::Tensor;
  EXPECT_LT(maxDev(fg::normalField(l, 1.0, 0.8, fg::IndexMap(4, 3)).NN[0][0], 1.0), 1e-14);
}

// Square hole: x <-> z symmetry, trace 1, N N^T = n n^T in the middle of the edges.
TEST(NormalField2D, SquareSymmetryTraceAndEdges) {
  const int G = 256;
  const auto l = square(0.45, 0.45, 0.25, G);
  const auto R = fg::normalFieldRaster(l, 0.9, 0.9, G, G);
  double asym = 0, tr = 0;
  for (int i = 0; i < G; ++i)
    for (int k = 0; k < G; ++k) {
      const auto a = static_cast<std::size_t>(i * G + k), b = static_cast<std::size_t>(k * G + i);
      asym = std::max(asym, std::abs(R.nxx[a] - R.nzz[b]) + std::abs(R.nxz[a] - R.nxz[b]));
      tr = std::max(tr, std::abs(R.nxx[a] + R.nzz[a] - 1.0));
    }
  EXPECT_LT(asym, 1e-12);
  EXPECT_LT(tr, 1e-12);
  // raster point next to the middle of the left edge (x = 0.325, z = 0.45)
  const int i = static_cast<int>(std::lround(0.325 / 0.9 * G)), k = G / 2;
  EXPECT_NEAR(R.nxx[static_cast<std::size_t>(i * G + k)], 1.0, 1e-12);
  EXPECT_EQ(R.region[static_cast<std::size_t>(k * G + k)], 0);  // centre is inside
  EXPECT_EQ(R.region[0], -1);
}

// Periodicity: a square across the cell corner equals the centred square shifted by half a period.
TEST(NormalField2D, PolygonAcrossCellBoundary) {
  const fg::IndexMap map(3, 3);
  const auto a = fg::normalField(square(0.0, 0.0, 0.3), 0.9, 0.9, map);
  const auto b = fg::normalField(square(0.45, 0.45, 0.3), 0.9, 0.9, map);
  double dev = 0;
  for (int i : {0, 2})
    for (int j : {0, 2})
      for (int p = -6; p <= 6; ++p)
        for (int q = -6; q <= 6; ++q)
          dev = std::max(dev, std::abs(a.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)](p, q) *
                                           ((p + q) % 2 == 0 ? 1.0 : -1.0) -
                                       b.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)](p, q)));
  EXPECT_LT(dev, 1e-10);
  EXPECT_EQ(fg::polygonRegion(square(0.0, 0.0, 0.3), 0.9, 0.9, 0.85, 0.88), 0);
  EXPECT_EQ(fg::polygonRegion(square(0.0, 0.0, 0.3), 0.9, 0.9, 0.45, 0.45), -1);
}

// ------------------------------------------------------------------ layer models (2D)
namespace {
fg::Stack stack2D(double dx, double dz) {
  fg::Stack st;
  st.dx_um = dx;
  st.dz_um = dz;
  st.materials["a"] = fg::Material::fromIndex("a", 1.0);
  st.materials["b"] = fg::Material::fromIndex("b", cplx(1.5, 0.1));
  return st;
}
/// Coefficient f_pq from a Laurent Toeplitz block: T(r(p,q), r(0,0)) = f_pq.
cplx coef(const fg::MatC& T, const fg::IndexMap& map, int p, int q) { return T(map.index(p, q), map.center()); }
}  // namespace

TEST(Polygon2DModel, ExactCoefficientsAndLaurent) {
  auto st = stack2D(0.9, 0.7);
  fg::LayerEntry le;
  auto l = square(0.3, 0.2, 0.2);
  std::reverse(l.inclusions[0].vertices.begin(), l.inclusions[0].vertices.end());  // clockwise input
  le.spec = l;
  le.factorization = fg::Factorization::Laurent;
  st.layers.push_back(le);
  const fg::IndexMap map(2, 2);
  const auto models = fg::detail::buildLayerModels(st, map, 10);
  const auto Q = models[0]->Q(0.0);
  const cplx e = cplx(1.5, 0.1) * cplx(1.5, 0.1);
  for (int p = -2; p <= 2; ++p)
    for (int q = -2; q <= 2; ++q) {
      const cplx ref = (p == 0 && q == 0 ? 1.0 : 0.0) +
                       (e - 1.0) * fg::fourier::rectangleCoefficient(0.2, 0.4, 0.1, 0.3, 0.9, 0.7, p, q);
      EXPECT_NEAR(std::abs(coef(Q(0, 0), map, p, q) - ref), 0.0, 1e-14);
    }
}

// Profile2D sinus_x reproduces the 1D sinus profile model (V9 prerequisite).
TEST(Profile2DModel, SinusXEqualsProfile1D) {
  auto st = stack2D(0.8, 0.6);
  fg::Profile1DLayer p1;
  p1.depth_um = 0.3;
  p1.x0_um = 0.1;
  p1.below = "b";
  p1.above = "a";
  fg::Profile2DLayer p2;
  p2.shape = fg::Profile2DLayer::Shape::SinusX;
  p2.depth_um = 0.3;
  p2.x0_um = 0.1;
  p2.below = "b";
  p2.above = "a";
  fg::LayerEntry e1, e2;
  e1.spec = p1;
  e2.spec = p2;
  st.layers = {e1, e2};
  const fg::IndexMap map(3, 2);
  const auto models = fg::detail::buildLayerModels(st, map, 10);
  for (double y : {0.02, 0.15, 0.27}) {
    const auto Qa = models[0]->Q(y), Qb = models[1]->Q(y);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) EXPECT_EQ((Qa(i, j) - Qb(i, j)).cwiseAbs().maxCoeff(), 0.0);
  }
}

// Egg crate: eps coefficients (exact in x, trapezoid in z) against brute-force quadrature.
TEST(Profile2DModel, EggCrateCoefficients) {
  auto st = stack2D(1.0, 0.8);
  fg::Profile2DLayer p;
  p.shape = fg::Profile2DLayer::Shape::EggCrate;
  p.depth_um = 0.4;
  p.below = "b";
  p.above = "a";
  fg::LayerEntry le;
  le.spec = p;
  le.factorization = fg::Factorization::Laurent;
  st.layers = {le};
  const fg::IndexMap map(2, 2);
  const auto models = fg::detail::buildLayerModels(st, map, 10);
  const double y = 0.13;
  const auto Q = models[0]->Q(y);
  const cplx e = cplx(1.5, 0.1) * cplx(1.5, 0.1);
  const int M = 200000;  // reference: exact along x (|x| < u), midpoint rule along z
  for (int pp : {0, 1, 2})
    for (int qq : {0, 2, -1}) {
      cplx s = 0;
      for (int k = 0; k < M; ++k) {
        const double z = 0.8 * (k + 0.5) / M, A = 0.2 * (1 + std::cos(2 * fg::PI * z / 0.8));
        if (A <= y) continue;
        const double u = std::acos(2 * y / A - 1) / (2 * fg::PI);
        const double cx = pp == 0 ? 2 * u : std::sin(2 * fg::PI * pp * u) / (fg::PI * pp);
        s += cx * std::exp(-fg::I * (2 * fg::PI * qq * z / 0.8));
      }
      s /= double(M);
      const cplx ref = (pp == 0 && qq == 0 ? 1.0 : 0.0) + (e - 1.0) * s;
      EXPECT_NEAR(std::abs(coef(Q(0, 0), map, pp, qq) - ref), 0.0, 1e-7) << pp << "," << qq;
    }
  // signed normal field (Plan eq:N2d): mean of N_y = 1/sqrt(1 + |grad g|^2)
  EXPECT_TRUE(models[0]->NN().hasSigned);
  double ny = 0;
  const int G = 400;
  for (int i = 0; i < G; ++i)
    for (int k = 0; k < G; ++k) {
      const double x = (i + 0.5) / G, z = 0.8 * (k + 0.5) / G, K = 2 * fg::PI, Kz = K / 0.8;
      const double gx = -0.1 * K * std::sin(K * x) * (1 + std::cos(Kz * z));
      const double gz = -0.1 * Kz * std::sin(Kz * z) * (1 + std::cos(K * x));
      ny += 1.0 / std::sqrt(1 + gx * gx + gz * gz);
    }
  EXPECT_NEAR(models[0]->NN().N[1](0, 0).real(), ny / (G * G), 1e-10);
}
