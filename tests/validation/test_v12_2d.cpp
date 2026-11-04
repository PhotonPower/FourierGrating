// V12 (2D): near fields of crossed gratings (Plan section "Nahfeld", M9; docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "fg/physics/anisotropic.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::Mat3C;
using fg::VecC;

namespace {
constexpr double DEG = fg::PI / 180.0;
constexpr double LAM = 0.6328;

/// Sy summed over all orders (Fourier coefficients F = ([Ex], [Ez], [Hx], [Hz])).
double flux(const VecC& F) {
  const Eigen::Index n = F.size() / 4;
  double s = 0;
  for (Eigen::Index r = 0; r < n; ++r)
    s += 0.5 * std::real(F(n + r) * std::conj(F(2 * n + r)) - F(r) * std::conj(F(3 * n + r)));
  return s;
}

fg::SolverOptions opts(fg::SolverKind k, int Nx, int Nz) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = Nx;
  o.Nz = Nz;
  o.rcwaSlices = 8;
  o.rkTol = 1e-9;
  return o;
}

fg::LayerEntry entry(fg::LayerSpec s, fg::Factorization f = fg::Factorization::Auto) {
  fg::LayerEntry e;
  e.spec = std::move(s);
  e.factorization = f;
  return e;
}

fg::Stack base(double dx, double dz) {
  fg::Stack g;
  g.dx_um = dx;
  g.dz_um = dz;
  g.materials["air"] = fg::Material::fromIndex("air", 1.0);
  g.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
  g.materials["bar"] = fg::Material::fromIndex("bar", 2.0);
  g.superstrate = "air";
  g.substrate = "glass";
  return g;
}

fg::Polygon2DLayer pillar(double t, const std::string& m) {
  fg::Polygon2DLayer p;
  p.thickness_um = t;
  p.background = "air";
  p.inclusions.push_back({{{0.2, 0.25}, {0.6, 0.25}, {0.6, 0.75}, {0.2, 0.75}}, m});
  return p;
}

/// Crossed stack: pillars (y-invariant), homogeneous layer, egg crate (y-dependent).
fg::Stack crossedStack() {
  auto g = base(1.0, 1.0);
  fg::Profile2DLayer e;
  e.shape = fg::Profile2DLayer::Shape::EggCrate;
  e.depth_um = 0.2;
  e.below = "bar";
  e.above = "air";
  e.slices = 4;
  g.layers = {entry(pillar(0.3, "bar")), entry(fg::HomogeneousLayer{0.1, "bar"}), entry(e)};
  return g;
}
}  // namespace

// z-invariant structures: the 2D near field equals the 1D (conical) near field; horizontal maps
// are the 1D line times exp(i gamma0 (z - z0)).
TEST(V12_2D, ZInvariantEqualsOneD) {
  auto g1 = base(1.0, 0.0);
  g1.layers = {entry(fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"bar", "air"}})};
  fg::Profile1DLayer p1;
  p1.depth_um = 0.2;
  p1.below = "bar";
  p1.above = "air";
  g1.layers.push_back(entry(p1));
  auto g2 = base(1.0, 0.8);
  fg::Polygon2DLayer st;
  st.thickness_um = 0.4;
  st.background = "air";
  st.inclusions.push_back({{{0.0, 0.0}, {0.45, 0.0}, {0.45, 0.8}, {0.0, 0.8}}, "bar"});
  fg::Profile2DLayer p2;
  p2.shape = fg::Profile2DLayer::Shape::SinusX;
  p2.depth_um = 0.2;
  p2.below = "bar";
  p2.above = "air";
  g2.layers = {entry(st), entry(p2)};
  const fg::Wavelength wl{LAM, 20 * DEG, 30 * DEG, 50 * DEG};
  fg::NearFieldRequest rq;
  rq.nx = 50;
  rq.ny = 30;
  rq.nz = 4;
  rq.y_min_um = -0.1;
  rq.y_max_um = 0.8;
  rq.horizontal_y_um = 0.25;
  rq.vertical_z_um = 0.3;
  for (auto rc : {fg::Reconstruction::Direct, fg::Reconstruction::FFF})
    for (int Nz : {0, 1}) {
      rq.reconstruction = rc;
      const auto a = fg::makeSolver(opts(fg::SolverKind::RCWA, 8, 0))->nearField(g1, wl, rq);
      const auto b = fg::makeSolver(opts(fg::SolverKind::RCWA, 8, Nz))->nearField(g2, wl, rq);
      ASSERT_EQ(a.vertical.E[0].size(), b.vertical.E[0].size());
      double dv = 0, dh = 0;
      for (int c = 0; c < 3; ++c)
        for (std::size_t k = 0; k < a.vertical.E[c].size(); ++k)
          dv = std::max({dv, std::abs(a.vertical.E[c][k] - b.vertical.E[c][k]),
                         std::abs(a.vertical.H[c][k] - b.vertical.H[c][k])});
      ASSERT_EQ(b.horizontal.E[0].size(), 50u * 4u);
      EXPECT_EQ(b.horizontal.axis2Name, "z");
      const double g0 = 2 * fg::PI / LAM * std::sin(20 * DEG) * std::sin(30 * DEG);
      for (std::size_t j = 0; j < 4; ++j) {
        const cplx ph = std::exp(fg::I * (g0 * (b.horizontal.axis2_um[j] - rq.vertical_z_um)));
        for (int c = 0; c < 3; ++c)
          for (std::size_t k = 0; k < 50; ++k)
            dh = std::max({dh, std::abs(a.horizontal.E[c][k] * ph - b.horizontal.E[c][j * 50 + k]),
                           std::abs(a.horizontal.H[c][k] * ph - b.horizontal.H[c][j * 50 + k])});
      }
      EXPECT_LT(dv, 1e-9) << "Nz=" << Nz;
      EXPECT_LT(dh, 1e-9) << "Nz=" << Nz;
    }
}

// Tangential fields (F) continuous at all element and slice boundaries, RCWA and RK.
TEST(V12_2D, ContinuityAtAllBoundaries) {
  const auto g = crossedStack();
  std::vector<double> yy;
  const double d = 1e-9;
  for (double b : {0.0, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6}) {
    yy.push_back(b - d);
    yy.push_back(b + d);
  }
  for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
    const auto fs = fg::makeSolver(opts(kind, 2, 2))->fieldCoefficients(g, {LAM, 20 * DEG, 30 * DEG, 40 * DEG}, yy);
    for (std::size_t i = 0; i + 1 < fs.samples.size(); i += 2)
      EXPECT_LT((fs.samples[i].F - fs.samples[i + 1].F).norm() / fs.samples[i].F.norm(), 1e-6) << "y=" << yy[i];
  }
}

// Lossless crossed stack: Sy summed over the orders is the same at every height (= transmitted power).
TEST(V12_2D, FluxConstantInLosslessStack) {
  auto g = base(1.0, 1.0);
  g.layers = {entry(pillar(0.3, "bar")), entry(fg::HomogeneousLayer{0.1, "bar"})};
  const fg::Wavelength wl{LAM, 20 * DEG, 25 * DEG, 40 * DEG};
  const std::vector<double> ys{-0.3, 0.05, 0.2, 0.35, 0.45, 0.7};
  auto solver = fg::makeSolver(opts(fg::SolverKind::RCWA, 3, 3));
  const auto fs = solver->fieldCoefficients(g, wl, ys);
  const auto res = solver->solve(g, wl);
  double sumT = 0;
  for (const auto& o : res.transmitted) sumT += o.efficiency;
  const double f0 = flux(fs.samples[0].F);
  for (const auto& s : fs.samples) EXPECT_NEAR(flux(s.F), f0, 1e-12);
  EXPECT_NEAR(-f0 / (0.5 * std::cos(20 * DEG)), sumT, 1e-10);
}

// RK field = RCWA field inside a y-invariant 2D layer.
TEST(V12_2D, RKFieldEqualsRCWAField) {
  auto g = base(1.0, 1.0);
  g.layers = {entry(pillar(0.3, "bar"))};
  const std::vector<double> ys{0.05, 0.15, 0.29};
  const fg::Wavelength wl{LAM, 20 * DEG, 30 * DEG, 40 * DEG};
  const auto a = fg::makeSolver(opts(fg::SolverKind::RCWA, 2, 2))->fieldCoefficients(g, wl, ys);
  const auto b = fg::makeSolver(opts(fg::SolverKind::RK, 2, 2))->fieldCoefficients(g, wl, ys);
  for (std::size_t i = 0; i < ys.size(); ++i)
    EXPECT_LT((a.samples[i].F - b.samples[i].F).norm() / a.samples[i].F.norm(), 1e-7) << "y=" << ys[i];
}

// Horizontal map through a pillar: D_x = eps E_x across the wall x = 0.2 um (row z = 0.5 um).
// Isotropic (tensor NF, FFF_sym) and anisotropic (cobalt-like tensor, FFF_B) reconstruction.
TEST(V12_2D, FFFReconstructionKeepsNormalDContinuous) {
  for (bool aniso : {false, true}) {
    auto g = base(1.0, 1.0);
    Mat3C eps = 4.0 * Mat3C::Identity();
    if (aniso) {
      eps(0, 2) = cplx(0.3, 0.1);
      eps(2, 0) = -eps(0, 2);
      eps(1, 1) = 3.0;
      g.materials["bar"] = fg::Material::fromTensor("bar", eps);
    }
    g.layers = {entry(pillar(0.3, "bar"))};
    fg::NearFieldRequest rq;
    rq.horizontal_y_um = 0.15;
    rq.ny = 2;
    rq.nx = 1000;
    rq.nz = 2;  // z = 0 and z = 0.5 um
    double jump[2];
    int k = 0;
    for (auto rc : {fg::Reconstruction::Direct, fg::Reconstruction::FFF}) {
      rq.reconstruction = rc;
      const auto nf = fg::makeSolver(opts(fg::SolverKind::RCWA, 6, 6))->nearField(g, {LAM, 20 * DEG, 0, 90 * DEG}, rq);
      const auto& h = nf.horizontal;
      ASSERT_EQ(h.E[0].size(), 2000u);
      const std::size_t row = 1000, iL = row + 199, iR = row + 201;  // x = 0.199 (air), 0.201 (pillar)
      const cplx DL = h.E[0][iL];
      const cplx DR = eps(0, 0) * h.E[0][iR] + eps(0, 1) * h.E[1][iR] + eps(0, 2) * h.E[2][iR];
      jump[k++] = std::abs(DL - DR) / std::abs(DL);
      EXPECT_NEAR(h.material[iR], 4.0, 1e-12);
      EXPECT_NEAR(h.material[iL], 1.0, 1e-12);
    }
    EXPECT_GT(jump[0], 1.0) << "aniso=" << aniso;   // direct: 293 %; FFF: 5 % (iso), 3 % (aniso)
    EXPECT_LT(jump[1], 0.1) << "aniso=" << aniso;   // FFF-consistent: D_x continuous
  }
}

// Shapes of the 2D cuts and the written files (meta: axis names).
TEST(V12_2D, CutsAndFiles) {
  const auto g = crossedStack();
  fg::NearFieldRequest rq;
  rq.nx = 24;
  rq.ny = 10;
  rq.nz = 16;
  rq.periods = 2;
  rq.horizontal_y_um = 0.5;
  const auto nf = fg::makeSolver(opts(fg::SolverKind::RCWA, 2, 2))->nearField(g, {LAM, 10 * DEG, 20 * DEG, 90 * DEG}, rq);
  EXPECT_EQ(nf.vertical.E[0].size(), 24u * 10u);
  EXPECT_EQ(nf.horizontal.E[2].size(), 24u * 16u);
  EXPECT_EQ(nf.horizontal.axis2_um.size(), 16u);
  EXPECT_NEAR(nf.horizontal.axis2_um.back(), 2.0 * 15 / 16, 1e-12);
  // horizontal map row z = 0 equals the vertical cut (z0 = 0) at the same height when y0 is a sample
  const auto dir = std::filesystem::temp_directory_path() / "fg_nf_test_2d";
  fg::NearFieldWriter::write(nf.horizontal, dir.string());
  std::ifstream in(dir / "meta.json");
  const auto meta = nlohmann::json::parse(in);
  EXPECT_EQ(meta["cut"], "horizontal");
  EXPECT_EQ(meta["axes"]["axis2"]["name"], "z");
  EXPECT_EQ(meta["axes"]["axis2"]["values_um"].size(), 16u);
  EXPECT_TRUE(std::filesystem::exists(dir / "material.npy"));
}
