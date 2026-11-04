// V12: near fields (Plan section "Nahfeld", docs/validation.md).
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::VecC;

namespace {
constexpr double DEG = fg::PI / 180.0;
constexpr double LAM = 0.6328;

double flux(const VecC& F) {
  const Eigen::Index n = F.size() / 4;
  double s = 0;
  for (Eigen::Index r = 0; r < n; ++r)
    s += 0.5 * std::real(F(n + r) * std::conj(F(2 * n + r)) - F(r) * std::conj(F(3 * n + r)));
  return s;
}

fg::SolverOptions opts(fg::SolverKind k, int N) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rcwaSlices = 20;
  o.rkTol = 1e-10;
  return o;
}

fg::Stack lamellarStack(bool withProfile) {
  fg::Stack g;
  g.dx_um = 1.0;
  g.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  g.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
  g.materials["bar"] = fg::Material::fromIndex("bar", 2.0);
  fg::LayerEntry l1, l2;
  l1.spec = fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"bar", "superstrate"}};
  l2.spec = fg::HomogeneousLayer{0.1, "bar"};
  g.layers = {l1, l2};
  if (withProfile) {
    fg::Profile1DLayer p;
    p.depth_um = 0.3;
    p.below = "bar";
    p.above = "superstrate";
    fg::LayerEntry l3;
    l3.spec = p;
    g.layers.push_back(l3);
  }
  return g;
}
}  // namespace

TEST(V12_NearField, AiryLayerAnalyticTE) {
  const cplx n1(2.1, 0.05);
  const double t = 0.37, n2 = 1.5, th = 35 * DEG, k0 = 2 * fg::PI / LAM;
  fg::Stack st;
  st.dx_um = 0.7;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("g", n2);
  st.materials["L"] = fg::Material::fromIndex("L", n1);
  fg::LayerEntry le;
  le.spec = fg::HomogeneousLayer{t, "L"};
  st.layers.push_back(le);
  const std::vector<double> ys{-0.2, -0.05, 0.0, 0.1, 0.2, 0.37, 0.5, 0.8};
  // analytic: Ez = tau exp(-i b2 y) below, transfer matrix of (Ez, Hx) in the layer, A/B above
  const double a = std::sin(th);
  const cplx b0 = std::sqrt(cplx(1 - a * a)), b1 = std::sqrt(n1 * n1 - a * a), b2 = std::sqrt(cplx(n2 * n2 - a * a));
  auto prop = [&](double yn, cplx& Ez, cplx& Hx) {
    const cplx E0 = 1.0, H0 = -b2;
    Ez = std::cos(b1 * yn) * E0 + fg::I * std::sin(b1 * yn) / b1 * H0;
    Hx = fg::I * b1 * std::sin(b1 * yn) * E0 + std::cos(b1 * yn) * H0;
  };
  cplx Et, Ht;
  prop(t * k0, Et, Ht);
  const cplx A = (b0 * Et - Ht) / (2.0 * b0), B = (b0 * Et + Ht) / (2.0 * b0);
  const cplx sc = -1.0 / A;  // incident E_z = A_z = -1 for delta = 0, psi = 0
  for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
    const auto fs = fg::makeSolver(opts(kind, 2))->fieldCoefficients(st, {LAM, th, 0, 0}, ys);
    const int n = static_cast<int>(fs.alpha.size()), c = n / 2;
    for (std::size_t i = 0; i < ys.size(); ++i) {
      const double y = ys[i] * k0;
      cplx ref;
      if (ys[i] < 0)
        ref = sc * std::exp(-fg::I * b2 * y);
      else if (ys[i] <= t) {
        cplx Ez, Hx;
        prop(y, Ez, Hx);
        ref = sc * Ez;
      } else {
        const double d = y - t * k0;
        ref = sc * (A * std::exp(-fg::I * b0 * d) + B * std::exp(fg::I * b0 * d));
      }
      EXPECT_NEAR(std::abs(fs.samples[i].F(n + c) - ref), 0.0, 1e-12) << "y=" << ys[i];
    }
  }
}

TEST(V12_NearField, ContinuityAtAllSliceBoundaries) {
  const auto g = lamellarStack(true);
  std::vector<double> yy;
  const double d = 1e-9;
  for (double b : {0.0, 0.4, 0.5, 0.8}) {
    yy.push_back(b - d);
    yy.push_back(b + d);
  }
  for (int k = 1; k < 20; ++k) {  // RCWA staircase boundaries in the profile layer
    yy.push_back(0.5 + 0.3 * k / 20 - d);
    yy.push_back(0.5 + 0.3 * k / 20 + d);
  }
  for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK})
    for (double delta : {0.0, 90.0}) {
      const auto fs = fg::makeSolver(opts(kind, 10))->fieldCoefficients(g, {LAM, 20 * DEG, 0, delta * DEG}, yy);
      for (std::size_t i = 0; i + 1 < fs.samples.size(); i += 2)
        EXPECT_LT((fs.samples[i].F - fs.samples[i + 1].F).norm() / fs.samples[i].F.norm(), 1e-6)
            << "y=" << yy[i] << " delta=" << delta;
    }
}

TEST(V12_NearField, FluxConstantInLosslessLamellarStack) {
  const auto g = lamellarStack(false);
  const fg::Wavelength wl{LAM, 20 * DEG, 25 * DEG, 40 * DEG};
  const std::vector<double> ys{-0.3, 0.05, 0.2, 0.39, 0.45, 0.7};
  for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
    auto solver = fg::makeSolver(opts(kind, 10));
    const auto fs = solver->fieldCoefficients(g, wl, ys);
    const auto res = solver->solve(g, wl);
    double sumT = 0;
    for (const auto& o : res.transmitted) sumT += o.efficiency;
    const double f0 = flux(fs.samples[0].F);
    for (const auto& s : fs.samples) EXPECT_NEAR(flux(s.F), f0, kind == fg::SolverKind::RCWA ? 1e-12 : 1e-8);
    EXPECT_NEAR(-f0 / (0.5 * std::cos(20 * DEG)), sumT, 1e-10);  // |S_inc| = cos(phi)/2 for |A| = 1
  }
}

TEST(V12_NearField, RKFieldEqualsRCWAFieldInLamellarLayer) {
  const auto g = lamellarStack(false);
  const std::vector<double> ys{0.05, 0.2, 0.35, 0.45};
  for (double delta : {0.0, 90.0}) {
    const fg::Wavelength wl{LAM, 20 * DEG, 0, delta * DEG};
    const auto a = fg::makeSolver(opts(fg::SolverKind::RCWA, 10))->fieldCoefficients(g, wl, ys);
    const auto b = fg::makeSolver(opts(fg::SolverKind::RK, 10))->fieldCoefficients(g, wl, ys);
    for (std::size_t i = 0; i < ys.size(); ++i)
      EXPECT_LT((a.samples[i].F - b.samples[i].F).norm() / a.samples[i].F.norm(), 1e-7) << "y=" << ys[i];
  }
}

TEST(V12_NearField, FFFReconstructionKeepsNormalDContinuous) {
  // TM, lamellar wall at x = 0.45 um (bar eps = 4 | air): D_x = eps E_x must be continuous
  const auto g = lamellarStack(false);
  fg::NearFieldRequest rq;
  rq.horizontal_y_um = 0.2;
  rq.ny = 3;
  rq.nx = 1000;
  double jump[2];
  int k = 0;
  for (auto rc : {fg::Reconstruction::Direct, fg::Reconstruction::FFF}) {
    rq.reconstruction = rc;
    const auto nf = fg::makeSolver(opts(fg::SolverKind::RCWA, 15))->nearField(g, {LAM, 20 * DEG, 0, 90 * DEG}, rq);
    const auto& h = nf.horizontal;
    ASSERT_EQ(h.E[0].size(), 1000u);
    const cplx DL = 4.0 * h.E[0][449], DR = h.E[0][451];
    jump[k++] = std::abs(DL - DR) / std::abs(DR);
  }
  EXPECT_GT(jump[0], 1.0);   // direct synthesis smooths E_x across the wall (wrong)
  EXPECT_LT(jump[1], 0.1);   // FFF-consistent: D_x continuous up to the grid spacing
}

TEST(V12_NearField, CutsAndFilesHaveConsistentShapes) {
  const auto g = lamellarStack(true);
  fg::NearFieldRequest rq;
  rq.nx = 64;
  rq.ny = 40;
  rq.periods = 2;
  rq.y_min_um = -0.2;
  rq.y_max_um = 1.0;
  const auto nf = fg::makeSolver(opts(fg::SolverKind::RCWA, 6))->nearField(g, {LAM, 10 * DEG, 0, 90 * DEG}, rq);
  EXPECT_EQ(nf.vertical.E[0].size(), 64u * 40u);
  EXPECT_EQ(nf.vertical.axis2_um.size(), 40u);
  EXPECT_EQ(nf.horizontal.E[2].size(), 64u);
  EXPECT_TRUE(nf.horizontal.axis2_um.empty());
  const auto dir = std::filesystem::temp_directory_path() / "fg_nf_test";
  fg::NearFieldWriter::write(nf.vertical, dir.string());
  std::ifstream in(dir / "meta.json");
  const auto meta = nlohmann::json::parse(in);
  EXPECT_EQ(meta["cut"], "vertical");
  EXPECT_EQ(meta["axes"]["axis2"]["values_um"].size(), 40u);
  EXPECT_TRUE(std::filesystem::exists(dir / "Ex.npy"));
  EXPECT_EQ(std::filesystem::file_size(dir / "Ex.npy") % 64u, (64u * 40u * 16u) % 64u);
}
