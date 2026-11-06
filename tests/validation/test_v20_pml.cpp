// Validation M18 (Plan section "M18: Tests und Abnahme", docs/validation.md): PML and isolated
// structures.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;

namespace {

constexpr double DEG = fg::PI / 180.0;

fg::LayerEntry entry(fg::LayerSpec s, fg::Factorization f = fg::Factorization::Auto) {
  fg::LayerEntry e;
  e.spec = std::move(s);
  e.factorization = f;
  return e;
}

fg::Result run(const fg::Stack& st, const fg::Wavelength& wl, int Nx, int Nz = 0,
               fg::SolverKind k = fg::SolverKind::RCWA) {
  fg::SolverOptions o;
  o.kind = k;
  o.Nx = Nx;
  o.Nz = Nz;
  o.rkTol = 1e-10;
  return fg::makeSolver(o)->solve(st, wl);
}

fg::Stack base(double dx) {
  fg::Stack st;
  st.dx_um = dx;
  st.materials["superstrate"] = fg::Material::fromIndex("air", 1.0);
  st.materials["substrate"] = fg::Material::fromIndex("glass", 1.5);
  st.materials["film"] = fg::Material::fromIndex("film", cplx(2.1, 0.0));
  st.materials["ridge"] = fg::Material::fromIndex("ridge", cplx(1.5, 0.0));
  return st;
}

fg::Stack pml(fg::Stack st, double w, double beta = 1.0) {
  st.aperiodic.enabled = true;
  st.aperiodic.x.width_um = w;
  st.aperiodic.x.beta = beta;
  return st;
}

/// Isolated ridge of width a, height h, centred in a cell of width L.
fg::Stack ridge(double L, double a, double h, double w, double beta = 1.0) {
  auto st = pml(base(L), w, beta);
  st.layers = {entry(fg::Lamellar1DLayer{h, {0.0, 0.5 * (L - a), 0.5 * (L + a), L}, {"superstrate", "ridge", "superstrate"}})};
  return st;
}

const fg::Wavelength TE{0.6328, 0, 0, 0}, TM{0.6328, 0, 0, 90 * DEG};

}  // namespace

// ---------------------------------------------------------------- V20a
// Homogeneous stack with PML at normal incidence: the constant field is exact in stretched coordinates
// (c_nm = delta_n0 delta_m0), so the field equals Fresnel/Airy (periodic N = 0) and nothing scatters.
TEST(M18_Pml, V20aHomogeneousStackIsFresnel) {
  auto per = base(2.0);
  per.layers = {entry(fg::HomogeneousLayer{0.37, "film"})};
  const auto st = pml(per, 0.3);
  const std::vector<double> ys{-0.2, 0.1, 0.3, 0.5};
  for (const auto& wl : {TE, TM}) {
    fg::SolverOptions o;
    o.Nx = 0;
    const auto ref = fg::makeSolver(o)->fieldCoefficients(per, wl, ys);
    o.Nx = 20;
    const auto f = fg::makeSolver(o)->fieldCoefficients(st, wl, ys);
    const int nf = 41, c = 20;
    for (std::size_t k = 0; k < ys.size(); ++k)
      for (int comp = 0; comp < 4; ++comp)
        for (int r = 0; r < nf; ++r) {
          const cplx v = f.samples[k].F(comp * nf + r), e = r == c ? ref.samples[k].F(comp) : cplx(0.0);
          EXPECT_LT(std::abs(v - e), 1e-10) << "y = " << ys[k] << " component " << comp << " order " << r - c;
        }
    const auto r = run(st, wl, 20);
    ASSERT_TRUE(r.scattering.present);
    EXPECT_LT(std::abs(r.scattering.up) + std::abs(r.scattering.down) + std::abs(r.scattering.absorbed), 1e-12);
    EXPECT_TRUE(r.reflected.empty() && r.transmitted.empty());
  }
}

// ---------------------------------------------------------------- V20b
// Isolated dielectric ridge (n = 1.5, 0.5 um x 0.3 um) on glass, normal incidence, TE and TM. Measured
// (docs/validation.md): cell 4 / 6 / 8 x ridge width (w = 0.15 L) -> up 0.00909 / 0.01030 / 0.01067 (TE), down
// 0.9072 / 0.9101 / 0.9104; w, beta change the widths by < 2e-5; lossless: |absorbed| <= 2e-4 (the
// truncated PML system conserves energy only up to the discretization error).
TEST(M18_Pml, V20bIsolatedRidge) {
  auto N = [](double L) { return static_cast<int>(12 * L / 0.6328); };
  for (const auto& wl : {TE, TM}) {
    const char* pol = wl.delta == 0 ? "TE" : "TM";
    const auto r6 = run(ridge(3.0, 0.5, 0.3, 0.3), wl, N(3.0));  // physical region 2.4 um
    const auto r8 = run(ridge(4.0, 0.5, 0.3, 0.4), wl, N(4.0));
    EXPECT_LT(std::abs(r8.scattering.up - r6.scattering.up), 1e-3) << pol;
    EXPECT_LT(std::abs(r8.scattering.down - r6.scattering.down), 1e-3) << pol;
    EXPECT_GT(r6.scattering.up, 0.004) << pol;
    EXPECT_GT(r6.scattering.down, 0.7) << pol;
    // the same physical region (2.4 um) with other zone widths and strengths
    for (const auto& [w, beta] : {std::pair{0.3, 4.0}, {0.45, 1.0}, {0.6, 0.5}}) {
      const double L = 2.4 + 2 * w;
      const auto r = run(ridge(L, 0.5, 0.3, w, beta), wl, N(L));
      EXPECT_LT(std::abs(r.scattering.up - r6.scattering.up), 1e-4) << pol << " w = " << w << " beta = " << beta;
      EXPECT_LT(std::abs(r.scattering.down - r6.scattering.down), 1e-4) << pol << " w = " << w << " beta = " << beta;
      EXPECT_LT(std::abs(r.scattering.absorbed), 3e-4) << pol << " w = " << w << " beta = " << beta;
    }
    EXPECT_LT(std::abs(r6.scattering.absorbed), 3e-4) << pol;
    EXPECT_LT(std::abs(r8.scattering.absorbed), 3e-4) << pol;
  }
}

// An absorbing ridge: the box flux gives the absorbed power (positive, scattering unchanged in sign).
TEST(M18_Pml, AbsorbingRidge) {
  auto st = ridge(3.0, 0.5, 0.3, 0.45);
  st.materials["ridge"] = fg::Material::fromIndex("ridge", cplx(1.5, 0.1));
  const auto r = run(st, TE, 56);
  EXPECT_GT(r.scattering.absorbed, 0.05);
  EXPECT_GT(r.scattering.up, 0.0);
}

// ---------------------------------------------------------------- configuration and preconditions
TEST(M18_Pml, Configuration) {
  const std::string j = R"({"dimension": "1D", "period_um": {"x": 3.0}, "truncation": {"Nx": 40},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "bar": {"n": 1.5}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"lambda_um": 0.6328, "phi_deg": 0},
    "layers": [{"type": "lamellar1d", "thickness_um": 0.3, "breakpoints_um": [0.0, 1.25, 1.75, 3.0], "materials": ["air", "bar", "air"]}],
    "aperiodic": {"pml": {"x": {"width_um": 0.45, "beta": 2.0}}, "excitation": "plane_wave"}})";
  const auto cfg = fg::io::parseConfig(j);
  EXPECT_TRUE(cfg.stack.aperiodic.enabled);
  EXPECT_DOUBLE_EQ(cfg.stack.aperiodic.x.width_um, 0.45);
  EXPECT_DOUBLE_EQ(cfg.stack.aperiodic.x.beta, 2.0);
  EXPECT_FALSE(cfg.stack.aperiodic.z.active());
  // the ridge reaches into the zone
  auto bad = j;
  bad.replace(bad.find("1.25, 1.75"), 10, "0.25, 1.75");
  EXPECT_THROW(fg::io::parseConfig(bad), fg::Error);
  // z PML in 1D, unknown excitation
  bad = j;
  bad.replace(bad.find("\"beta\": 2.0}}"), 13, "\"beta\": 2.0}, \"z\": {\"width_um\": 0.1}}");
  EXPECT_THROW(fg::io::parseConfig(bad), fg::Error);
  bad = j;
  bad.replace(bad.find("\"plane_wave\""), 12, "\"dipole\"");
  EXPECT_THROW(fg::io::parseConfig(bad), fg::Error);
  // oblique incidence in the total-field formulation, PML combined with ASR
  auto st = ridge(3.0, 0.5, 0.3, 0.45);
  EXPECT_THROW(run(st, fg::Wavelength{0.6328, 10 * DEG, 0, 0}, 20), fg::Error);
  st.transform.x.type = fg::AxisTransformSpec::Type::ASR;
  st.transform.x.eta = 0.5;
  EXPECT_THROW(run(st, TE, 20), fg::Error);
}

// The background stack: lamellae -> the material at x = 0, a flat-in-the-zone profile -> two layers.
TEST(M18_Pml, BackgroundStack) {
  auto st = ridge(3.0, 0.5, 0.3, 0.45);
  fg::Profile1DLayer p;
  p.shape = fg::Profile1DLayer::Shape::PiecewiseLinear;
  p.depth_um = 0.2;
  p.xs_um = {0.0, 1.2, 1.5, 1.8};
  p.ys_um = {0.05, 0.05, 0.2, 0.05};
  p.below = "ridge";
  p.above = "superstrate";
  st.layers.push_back(entry(p));
  fg::checkPmlZones(st);
  const auto bg = fg::backgroundStack(st);
  ASSERT_EQ(bg.layers.size(), 3u);
  const auto& l0 = std::get<fg::HomogeneousLayer>(bg.layers[0].spec);
  EXPECT_EQ(l0.material, "superstrate");
  EXPECT_DOUBLE_EQ(l0.thickness_um, 0.3);
  const auto& l1 = std::get<fg::HomogeneousLayer>(bg.layers[1].spec);
  const auto& l2 = std::get<fg::HomogeneousLayer>(bg.layers[2].spec);
  EXPECT_EQ(l1.material, "ridge");
  EXPECT_NEAR(l1.thickness_um, 0.05, 1e-15);
  EXPECT_EQ(l2.material, "superstrate");
  EXPECT_NEAR(l2.thickness_um, 0.15, 1e-15);
  // not flat in the zone
  std::get<fg::Profile1DLayer>(st.layers[1].spec).ys_um[0] = 0.06;
  EXPECT_THROW(fg::checkPmlZones(st), fg::Error);
}
