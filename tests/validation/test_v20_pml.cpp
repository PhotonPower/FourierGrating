// Validation M18 (Plan section "M18: Tests und Abnahme", docs/validation.md): PML and isolated
// structures.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <tuple>
#include <utility>

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
  EXPECT_NO_THROW(run(st, fg::Wavelength{0.6328, 10 * DEG, 0, 0}, 20));  // scattered-field formulation
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

double maxDiff(const fg::Result& a, const fg::Result& b) {
  double d = std::abs(a.absorption - b.absorption);
  EXPECT_EQ(a.reflected.size(), b.reflected.size());
  EXPECT_EQ(a.transmitted.size(), b.transmitted.size());
  for (std::size_t k = 0; k < a.reflected.size() && k < b.reflected.size(); ++k)
    d = std::max(d, std::abs(a.reflected[k].efficiency - b.reflected[k].efficiency));
  for (std::size_t k = 0; k < a.transmitted.size() && k < b.transmitted.size(); ++k)
    d = std::max(d, std::abs(a.transmitted[k].efficiency - b.transmitted[k].efficiency));
  return d;
}

// ---------------------------------------------------------------- V20d
// Scattered-field formulation without PML: background + scattered field (alpha_n = alpha0 + n K) equals
// the total-field computation (Plan eq:scat_source, eq:scat_ode; sigma = (M - M_bg) F_bg is exactly the
// difference of the two discrete systems). Measured: RCWA 1e-15, RK 3e-12.
TEST(M18_Scattered, V20dPeriodicEqualsTotalField) {
  auto lam = base(1.0);
  lam.layers = {entry(fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"film", "superstrate"}})};
  const fg::Wavelength te{0.6328, 20 * DEG, 0, 0}, tm{0.6328, 20 * DEG, 0, 90 * DEG},
      con{0.6328, 20 * DEG, 15 * DEG, 30 * DEG};
  for (const auto& wl : {te, tm, con})
    for (auto kind : {fg::SolverKind::RCWA, fg::SolverKind::RK}) {
      fg::SolverOptions o;
      o.kind = kind;
      o.Nx = 10;
      o.rkTol = 1e-10;
      const auto a = fg::makeSolver(o)->solve(lam, wl);
      o.scatteredField = true;
      const auto b = fg::makeSolver(o)->solve(lam, wl);
      EXPECT_LT(maxDiff(a, b), kind == fg::SolverKind::RCWA ? 1e-12 : 1e-10) << static_cast<int>(kind);
    }
  // metallic profile (RCWA staircase, background split at the profile height)
  auto prof = base(0.8);
  fg::Profile1DLayer p;
  p.depth_um = 0.2;
  p.x0_um = 0.1;
  p.below = "metal";
  p.above = "superstrate";
  prof.materials["metal"] = fg::Material::fromIndex("metal", cplx(0.6, 2.5));
  prof.layers = {entry(p)};
  fg::SolverOptions o;
  o.Nx = 8;
  const auto a = fg::makeSolver(o)->solve(prof, con);
  o.scatteredField = true;
  EXPECT_LT(maxDiff(a, fg::makeSolver(o)->solve(prof, con)), 1e-10);
  // Contour2D pillars (2D), conical
  fg::Stack pil = base(0.8);
  pil.dz_um = 0.7;
  fg::Contour2DLayer c;
  c.thickness_um = 0.3;
  c.background = "superstrate";
  fg::geometry::Region r;
  r.outer = fg::geometry::Contour2D::circle({0.4, 0.35}, 0.2);
  r.material = "film";
  r.normalize(1e-9);
  c.regions = {r};
  pil.layers = {entry(c)};
  o = {};
  o.Nx = 3;
  o.Nz = 3;
  const auto pa = fg::makeSolver(o)->solve(pil, con);
  o.scatteredField = true;
  EXPECT_LT(maxDiff(pa, fg::makeSolver(o)->solve(pil, con)), 1e-10);
}

// Near field of the scattered-field formulation (background added back) = total-field near field.
TEST(M18_Scattered, NearFieldEqualsTotalField) {
  auto lam = base(1.0);
  lam.layers = {entry(fg::Lamellar1DLayer{0.4, {0.0, 0.45, 1.0}, {"film", "superstrate"}})};
  const fg::Wavelength con{0.6328, 20 * DEG, 15 * DEG, 30 * DEG};
  fg::NearFieldRequest req;
  req.nx = 9;
  req.ny = 7;
  req.y_min_um = -0.2;
  req.y_max_um = 0.6;
  req.horizontal_y_um = 0.2;
  fg::SolverOptions o;
  o.Nx = 10;
  const auto a = fg::makeSolver(o)->nearField(lam, con, req);
  o.scatteredField = true;
  const auto b = fg::makeSolver(o)->nearField(lam, con, req);
  double d = 0;
  for (int c = 0; c < 3; ++c)
    for (std::size_t k = 0; k < a.vertical.E[c].size(); ++k)
      d = std::max({d, std::abs(a.vertical.E[c][k] - b.vertical.E[c][k]), std::abs(a.vertical.H[c][k] - b.vertical.H[c][k])});
  EXPECT_LT(d, 1e-10);
  // scattered part: total minus the background of the incidence order
  req.scattered = true;
  const auto sc = fg::makeSolver(o)->nearField(lam, con, req);
  ASSERT_EQ(sc.vertical.E[0].size(), a.vertical.E[0].size());
  double ds = 0;
  for (std::size_t k = 0; k < a.vertical.E[0].size(); ++k) ds = std::max(ds, std::abs(a.vertical.E[0][k] - sc.vertical.E[0][k]));
  EXPECT_GT(ds, 1e-3);  // the background is not zero
}

// ---------------------------------------------------------------- V20e
// Normal incidence with PML: the scattered-field formulation equals total field minus background
// (measured: identical to 1e-10 in all outputs).
TEST(M18_Scattered, V20eNormalIncidence) {
  for (const auto& wl : {TE, TM}) {
    auto st = ridge(3.0, 0.5, 0.3, 0.3);
    fg::SolverOptions o;
    o.Nx = 40;
    const auto a = fg::makeSolver(o)->solve(st, wl);
    o.scatteredField = true;
    const auto b = fg::makeSolver(o)->solve(st, wl);
    EXPECT_NEAR(a.scattering.up, b.scattering.up, 1e-8);
    EXPECT_NEAR(a.scattering.down, b.scattering.down, 1e-8);
    EXPECT_NEAR(a.scattering.lateral, b.scattering.lateral, 1e-8);
    EXPECT_NEAR(a.scattering.absorbed, b.scattering.absorbed, 1e-8);
  }
}

// ---------------------------------------------------------------- V20f
/// Trench of width a in a film (n = 2.1, thickness h) on glass, centred in a cell of width L.
fg::Stack trench(double L, double a, double h, double w, double beta = 1.0) {
  auto st = pml(base(L), w, beta);
  st.layers = {entry(fg::Lamellar1DLayer{h, {0.0, 0.5 * (L - a), 0.5 * (L + a), L}, {"film", "superstrate", "film"}})};
  return st;
}

// Ridge and trench under oblique incidence (scattered-field formulation), TE/TM, phi = 30, 60 deg.
// The total scattered power (up + down + lateral) converges in the cell width (measured 1e-4 over
// 2.4 / 3.2 / 4.0 / 4.8 um windows); the split into up and down converges only like 1/L (grazing power
// crosses the planes outside the window). w, beta change all outputs at a fixed window by < 5e-4;
// |absorbed| < 1e-3 (lossless; trench TM 7.6e-4).
TEST(M18_Scattered, V20fRidgeAndTrenchOblique) {
  auto N = [](double L) { return static_cast<int>(12 * L / 0.6328); };
  for (int kind = 0; kind < 2; ++kind)
    for (double delta : {0.0, 90.0})
      for (double phi : {30.0, 60.0}) {  // phi = 0: V20b, V20e
        const fg::Wavelength wl{0.6328, phi * DEG, 0, delta * DEG};
        auto make = [&](double win, double w, double beta) {
          const double L = win + 2 * w;
          return run(kind ? trench(L, 0.5, 0.2, w, beta) : ridge(L, 0.5, 0.3, w, beta), wl, N(L));
        };
        const std::string what =
            std::string(kind ? "trench " : "ridge ") + (delta == 0 ? "TE " : "TM ") + std::to_string(phi);
        const auto a = make(3.2, 0.4, 1.0), b = make(4.0, 0.5, 1.0);
        EXPECT_LT(std::abs(a.scattering.total() - b.scattering.total()), 1e-3) << what;
        const auto c = make(2.4, 0.3, 1.0), d = make(2.4, 0.6, 4.0);
        EXPECT_LT(std::abs(c.scattering.up - d.scattering.up), 1e-3) << what;
        EXPECT_LT(std::abs(c.scattering.down - d.scattering.down), 1e-3) << what;
        EXPECT_LT(std::abs(c.scattering.total() - d.scattering.total()), 1e-3) << what;
        for (const auto* r : {&a, &b, &c, &d}) EXPECT_LT(std::abs(r->scattering.absorbed), 1e-3) << what;
        EXPECT_GT(a.scattering.total(), 0.1) << what;
      }
  // RK with sources = RCWA (y-invariant lamellae)
  const fg::Wavelength wl{0.6328, 30 * DEG, 0, 90 * DEG};
  const auto st = ridge(3.0, 0.5, 0.3, 0.3);
  const auto rc = run(st, wl, 40), rk = run(st, wl, 40, 0, fg::SolverKind::RK);
  EXPECT_NEAR(rc.scattering.up, rk.scattering.up, 1e-7);
  EXPECT_NEAR(rc.scattering.down, rk.scattering.down, 1e-7);
}

// V20f reciprocity: far-field amplitude f(alpha_s; alpha_i) = f(-alpha_i; -alpha_s) for an asymmetric
// structure (film and ridge side by side). f(alpha) = beta A(alpha), A from the scattered field on a line
// above the structure (angular spectrum over the physical window with a cos^2 taper; incident phase at the
// cell centre and the top of the stack). The truncated line limits the accuracy: measured complex deviation
// 2.5-6.7 % at L = 8 um (TE/TM, two angle pairs), falling with L (4 um: 15-40 %).
TEST(M18_Scattered, V20fReciprocity) {
  const double L = 8.0, w = 0.5, y0 = 0.4, top = 0.3, k = 2 * fg::PI / 0.6328, c = 0.5 * L;
  auto st = pml(base(L), w);
  st.layers = {entry(fg::Lamellar1DLayer{top, {0.0, c - 0.3, c, c + 0.2, L}, {"superstrate", "film", "ridge", "superstrate"}})};
  fg::SolverOptions o;
  o.Nx = static_cast<int>(12 * L / 0.6328);
  fg::NearFieldRequest req;
  req.nx = static_cast<int>(80 * L);
  req.ny = 3;
  req.y_min_um = 0.0;
  req.y_max_um = 0.5;
  req.horizontal_y_um = y0;
  req.reconstruction = fg::Reconstruction::Direct;
  req.scattered = true;
  const auto solver = fg::makeSolver(o);
  // f(alpha) of the scattered field for incidence angle phi (deg)
  auto farField = [&](double phi, double delta, double alpha) {
    const auto nf = solver->nearField(st, fg::Wavelength{0.6328, phi * DEG, 0, delta * DEG}, req);
    const auto& F = delta == 0 ? nf.horizontal.E[2] : nf.horizontal.H[2];
    const auto& x = nf.horizontal.axis1_um;
    const double a0 = k * std::sin(phi * DEG), b0 = k * std::cos(phi * DEG), dx = x[1] - x[0], hw = 0.5 * L - w;
    const double b = std::sqrt(k * k - alpha * alpha);
    cplx A = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double xs = x[i] - c;
      if (std::abs(xs) >= hw) continue;
      const double taper = std::pow(std::cos(0.5 * fg::PI * xs / hw), 2);
      A += F[i] * taper * std::exp(cplx(0, -alpha * xs));
    }
    // incident amplitude referred to (x = c, y = top)
    return b * A * dx / (2 * fg::PI) * std::exp(cplx(0, -b * y0)) * std::exp(cplx(0, -a0 * c - b0 * top));
  };
  for (double delta : {0.0, 90.0}) {
    const double ti = 20, ts = -35;
    const cplx f1 = farField(ti, delta, k * std::sin(ts * DEG)), f2 = farField(-ts, delta, -k * std::sin(ti * DEG));
    std::printf("V20f reciprocity %s: |f1| = %.5f, |f2| = %.5f, |f1 - f2| / |f1| = %.3e\n", delta == 0 ? "TE" : "TM",
                std::abs(f1), std::abs(f2), std::abs(f1 - f2) / std::abs(f1));
    EXPECT_LT(std::abs(f1 - f2), 0.1 * std::abs(f1)) << delta;
  }
}

// Near field of the scattered field with PML under oblique incidence: finite, and decaying into the
// zones (the scattered field is outgoing).
TEST(M18_Scattered, NearFieldOblique) {
  const auto st = ridge(3.0, 0.5, 0.3, 0.3);
  fg::NearFieldRequest req;
  req.nx = 60;
  req.ny = 5;
  req.y_min_um = 0.35;
  req.y_max_um = 0.6;
  req.horizontal_y_um = 0.4;
  req.scattered = true;
  fg::SolverOptions o;
  o.Nx = 40;
  const auto nf = fg::makeSolver(o)->nearField(st, fg::Wavelength{0.6328, 30 * DEG, 0, 0}, req);
  const auto& Ez = nf.horizontal.E[2];
  ASSERT_EQ(Ez.size(), 60u);
  double inner = 0, edge = 0;
  for (std::size_t k = 0; k < Ez.size(); ++k) {
    ASSERT_TRUE(std::isfinite(std::abs(Ez[k])));
    const double x = 3.0 * static_cast<double>(k) / 60;
    if (x > 1.0 && x < 2.0) inner = std::max(inner, std::abs(Ez[k]));
    if (x < 0.05 || x > 2.95) edge = std::max(edge, std::abs(Ez[k]));
  }
  EXPECT_GT(inner, 1e-3);
  EXPECT_LT(edge, 0.2 * inner);
}

/// Slab waveguide along y: core n = 2.0 of width a centred in a cell of width L (air cladding).
fg::LayerEntry slab(double L, double a) {
  return entry(fg::Lamellar1DLayer{1.0, {0.0, 0.5 * (L - a), 0.5 * (L + a), L}, {"superstrate", "core", "superstrate"}});
}
fg::Stack waveguideEnd(double L, double a, double w, double beta = 1.0) {
  auto st = pml(base(L), w, beta);
  st.materials["core"] = fg::Material::fromIndex("core", cplx(2.0, 0.0));
  st.substrate = "superstrate";  // air below the facet
  st.superstrateStructure = slab(L, a);
  st.aperiodic.modeExcitation = true;
  return st;
}

// ---------------------------------------------------------------- V20c (modal excitation)
// A waveguide continued by the same waveguide: no reflection, full transmission of the incident guided
// mode (TE0 and TM0; measured 3e-6, the PML perturbs the guided modes slightly).
TEST(M18_Modal, ThroughWaveguide) {
  for (int mode : {0, 1}) {
    auto st = waveguideEnd(3.0, 0.4, 0.4);
    st.aperiodic.mode = mode;
    st.substrateStructure = slab(3.0, 0.4);
    st.layers = {slab(3.0, 0.4)};
    std::get<fg::Lamellar1DLayer>(st.layers[0].spec).thickness_um = 0.5;
    const auto r = run(st, TE, 40);
    ASSERT_TRUE(r.modal.present);
    EXPECT_NEAR(r.modal.incident.teFraction, mode == 0 ? 1.0 : 0.0, 1e-6);
    for (const auto& m : r.modal.reflected) EXPECT_LT(m.power, 1e-10);
    ASSERT_GT(r.modal.transmitted.size(), static_cast<std::size_t>(mode));
    EXPECT_NEAR(r.modal.transmitted[static_cast<std::size_t>(mode)].power, 1.0, 1e-5);
    EXPECT_NEAR(std::abs(r.modal.transmitted[static_cast<std::size_t>(mode)].nEff - r.modal.incident.nEff), 0.0, 1e-10);
    EXPECT_LT(std::abs(r.modal.radiated), 1e-5);
  }
}

// Abrupt end of a slab waveguide (core n = 2.0, 0.4 um, air), TE0 incident (Plan V20c). Measured: R(TE0)
// = 0.15999, R(TE2, n_eff = 1.067) = 0.0191, radiated 0.8210; converged in N (30 / 50 / 80), cell width
// and w, beta to 1e-4. The comparison with Silberstein et al. (JOSA A 18, 2865, 2001) is open: the value
// could not be checked bibliographically (docs/validation.md).
TEST(M18_Modal, V20cAbruptEnd) {
  const auto a = run(waveguideEnd(3.0, 0.4, 0.4), TE, 50);
  ASSERT_TRUE(a.modal.present);
  ASSERT_GE(a.modal.reflected.size(), 5u);
  const double R0 = a.modal.reflected[0].power;
  EXPECT_NEAR(R0, 0.15999, 1e-4);
  EXPECT_NEAR(a.modal.incident.nEff.real(), 1.9047, 1e-3);
  for (const auto& [L, w, beta, N] : {std::tuple{3.0, 0.4, 1.0, 80}, {4.0, 0.5, 1.0, 70}, {3.0, 0.6, 3.0, 50}}) {
    const auto b = run(waveguideEnd(L, 0.4, w, beta), TE, N);
    EXPECT_NEAR(b.modal.reflected[0].power, R0, 1e-4) << L << " " << w << " " << beta << " " << N;
    EXPECT_NEAR(b.modal.radiated, a.modal.radiated, 1e-3) << L << " " << w << " " << beta << " " << N;
  }
  // symmetry: no TM and no odd TE mode
  for (const auto& m : a.modal.reflected)
    if (m.teFraction < 0.5 || m.index == 2) EXPECT_LT(m.power, 1e-8) << m.index;
  EXPECT_TRUE(a.modal.transmitted.empty());  // homogeneous substrate: no guided modes
  EXPECT_GT(a.modal.radiated, 0.7);
}

TEST(M18_Modal, Configuration) {
  const std::string j = R"({"dimension": "1D", "period_um": {"x": 3.0}, "truncation": {"Nx": 40},
    "materials": {"air": {"n": 1.0}, "core": {"n": 2.0}},
    "superstrate": "air", "substrate": "air",
    "incidence": {"lambda_um": 0.6328},
    "superstrate_structure": {"type": "lamellar1d", "breakpoints_um": [0.0, 1.3, 1.7, 3.0], "materials": ["air", "core", "air"]},
    "aperiodic": {"pml": {"x": {"width_um": 0.4}}, "excitation": {"mode": 0}}})";
  const auto cfg = fg::io::parseConfig(j);
  EXPECT_TRUE(cfg.stack.aperiodic.modeExcitation);
  ASSERT_TRUE(cfg.stack.superstrateStructure.has_value());
  EXPECT_FALSE(cfg.stack.substrateStructure.has_value());
  fg::SolverOptions o;
  o.Nx = 40;
  const auto r = fg::makeSolver(o)->solve(cfg.stack, cfg.scan.front());
  EXPECT_NEAR(r.modal.reflected[0].power, 0.15999, 2e-4);
  // the core reaches into the zone
  auto bad = j;
  bad.replace(bad.find("1.3, 1.7"), 8, "0.3, 1.7");
  EXPECT_THROW(fg::io::parseConfig(bad), fg::Error);
}

// ---------------------------------------------------------------- 2D (V20g)
/// Contact hole: round hole (radius r) in a film (n = 2.1, thickness h) on glass, cell L x L, PML in x and z.
fg::Stack contactHole(double L, double r, double h, double w, double beta = 1.0) {
  auto st = pml(base(L), w, beta);
  st.dz_um = L;
  st.aperiodic.z.width_um = w;
  st.aperiodic.z.beta = beta;
  fg::Contour2DLayer c;
  c.thickness_um = h;
  c.background = "film";
  fg::geometry::Region reg;
  reg.outer = fg::geometry::Contour2D::circle({0.5 * L, 0.5 * L}, r);
  reg.material = "superstrate";
  reg.normalize(1e-9);
  c.regions = {reg};
  st.layers = {entry(c)};
  return st;
}

// A z-invariant trench in a 2D cell (PML in x, periodic in z) reproduces the 1D trench: 2D outputs are
// cross sections over the z period, i.e. dz times the 1D widths (normal incidence: total field; 30 deg:
// scattered-field formulation in 2D).
TEST(M18_Pml2D, ZInvariantTrenchEquals1D) {
  for (double phi : {0.0, 30.0}) {
    const fg::Wavelength wl{0.6328, phi * DEG, 0, 0};
    auto st1 = trench(3.0, 0.5, 0.2, 0.3);
    st1.layers[0].factorization = fg::Factorization::Laurent;
    const auto r1 = run(st1, wl, 40);
    auto st = pml(base(3.0), 0.3);
    st.dz_um = 0.6;
    fg::Polygon2DLayer p;
    p.thickness_um = 0.2;
    p.background = "film";
    p.inclusions.push_back({{{1.25, 0.0}, {1.75, 0.0}, {1.75, 0.6}, {1.25, 0.6}}, "superstrate"});
    st.layers = {entry(p, fg::Factorization::Laurent)};
    const auto r2 = run(st, wl, 40, 1);
    EXPECT_NEAR(r2.scattering.up, 0.6 * r1.scattering.up, 1e-10) << phi;
    EXPECT_NEAR(r2.scattering.down, 0.6 * r1.scattering.down, 1e-10) << phi;
    EXPECT_NEAR(r2.scattering.lateral, 0.6 * r1.scattering.lateral, 1e-10) << phi;
    EXPECT_NEAR(r2.scattering.absorbed, 0.6 * r1.scattering.absorbed, 1e-10) << phi;
  }
}

// V20g: contact hole (round hole r = 0.25 um in a 0.3 um film, n = 2.1, on glass), PML in x and z. By the
// symmetry of hole and cell, x- and z-polarized normal incidence scatter identically (checks the z-PML
// against the x-PML). Convergence in N (measured, cell 2 um): total 0.501 / 0.568 / 0.570 at N = 4 / 6 / 8,
// |absorbed| 3e-2 / 1e-2 / 1e-3; the study at larger N belongs on the target system (docs/validation.md).
// Contour2D = CAD body: tests/cad/test_cad.cpp (CAD_V20g).
TEST(M18_Pml2D, V20gContactHoleSymmetry) {
  const auto st = contactHole(2.0, 0.25, 0.3, 0.3);
  const auto a = run(st, TE, 5, 5), b = run(st, TM, 5, 5);
  EXPECT_NEAR(a.scattering.up, b.scattering.up, 1e-10);
  EXPECT_NEAR(a.scattering.down, b.scattering.down, 1e-10);
  EXPECT_NEAR(a.scattering.lateral, b.scattering.lateral, 1e-10);
  EXPECT_GT(a.scattering.total(), 0.3);
  // oblique incidence (scattered field, Bloch-free in x) runs and scatters a similar amount
  const auto c = run(st, fg::Wavelength{0.6328, 30 * DEG, 0, 0}, 5, 5);
  EXPECT_GT(c.scattering.total(), 0.3);
  EXPECT_LT(std::abs(c.scattering.absorbed), 0.05);
}
