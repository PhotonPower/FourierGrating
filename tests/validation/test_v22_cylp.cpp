// V22 (M21, Plan section "Tests und Abnahme (M21)"): P-fold rotational symmetry with coupled harmonics.
// V22b: guided modes of an eccentric fiber core (P = 1) equal those of the centred core (Plan "Befund Prototyp").
// V22a: rotationally symmetric structures through the class path equal M20; V22d: selection rules; V22e: mirror
// symmetry (classes m0 and -m0).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "fg/cyl/coupled.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"
#include "fg/io/config.hpp"
#include "fg/solvers/solver.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;
namespace cyl = fg::cyl;

namespace {

/// Propagation constants beta > sqrt(eps_clad) of the guided modes, descending (M^2 = diag(PR, RP), eigenvalues
/// -beta^2 of PR, Plan eq:cyl_M; as in V21b).
std::vector<double> guidedBetas(const VecC& k, const cyl::CylQ& Q, double epsClad) {
  const MatC M = cyl::propagationMatrix(k, Q);
  const Eigen::Index n2 = M.rows() / 2;
  const MatC PR = M.topRightCorner(n2, n2) * M.bottomLeftCorner(n2, n2);
  const VecC mu = PR.eigenvalues();
  std::vector<double> b;
  for (Eigen::Index i = 0; i < mu.size(); ++i) {
    const cplx beta = std::sqrt(-mu(i));
    if (std::abs(beta.imag()) < 1e-6 * std::abs(beta) && beta.real() > std::sqrt(epsClad) + 1e-6) b.push_back(beta.real());
  }
  std::sort(b.rbegin(), b.rend());
  return b;
}

fg::Result solveConfig(const std::string& json) {
  fg::io::RunConfig cfg = fg::io::parseConfig(json);
  return fg::makeSolver(cfg.solver)->solve(cfg.stack, cfg.scan.front());
}

/// Pillar of n = 1.45 (height 0.3 um) on glass in air; `symmetry` and `layer` as JSON, the rest fixed.
std::string pillar(const std::string& cylindrical, const std::string& incidence, const std::string& layer) {
  return R"({"cylindrical": {)" + cylindrical + R"(, "spectrum": {"kappa_max": 3}},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "p": {"n": 1.45}},
    "superstrate": "air", "substrate": "glass", "incidence": )" + incidence + R"(, "layers": [)" + layer + "]}";
}

const char* kBeam = R"({"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.25, "polarization": "x"})";
const char* kOblique = R"({"type": "plane_wave", "lambda_um": 1.0, "phi_deg": 30, "psi_deg": 0, "delta_deg": 30})";
const char* kRings = R"({"type": "rings", "thickness_um": 0.3, "radii_um": [0.3], "materials": ["p", "air"]})";
const char* kCircle = R"({"type": "contour2d", "thickness_um": 0.3, "background": "air",
    "regions": [{"material": "p", "circle": {"center_um": [0, 0], "radius_um": 0.3}}]})";
const char* kSquare = R"({"type": "polygon2d", "thickness_um": 0.3, "background": "air",
    "inclusions": [{"material": "p", "rectangle_um": {"center": [0, 0], "size": [0.5, 0.5]}}]})";

/// Value of harmonic m in a result (per-harmonic list), NaN if absent.
double perHarmonic(const fg::CylindricalResult& r, const std::vector<double>& v, int m) {
  for (std::size_t i = 0; i < r.harmonics.size(); ++i)
    if (r.harmonics[i] == m) return v[i];
  return std::nan("");
}

}  // namespace

// Rings layer (Li in r) through the class path: the class is block-diagonal with the M20 blocks, so R and T of every
// harmonic equal M20 (1e-12, Plan); the harmonics without incidence stay dark.
TEST(V22a, RingsThroughTheClassPathEqualM20) {
  const auto ref = solveConfig(pillar(R"("factorization": "li_r")", kBeam, kRings)).cylindrical;
  for (const int P : {1, 2, 3}) {
    const auto c = solveConfig(pillar(R"("factorization": "li_r", "symmetry": {"rotation": )" + std::to_string(P) +
                                          R"(, "L": 1})",
                                      kBeam, kRings))
                       .cylindrical;
    ASSERT_TRUE(c.present);
    EXPECT_EQ(c.rotation, P);
    std::printf("V22a rings, P = %d: R = %.15f T = %.15f (M20 %.15f %.15f), %zu harmonics\n", P, c.reflected,
                c.transmitted, ref.reflected, ref.transmitted, c.harmonics.size());
    EXPECT_NEAR(c.reflected, ref.reflected, 1e-12) << P;
    EXPECT_NEAR(c.transmitted, ref.transmitted, 1e-12) << P;
    for (const int m : {-1, 1}) {
      EXPECT_NEAR(perHarmonic(c, c.reflectedPerHarmonic, m), perHarmonic(ref, ref.reflectedPerHarmonic, m), 1e-12);
      EXPECT_NEAR(perHarmonic(c, c.transmittedPerHarmonic, m), perHarmonic(ref, ref.transmittedPerHarmonic, m), 1e-12);
    }
    for (std::size_t i = 0; i < c.harmonics.size(); ++i)
      if (std::abs(c.harmonics[i]) != 1) EXPECT_LT(c.reflectedPerHarmonic[i] + c.transmittedPerHarmonic[i], 1e-14);
  }
}

// Centred circle as contour2d (angular coefficients, Laurent, coupledQ and coupledSource by quadrature) against the
// M20 rings layer (Laurent, Lommel), oblique plane wave: cross sections per harmonic and in total.
TEST(V22a, CircleThroughTheClassPathEqualsM20) {
  const auto ref = solveConfig(pillar(R"("factorization": "laurent")", kOblique, kRings)).cylindrical;
  const auto c =
      solveConfig(pillar(R"("factorization": "laurent", "symmetry": {"rotation": 3, "L": 2})", kOblique, kCircle)).cylindrical;
  ASSERT_TRUE(c.present && c.planeWave);
  std::printf("V22a circle, P = 3, L = %d: sigma up %.15f down %.15f (M20 %.15f %.15f)\n", c.classL, c.scatteredUp_um2,
              c.scatteredDown_um2, ref.scatteredUp_um2, ref.scatteredDown_um2);
  // P = 3, L = 2 holds every harmonic |m| <= m_cut = 7 of M20; measured 3e-13 um^2 (rounding; L = 5: 4e-14)
  const double scale = ref.scattering_um2;
  EXPECT_NEAR(c.scatteredUp_um2, ref.scatteredUp_um2, 1e-12);
  EXPECT_NEAR(c.scatteredDown_um2, ref.scatteredDown_um2, 1e-12);
  double worst = 0;
  for (std::size_t i = 0; i < ref.harmonics.size(); ++i) {
    const int m = ref.harmonics[i];
    worst = std::max({worst, std::abs(perHarmonic(c, c.reflectedPerHarmonic, m) - ref.reflectedPerHarmonic[i]),
                      std::abs(perHarmonic(c, c.transmittedPerHarmonic, m) - ref.transmittedPerHarmonic[i])});
  }
  std::printf("V22a circle: largest deviation per harmonic %.1e (sigma %.6f)\n", worst, scale);
  EXPECT_LT(worst, 1e-12);
}

// Selection rules: in one class of P = 1 the Laurent blocks of a square couple only harmonics that differ by a
// multiple of 4 (the others vanish to rounding); a rotation of the incidence azimuth by 2 pi / 4 leaves the cross
// sections unchanged (1e-12) and rotates the far field.
TEST(V22d, SelectionRulesOfASquare) {
  cyl::SpectrumOptions o;
  o.kmax = 4;
  o.Rs = 2.0;
  o.branchPoints = {1.0, 1.5};
  const auto grid = cyl::spectralGrid(o);
  fg::geometry::Region sq;
  sq.outer = fg::geometry::Contour2D::polygon({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}});
  sq.material = "p";
  sq.normalize(1e-9);
  const cyl::HarmonicMap map{0, 1, 4};
  const auto Q = cyl::coupledQ(grid, map, {{sq}, {2.1025}, 1.0});
  const Eigen::Index n = grid.k.size();
  double coupled = 0, forbidden = 0;
  for (int b = 0; b < map.blocks(); ++b)
    for (int bp = 0; bp < map.blocks(); ++bp) {
      if (b == bp) continue;
      double& t = (b - bp) % 4 == 0 ? coupled : forbidden;
      t = std::max(t, Q.Qyy.block(b * n, bp * n, n, n).cwiseAbs().maxCoeff());
    }
  std::printf("V22d: largest block coupling l - l' in 4Z %.3e, otherwise %.1e\n", coupled, forbidden);
  EXPECT_GT(coupled, 1e-3);
  EXPECT_LT(forbidden, 1e-14 * coupled);

  const std::string sym = R"("factorization": "laurent", "symmetry": {"rotation": 4, "L": 1}, "far_field": {"n_phi": 8})";
  const auto a = solveConfig(pillar(sym, kOblique, kSquare)).cylindrical;
  std::string rotated = kOblique;
  rotated.replace(rotated.find(R"("psi_deg": 0)"), 12, R"("psi_deg": 90)");
  const auto b = solveConfig(pillar(sym, rotated, kSquare)).cylindrical;
  std::printf("V22d: sigma up %.15f / %.15f, down %.15f / %.15f (psi 0 / 90), L = %d\n", a.scatteredUp_um2,
              b.scatteredUp_um2, a.scatteredDown_um2, b.scatteredDown_um2, a.classL);
  EXPECT_NEAR(a.scatteredUp_um2, b.scatteredUp_um2, 1e-12 * a.scattering_um2);
  EXPECT_NEAR(a.scatteredDown_um2, b.scatteredDown_um2, 1e-12 * a.scattering_um2);
  ASSERT_EQ(a.farField.size(), b.farField.size());
  double ffMax = 0, ffDev = 0;
  for (std::size_t s = 0; s < a.farField.size(); ++s) {
    const auto &fa = a.farField[s], &fb = b.farField[s];
    const std::size_t np = fa.phi_deg.size();
    for (std::size_t t = 0; t < fa.theta_deg.size(); ++t)
      for (std::size_t l = 0; l < np; ++l) {  // phi_l + 90 deg = phi_{l + np/4}
        const std::size_t i = t * np + l, j = t * np + (l + np / 4) % np;
        ffMax = std::max(ffMax, fa.s[i] + fa.p[i]);
        ffDev = std::max({ffDev, std::abs(fa.s[i] - fb.s[j]), std::abs(fa.p[i] - fb.p[j])});
      }
  }
  std::printf("V22d: far field rotated by 90 deg, deviation %.1e of %.3e\n", ffDev, ffMax);
  EXPECT_LT(ffDev, 1e-11 * ffMax);
}

// Mirror symmetry phi -> -phi of the square: circular_plus (class 1) and circular_minus (class -1) go over into each
// other, harmonic m of the one equals -m of the other (1e-12).
TEST(V22e, MirrorSymmetryMapsTheClasses) {
  const std::string sym = R"("factorization": "laurent", "symmetry": {"rotation": 4, "L": 2})";
  auto beam = [](const char* pol) {
    return std::string(R"({"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.25, "polarization": ")") + pol +
           R"("})";
  };
  const auto p = solveConfig(pillar(sym, beam("circular_plus"), kSquare)).cylindrical;
  const auto m = solveConfig(pillar(sym, beam("circular_minus"), kSquare)).cylindrical;
  std::printf("V22e: R = %.15f / %.15f, T = %.15f / %.15f (circular plus / minus), L = %d\n", p.reflected, m.reflected,
              p.transmitted, m.transmitted, p.classL);
  EXPECT_NEAR(p.reflected, m.reflected, 1e-12);
  EXPECT_NEAR(p.transmitted, m.transmitted, 1e-12);
  ASSERT_EQ(p.harmonics.size(), m.harmonics.size());
  double dev = 0, coupled = 0;
  for (std::size_t i = 0; i < p.harmonics.size(); ++i) {
    const int h = p.harmonics[i];
    dev = std::max({dev, std::abs(p.reflectedPerHarmonic[i] - perHarmonic(m, m.reflectedPerHarmonic, -h)),
                    std::abs(p.transmittedPerHarmonic[i] - perHarmonic(m, m.transmittedPerHarmonic, -h))});
    if (h != 1) coupled += p.reflectedPerHarmonic[i] + p.transmittedPerHarmonic[i];
  }
  std::printf("V22e: harmonics m <-> -m deviation %.1e; power scattered into m != 1: %.3e\n", dev, coupled);
  EXPECT_LT(dev, 1e-12);
  EXPECT_GT(coupled, 1e-6);  // the square does couple (m = -3, 5, ...)
}

// Configuration (Plan sec:cylP "CAD und Konfiguration"): polygon2d / contour2d need symmetry.rotation; a cross
// section that is not P-fold symmetric is an error naming the order found; invalid symmetry blocks are rejected.
TEST(V22Config, SymmetryChecksAndErrors) {
  auto message = [](const std::string& json) {
    try {
      (void)solveConfig(json);
    } catch (const std::exception& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  EXPECT_NE(message(pillar(R"("factorization": "laurent")", kBeam, kSquare)).find("rotation"), std::string::npos);
  const char* rect = R"({"type": "polygon2d", "thickness_um": 0.3, "background": "air",
      "inclusions": [{"material": "p", "rectangle_um": {"center": [0, 0], "size": [0.5, 0.4]}}]})";
  const std::string e = message(pillar(R"("symmetry": {"rotation": 4})", kBeam, rect));
  EXPECT_NE(e.find("not 4-fold"), std::string::npos) << e;
  EXPECT_NE(e.find("largest order found: 2"), std::string::npos) << e;
  const char* offAxis = R"({"type": "polygon2d", "thickness_um": 0.3, "background": "air",
      "inclusions": [{"material": "p", "rectangle_um": {"center": [0.1, 0], "size": [0.5, 0.5]}}]})";
  EXPECT_NE(message(pillar(R"("symmetry": {"rotation": 4})", kBeam, offAxis)).find("largest order found: 1"),
            std::string::npos);
  EXPECT_NE(message(pillar(R"("symmetry": {"rotation": 0})", kBeam, kSquare)).find("rotation must be"),
            std::string::npos);
  EXPECT_NE(message(pillar(R"("symmetry": {"rotation": 4, "M": 2})", kBeam, kSquare)).find("unknown key"),
            std::string::npos);
  EXPECT_NE(message(pillar(R"("symmetry": "helical")", kBeam, kSquare)).find("axial"), std::string::npos);
}

TEST(V22b, EccentricFiberCoreHasTheModesOfTheCentredCore) {
  // core eps = 4, radius a = 2 (normalized, as V21b), cladding eps = 1; Laurent, k_max = 3, one grid for all cases
  // (R_s = a + 0.6); class m0 = 1 of P = 1 holds HE11 at m = +1 and m = -1
  const double eps1 = 4, eps2 = 1, a = 2;
  cyl::SpectrumOptions o;
  o.kmax = 3;
  o.Rs = a + 0.6;
  o.branchPoints = {1.0, 2.0};
  const auto grid = cyl::spectralGrid(o);
  const auto centred = guidedBetas(grid.k, cyl::cylQ(grid, 1, {{a}, {eps1}, eps2}, cyl::CylFactorization::Laurent), eps2);
  ASSERT_FALSE(centred.empty());
  const double he11 = centred.front();
  std::printf("V22b: N_k = %ld, centred HE11 beta = %.10f\n", static_cast<long>(grid.k.size()), he11);
  struct Case {
    double x0;
    int L;
    double tol;
  };
  // measured: 8.4e-5, 3.0e-9, 6.8e-7, 1.6e-10 (prototype: 8e-5, 4e-8, 7e-7, 5e-8); acceptance 1e-7 at L = 7
  for (const Case c : {Case{0.3, 3, 2e-4}, Case{0.3, 5, 1e-8}, Case{0.6, 5, 2e-6}, Case{0.6, 7, 1e-7}}) {
    fg::geometry::Region core;
    core.outer = fg::geometry::Contour2D::circle({c.x0, 0.0}, a);
    core.material = "core";
    core.normalize(1e-9);
    const cyl::HarmonicMap map{1, 1, c.L};
    const auto Q = cyl::coupledQ(grid, map, {{core}, {eps1}, eps2});
    const auto b = guidedBetas(cyl::tiledK(grid.k, map.blocks()), Q, eps2);
    ASSERT_GE(b.size(), 2u) << c.x0 << " " << c.L;
    const double dev = std::max(std::abs(b[0] - he11), std::abs(b[1] - he11));
    std::printf("V22b: x0 = %.1f, L = %d: HE11 pair %.10f %.10f, deviation %.1e\n", c.x0, c.L, b[0], b[1], dev);
    EXPECT_LT(dev, c.tol) << c.x0 << " " << c.L;
  }
}
