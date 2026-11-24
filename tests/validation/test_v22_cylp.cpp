// V22 (M21, Plan section "Tests und Abnahme (M21)"): P-fold rotational symmetry with coupled harmonics.
// V22b: guided modes of an eccentric fiber core (P = 1) equal those of the centred core (Plan "Befund Prototyp").
// V22a: rotationally symmetric structures through the class path equal M20; V22c: square and hexagonal pillar
// against the Cartesian path with PML (M18); V22d: selection rules; V22e: mirror symmetry (classes m0 and -m0).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

#include <nlohmann/json.hpp>

#include "fg/adjoint/gradient.hpp"
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

bool serverTests() {
  const char* v = std::getenv("FG_SERVER_TESTS");
  return v && std::string(v) == "1";
}

/// JSON vertices of the pillar cross section centred at (c, c): square of side 0.6 um or regular hexagon of
/// circumradius 0.37 um (about the same area), 17 digits (the symmetry check works at 1e-9 of the radius).
std::string vertices(bool hexagon, double c) {
  std::string v = "[";
  const int n = hexagon ? 6 : 4;
  const double R = hexagon ? 0.37 : 0.3 * std::sqrt(2.0), a0 = hexagon ? 0.0 : 0.25 * fg::PI;
  for (int k = 0; k < n; ++k) {
    char b[96];
    const double a = a0 + 2 * fg::PI * k / n;
    std::snprintf(b, sizeof b, "%s[%.17g, %.17g]", k ? ", " : "", c + R * std::cos(a), c + R * std::sin(a));
    v += b;
  }
  return v + "]";
}

/// V22c pillar (height 0.3 um on glass, lambda = 1 um, normal incidence): material "p" is n = 1.45 or silver
/// (eps = -40 + 4i, roughly at 1 um).
std::string materials(bool silver) {
  return std::string(R"("materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "p": )") +
         (silver ? R"({"eps": [-40.0, 4.0]})" : R"({"n": 1.45})") + "},";
}

/// Threads for the server runs (Eigen's products with OpenMP).
int blasThreads() { return std::max(1, static_cast<int>(std::thread::hardware_concurrency())); }

/// Class path: P = 4 (square) or 6 (hexagon), Laurent, kappa_max, L; scattering and absorption cross sections in um^2
/// (the absorption from the volume integral, 0 for the dielectric pillar).
std::pair<double, double> v22cCylindrical(bool hexagon, bool silver, double kappa, int L) {
  const std::string head = materials(silver);
  const std::string json = R"({"threads": {"blas": )" + std::to_string(blasThreads()) +
                           R"(}, "cylindrical": {"factorization": "laurent", "symmetry": {"rotation": )" +
                           (hexagon ? "6" : "4") + R"(, "L": )" + std::to_string(L) +
                           R"(}, "spectrum": {"kappa_max": )" + std::to_string(kappa) + "}}, " + head +
                           R"( "superstrate": "air", "substrate": "glass",
      "incidence": {"type": "plane_wave", "lambda_um": 1.0, "phi_deg": 0, "psi_deg": 0, "delta_deg": 0},
      "layers": [{"type": "polygon2d", "thickness_um": 0.3, "background": "air",
                  "inclusions": [{"material": "p", "vertices_um": )" + vertices(hexagon, 0.0) + "}]}]}";
  const auto c = solveConfig(json).cylindrical;
  EXPECT_TRUE(c.present && c.planeWave);
  return {c.scattering_um2, c.absorption_um2};
}

/// Cartesian cell 2 x 2 um with PML (0.3 um) in x and z, truncation N, Laurent: up + down + lateral (um^2) and the
/// absorbed power in the box (a convergence indicator for the dielectric pillar).
std::pair<double, double> v22cCartesian(bool hexagon, bool silver, int N) {
  const std::string head = materials(silver);
  const std::string json = R"({"threads": {"blas": )" + std::to_string(blasThreads()) +
                           R"(}, "dimension": "2D", "solver": "rcwa", "truncation": {"Nx": )" + std::to_string(N) +
                           R"(, "Nz": )" + std::to_string(N) + R"(}, "period_um": {"x": 2.0, "z": 2.0}, )" + head +
                           R"( "superstrate": "air", "substrate": "glass",
      "incidence": {"lambda_um": 1.0, "phi_deg": 0, "psi_deg": 0, "delta_deg": 0},
      "aperiodic": {"pml": {"x": {"width_um": 0.3, "beta": 1.0}, "z": {"width_um": 0.3, "beta": 1.0}}},
      "layers": [{"type": "polygon2d", "thickness_um": 0.3, "background": "air", "factorization": "laurent",
                  "inclusions": [{"material": "p", "vertices_um": )" + vertices(hexagon, 1.0) + "}]}]}";
  const fg::Result r = solveConfig(json);
  EXPECT_TRUE(r.scattering.present);
  return {r.scattering.up + r.scattering.down + r.scattering.lateral, r.scattering.absorbed};
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

// Square pillar, coarse: class path kappa_max = 4, L = 1 (0.30063) against M18 N = 6 (0.30314). Both converge from
// above (class path L = 1, kappa_max 4/6/8: 0.30063/0.29951/0.29919, kappa_max 6, L = 2: 0.29852; M18 N = 6/8:
// 0.30314/0.30206 with |absorbed| up to 2e-3), the acceptance at 1e-3 runs on the target system (below).
TEST(V22c, SquarePillarCoarse) {
  const double cyl = v22cCylindrical(false, false, 4, 1).first;
  const auto [cart, absorbed] = v22cCartesian(false, false, 6);
  std::printf("V22c square, coarse: class path %.6f | M18 %.6f (absorbed %.1e)\n", cyl, cart, absorbed);
  EXPECT_NEAR(cyl, cart, 5e-3);
}

// Square (P = 4) and hexagonal (P = 6) dielectric pillar against M18 within its convergence (Plan: 1e-3, both
// Laurent): |class path - M18(N = 14)| <= max(1e-3, |M18(14) - M18(10)|); the class path at kappa_max = 8, L = 2 and
// 3. Silver pillar as a convergence study (output only; FFF stage 2). Target system: FG_SERVER_TESTS=1.
TEST(V22c, PillarsAgainstCartesianPml) {
  if (!serverTests()) GTEST_SKIP() << "target system only (FG_SERVER_TESTS=1)";
  for (const bool hexagon : {false, true}) {
    const char* name = hexagon ? "hexagon" : "square";
    const double c2 = v22cCylindrical(hexagon, false, 8, 2).first, c3 = v22cCylindrical(hexagon, false, 8, 3).first;
    const auto m10 = v22cCartesian(hexagon, false, 10), m14 = v22cCartesian(hexagon, false, 14);
    std::printf("V22c %s: class path L = 2/3: %.6f / %.6f | M18 N = 10/14: %.6f / %.6f (absorbed %.1e / %.1e)\n",
                name, c2, c3, m10.first, m14.first, m10.second, m14.second);
    EXPECT_NEAR(c3, m14.first, std::max(1e-3, std::abs(m14.first - m10.first))) << name;
    EXPECT_LT(std::abs(c3 - c2), 1e-3) << name;
  }
  for (const bool hexagon : {false, true}) {
    const char* name = hexagon ? "hexagon" : "square";
    const auto c4 = v22cCylindrical(hexagon, true, 4, 2), c8 = v22cCylindrical(hexagon, true, 8, 2);
    const auto m10 = v22cCartesian(hexagon, true, 10), m14 = v22cCartesian(hexagon, true, 14);
    std::printf("V22c silver %s: class path kappa_max = 4/8: sca %.6f / %.6f, abs %.6f / %.6f | M18 N = 10/14: "
                "up + down + lateral %.6f / %.6f, absorbed %.6f / %.6f\n",
                name, c4.first, c8.first, c4.second, c8.second, m10.first, m14.first, m10.second, m14.second);
    EXPECT_TRUE(std::isfinite(c8.first) && std::isfinite(c8.second) && std::isfinite(m14.first)) << name;
  }
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

/// The pillar of `pillar` with the lossy material n = 1.45 + 0.1i.
std::string lossy(std::string json) {
  const std::string from = R"("p": {"n": 1.45})";
  json.replace(json.find(from), from.size(), R"("p": {"eps": [2.0925, 0.29]})");
  return json;
}

// Volume integrals in the class path (Plan v2.42): Laurent keeps the energy discretely, so the rings absorption of a
// beam equals 1 - R - T and the optical theorem of a plane wave equals sigma_sca + sigma_abs up to rounding; turning
// the azimuth of incidence by 2 pi / P leaves sigma_abs unchanged (V22d).
TEST(V22Volume, LaurentKeepsTheEnergyBalance) {
  const std::string sym = R"("factorization": "laurent", "symmetry": {"rotation": 4, "L": 1})";
  const auto beam = solveConfig(lossy(pillar(sym, kBeam, kSquare))).cylindrical;
  std::printf("V22Volume: beam A = 1 - R - T = %.15f, rings %.15f\n", beam.absorbed, beam.absorbedRings);
  EXPECT_GT(beam.absorbedRings, 1e-3);
  EXPECT_NEAR(beam.absorbedRings, beam.absorbed, 1e-12);

  const auto a = solveConfig(lossy(pillar(sym, kOblique, kSquare))).cylindrical;
  std::printf("V22Volume: sca %.15f abs %.15f ext %.15f theorem %.15f um^2\n", a.scattering_um2, a.absorption_um2,
              a.extinction_um2, a.extinctionTheorem_um2);
  EXPECT_GT(a.absorption_um2, 1e-3);
  EXPECT_NEAR(a.extinctionTheorem_um2, a.extinction_um2, 1e-12 * a.extinction_um2);
  std::string rotated = kOblique;
  rotated.replace(rotated.find(R"("psi_deg": 0)"), 12, R"("psi_deg": 90)");
  const auto b = solveConfig(lossy(pillar(sym, rotated, kSquare))).cylindrical;
  std::printf("V22Volume: sigma_abs %.15f / %.15f (psi 0 / 90)\n", a.absorption_um2, b.absorption_um2);
  EXPECT_NEAR(a.absorption_um2, b.absorption_um2, 1e-12 * a.absorption_um2);
}

// Rotationally symmetric lossy layers through the class path equal M20: rings layer (Li in r, block-wise through the
// M20 volume integral) with a beam, and the centred circle (contour2d, coupled angular coefficients) against the rings
// layer with Laurent for an oblique plane wave (sigma_abs and the optical theorem).
TEST(V22Volume, RotationallySymmetricLayersEqualM20) {
  const auto ref = solveConfig(lossy(pillar(R"("factorization": "li_r")", kBeam, kRings))).cylindrical;
  const auto cls =
      solveConfig(lossy(pillar(R"("factorization": "li_r", "symmetry": {"rotation": 3, "L": 1})", kBeam, kRings)))
          .cylindrical;
  std::printf("V22Volume: rings absorption %.15f (M20) / %.15f (P = 3)\n", ref.absorbedRings, cls.absorbedRings);
  EXPECT_GT(ref.absorbedRings, 1e-3);
  EXPECT_NEAR(cls.absorbedRings, ref.absorbedRings, 1e-12);

  const auto rings = solveConfig(lossy(pillar(R"("factorization": "laurent")", kOblique, kRings))).cylindrical;
  const auto circle = solveConfig(lossy(pillar(R"("factorization": "laurent", "symmetry": {"rotation": 3, "L": 2})",
                                               kOblique, kCircle)))
                          .cylindrical;
  std::printf("V22Volume: abs %.15f / %.15f, theorem %.15f / %.15f um^2 (rings M20 / circle P = 3)\n",
              rings.absorption_um2, circle.absorption_um2, rings.extinctionTheorem_um2, circle.extinctionTheorem_um2);
  EXPECT_GT(rings.absorption_um2, 1e-3);
  EXPECT_NEAR(circle.absorption_um2, rings.absorption_um2, 1e-12);
  EXPECT_NEAR(circle.extinctionTheorem_um2, rings.extinctionTheorem_um2, 1e-12);
}

/// Near field of a configuration with the cuts of `block` (V21 style).
fg::NearFieldResult nearFieldOf(const std::string& json, const std::string& block) {
  const auto p = json.rfind('}');
  fg::io::RunConfig cfg = fg::io::parseConfig(json.substr(0, p) + R"(, "nearfield": {"enabled": true, )" + block + "}}");
  return fg::makeSolver(cfg.solver)->nearField(cfg.stack, cfg.scan.front(), cfg.field);
}

/// Largest |E| and largest deviation of E and H between two near fields (both cuts), and the largest material
/// deviation.
struct FieldDeviation {
  double scale = 0, dev = 0, material = 0;
};
FieldDeviation deviation(const fg::NearFieldResult& a, const fg::NearFieldResult& b) {
  FieldDeviation d;
  for (const bool vertical : {true, false}) {
    const fg::FieldCut &x = vertical ? a.vertical : a.horizontal, &y = vertical ? b.vertical : b.horizontal;
    EXPECT_EQ(x.E[0].size(), y.E[0].size());
    for (int c = 0; c < 3; ++c)
      for (std::size_t i = 0; i < x.E[c].size(); ++i) {
        d.scale = std::max(d.scale, std::abs(x.E[c][i]));
        d.dev = std::max({d.dev, std::abs(x.E[c][i] - y.E[c][i]), std::abs(x.H[c][i] - y.H[c][i])});
      }
    for (std::size_t i = 0; i < x.material.size(); ++i) d.material = std::max(d.material, std::abs(x.material[i] - y.material[i]));
  }
  return d;
}

/// Cuts within about R_s (beam: R_s = 2 w0 = 0.5 um, plane wave: R_s = 0.3 um, the pillar radius); beyond R_s the
/// synthesis loses accuracy (Plan v2.25) and amplifies the tiny differences of the spectra. No point lies on the wall
/// r = 0.3 um (rings count it to the inside, Region::contains to the outside).
std::string cuts(double extent) {
  return R"("horizontal_y_um": 0.15, "nx": 15, "nz": 15, "ny": 7, "y_range_um": [-0.2, 0.5], "extent_um": )" +
         std::to_string(extent);
}

// Near field in the class path (Plan v2.42): the class spectra split into their harmonics and synthesized as in M20.
// Rings layer (Li in r, FFF reconstruction per block) with a beam equals the M20 near field (1e-11); the centred
// circle (contour2d, Laurent, quadrature) with an oblique plane wave equals the rings layer (Laurent, Lommel integrals)
// to 1e-9 (measured 3e-10); the material map of the square follows its cross section.
TEST(V22NearField, ClassPathEqualsM20) {
  const auto ref = nearFieldOf(pillar(R"("factorization": "li_r")", kBeam, kRings), cuts(0.5));
  const auto cls =
      nearFieldOf(pillar(R"("factorization": "li_r", "symmetry": {"rotation": 3, "L": 1})", kBeam, kRings), cuts(0.5));
  const FieldDeviation r = deviation(ref, cls);
  std::printf("V22NearField: rings, beam: deviation %.1e of max |E| %.3f, material %.1e\n", r.dev, r.scale, r.material);
  EXPECT_LT(r.dev, 1e-11 * r.scale);
  EXPECT_EQ(r.material, 0.0);

  const auto rings = nearFieldOf(pillar(R"("factorization": "laurent")", kOblique, kRings), cuts(0.29));
  const auto circle =
      nearFieldOf(pillar(R"("factorization": "laurent", "symmetry": {"rotation": 3, "L": 2})", kOblique, kCircle), cuts(0.29));
  const FieldDeviation c = deviation(rings, circle);
  std::printf("V22NearField: circle, plane wave: deviation %.1e of max |E| %.3f, material %.1e\n", c.dev, c.scale,
              c.material);
  EXPECT_LT(c.dev, 1e-9 * c.scale);  // quadrature against Lommel integrals: 3e-10
  EXPECT_EQ(c.material, 0.0);

  // square 0.5 x 0.5 um: (0.2, 0.2) inside, (0.3, 0.1) outside at the height of the layer
  const auto sq = nearFieldOf(pillar(R"("factorization": "laurent", "symmetry": {"rotation": 4, "L": 1})", kBeam, kSquare),
                              R"("horizontal_y_um": 0.15, "nx": 7, "nz": 7, "ny": 1, "y_range_um": [0.15, 0.15], "extent_um": 0.3)");
  const auto& h = sq.horizontal;
  auto material = [&](int ix, int iz) { return h.material[static_cast<std::size_t>(iz * 7 + ix)]; };
  EXPECT_DOUBLE_EQ(material(5, 5), 1.45 * 1.45);  // (0.2, 0.2)
  EXPECT_DOUBLE_EQ(material(6, 4), 1.0);          // (0.3, 0.1)
  EXPECT_DOUBLE_EQ(material(3, 3), 1.45 * 1.45);  // axis
  for (int c3 = 0; c3 < 3; ++c3)
    for (const cplx v : h.E[c3]) EXPECT_TRUE(std::isfinite(std::abs(v)));
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

/// Lossy square pillar (P = 4) over a ring layer in a film, Li in r (the rings per block, the square Laurent), x-polarized
/// beam; `scattered`: in the scattered-field formulation.
std::string scatteredStack(bool scattered) {
  return R"({)" + std::string(scattered ? R"("scattered_field": true, )" : "") +
         R"("cylindrical": {"factorization": "li_r", "symmetry": {"rotation": 4, "L": 1}, "spectrum": {"kappa_max": 2}},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "f": {"n": 1.38}, "p": {"eps": [2.0925, 0.29]}, "q": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.25, "focus_y_um": 0.3, "polarization": "x"},
    "layers": [{"type": "rings", "thickness_um": 0.15, "radii_um": [0.2], "materials": ["q", "f"]},
               {"type": "polygon2d", "thickness_um": 0.25, "background": "air",
                "inclusions": [{"material": "p", "rectangle_um": {"center": [0, 0], "size": [0.4, 0.4]}}]}]})";
}

// Beam as scattered field in the class path (Plan v2.46, V21e for M21): the background stack (structured layers ->
// their outer material) is diagonal in the harmonics and solved per class on the grid; scattered field + background
// = total field for R, T, A_rings per harmonic and the near field.
TEST(V22Scattered, BeamScatteredFieldEqualsTotalField) {
  const fg::CylindricalResult t = solveConfig(scatteredStack(false)).cylindrical,
                              s = solveConfig(scatteredStack(true)).cylindrical;
  std::printf("V22Scattered: R = %.12f (%.1e)  T = %.12f (%.1e)  A_rings = %.12f (%.1e)\n", t.reflected,
              s.reflected - t.reflected, t.transmitted, s.transmitted - t.transmitted, t.absorbedRings,
              s.absorbedRings - t.absorbedRings);
  EXPECT_GT(t.absorbedRings, 1e-3);
  EXPECT_NEAR(s.reflected, t.reflected, 1e-10);
  EXPECT_NEAR(s.transmitted, t.transmitted, 1e-10);
  EXPECT_NEAR(s.absorbedRings, t.absorbedRings, 1e-10);
  ASSERT_EQ(s.harmonics, t.harmonics);
  for (std::size_t i = 0; i < t.reflectedPerHarmonic.size(); ++i) {
    EXPECT_NEAR(s.reflectedPerHarmonic[i], t.reflectedPerHarmonic[i], 1e-10) << t.harmonics[i];
    EXPECT_NEAR(s.transmittedPerHarmonic[i], t.transmittedPerHarmonic[i], 1e-10) << t.harmonics[i];
  }
  const std::string nf =
      R"("horizontal_y_um": 0.27, "nx": 9, "nz": 9, "ny": 7, "y_range_um": [-0.1, 0.5], "extent_um": 0.45)";
  const FieldDeviation d = deviation(nearFieldOf(scatteredStack(false), nf), nearFieldOf(scatteredStack(true), nf));
  std::printf("V22Scattered near field: deviation %.1e of max |E| %.3f\n", d.dev, d.scale);
  EXPECT_LT(d.dev, 1e-10 * d.scale);
  EXPECT_EQ(d.material, 0.0);
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

// V22g (Plan "Tests und Abnahme (M21)"): gradients in the class path (Gaussian beam, total field, RCWA) against central
// difference quotients, criterion of M19 1e-6 max(1, |fd|). The spectral grid stays fixed (substrate of the highest
// index, R_s = 2 w0 larger than the cross sections). P = 1: the vertices of an irregular contour2d pentagon (a
// regular one would make the difference quotients of the discrete map unreliable: its critical radii coincide and
// split under the step), its material, the thickness and the background; P = 2: semi-axis and angle of a centred
// ellipse (both keep the 2-fold symmetry) above a rings layer (Li in r, block-diagonal over the class)
namespace {

nlohmann::json v22gCase(int P) {
  auto j = nlohmann::json::parse(R"({
    "cylindrical": {"spectrum": {"kappa_max": 2, "dk_factor": 1}, "factorization": "li_r", "symmetry": {"L": 1}},
    "materials": {"air": {"n": 1.0}, "sub": {"n": 2.0}, "core": {"n": [1.7, 0.02]}, "film": {"n": 1.6}},
    "superstrate": "air", "substrate": "sub",
    "incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.2, "focus_y_um": 0.1,
                  "polarization": "x"},
    "layers": []})");
  j["cylindrical"]["symmetry"]["rotation"] = P;
  if (P == 1) {
    j["layers"] = nlohmann::json::array(
        {{{"type", "homogeneous"}, {"thickness_um", 0.08}, {"material", "film"}},
         {{"type", "contour2d"}, {"thickness_um", 0.15}, {"background", "film"},
          {"regions", {{{"material", "core"},
                        {"polygon_um", {{0.15, 0.03}, {0.05, 0.16}, {-0.12, 0.11}, {-0.1, -0.09}, {0.07, -0.13}}}}}}}});
  } else {
    j["layers"] = nlohmann::json::array(
        {{{"type", "rings"}, {"thickness_um", 0.08}, {"radii_um", {0.12}}, {"materials", {"film", "air"}}},
         {{"type", "contour2d"}, {"thickness_um", 0.15}, {"background", "air"},
          {"regions",
           {{{"material", "core"}, {"ellipse", {{"center_um", {0, 0}}, {"radii_um", {0.18, 0.1}}, {"angle_deg", 20}}}}}}}});
  }
  return j;
}

double v22gObjective(const fg::CylindricalResult& r) {
  double t1 = 0;
  for (std::size_t i = 0; i < r.harmonics.size(); ++i)
    if (r.harmonics[i] == 1) t1 = r.transmittedPerHarmonic[i];
  return r.reflected + 0.5 * r.transmitted + 2 * t1;
}

void checkV22g(int P, const std::vector<std::string>& params) {
  fg::Objective obj;
  obj.cylindrical = {{fg::CylRef{'R', true, 0}, 1.0}, {fg::CylRef{'T', true, 0}, 0.5}, {fg::CylRef{'T', false, 1}, 2.0}};
  const nlohmann::json j = v22gCase(P);
  const fg::GradientResult g = fg::gradient(j.dump(), obj, params);
  EXPECT_NEAR(g.value, v22gObjective(solveConfig(j.dump()).cylindrical), 1e-12) << P;
  for (std::size_t q = 0; q < params.size(); ++q) {
    const nlohmann::json::json_pointer p(params[q]);
    const double v = j.at(p).get<double>(), h = 1e-5 * std::max(0.1, std::abs(v));
    nlohmann::json a = j, b = j;
    a[p] = v + h;
    b[p] = v - h;
    const double fd =
        (v22gObjective(solveConfig(a.dump()).cylindrical) - v22gObjective(solveConfig(b.dump()).cylindrical)) / (2 * h);
    std::printf("V22g P = %d %-40s adjoint %+.10e  fd %+.10e  diff %.1e\n", P, params[q].c_str(), g.grad[q], fd,
                std::abs(g.grad[q] - fd));
    EXPECT_NEAR(g.grad[q], fd, 1e-6 * std::max(1.0, std::abs(fd))) << P << " " << params[q];
  }
}

}  // namespace

TEST(V22g, ContourVerticesOfOneClass) {
  checkV22g(1, {"/layers/1/regions/0/polygon_um/1/0", "/layers/1/regions/0/polygon_um/1/1",
                "/layers/1/regions/0/polygon_um/3/0", "/materials/core/n/0", "/materials/core/n/1",
                "/materials/film/n", "/layers/1/thickness_um"});
}

TEST(V22g, EllipseAboveRingsWithTwofoldSymmetry) {
  checkV22g(2, {"/layers/1/regions/0/ellipse/radii_um/0", "/layers/1/regions/0/ellipse/angle_deg",
                "/layers/0/radii_um/0", "/materials/film/n", "/materials/core/n/1"});
  // a parameter that breaks the 2-fold symmetry is rejected (the class path would ignore the coupling)
  fg::Objective obj;
  obj.cylindrical = {{fg::CylRef{'R', true, 0}, 1.0}};
  EXPECT_THROW(fg::gradient(v22gCase(2).dump(), obj, {"/layers/1/regions/0/ellipse/center_um/0"}), fg::Error);
}
