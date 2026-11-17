// Sectioned1DLayer with analytic line-section providers (Plan section "M23: Schichttyp Sectioned1DLayer",
// tests without Nurbforge): vertical walls = lamellar1d (Li rule), trapezoid = profile1d with Laurent,
// FFF convergence, energy, Chebyshev interpolation, anisotropic Laurent, coordinate transformation.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "../../src/solvers/layer_model.hpp"
#include "fg/core/error.hpp"
#include "fg/geometry/line_section.hpp"
#include "fg/solvers/solver.hpp"

using namespace fg;
using geometry::Contour2D;
using geometry::ProfileBody;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double d = 0.3, h = 0.1;
const double w = std::tan(15 * pi / 180);  // 15 degree walls

std::shared_ptr<const std::vector<ProfileBody>> bodies(const std::vector<geometry::Point2>& poly) {
  return std::make_shared<const std::vector<ProfileBody>>(
      std::vector<ProfileBody>{{{Contour2D::polygon(poly)}, "si", "Line"}});
}

std::vector<geometry::Point2> trapezoid() {
  return {{0.0866, 0.0}, {0.2134, 0.0}, {0.2134 - h * w, h}, {0.0866 + h * w, h}};
}

Stack base() {
  Stack st;
  st.dx_um = d;
  st.materials["air"] = Material::fromIndex("air", 1.0);
  st.materials["glass"] = Material::fromIndex("glass", 1.5);
  st.materials["si"] = Material::fromIndex("si", cplx(2.0, 0.0));
  st.superstrate = "air";
  st.substrate = "glass";
  return st;
}

Stack sectioned(const std::vector<geometry::Point2>& poly, Factorization f, int cheb = 0) {
  Stack st = base();
  Sectioned1DLayer l;
  l.provider = std::make_shared<ProfileLineProvider>(bodies(poly), 0.0, h, d, "air", geometry::LineSectionOptions{});
  l.background = "air";
  l.materials = {"air", "si"};
  l.chebyshevNodes = cheb;
  st.layers.push_back({l, f});
  return st;
}

Stack profile(Factorization f) {
  Stack st = base();
  Profile1DLayer p;
  p.depth_um = h;
  p.shape = Profile1DLayer::Shape::PiecewiseLinear;
  const auto t = trapezoid();
  p.xs_um = {0.0, t[0][0], t[3][0], t[2][0], t[1][0]};
  p.ys_um = {0.0, 0.0, h, h, 0.0};
  p.below = "si";
  p.above = "air";
  st.layers.push_back({p, f});
  return st;
}

Result run(const Stack& st, SolverKind k, int N, int slices = 40, double tol = 1e-8) {
  SolverOptions o;
  o.kind = k;
  o.Nx = N;
  o.rcwaSlices = slices;
  o.rkTol = tol;
  const Wavelength wl{0.633, 20 * pi / 180, 30 * pi / 180, 35 * pi / 180};  // conical
  return makeSolver(o)->solve(st, wl);
}

double maxDiff(const Result& a, const Result& b) {
  double m = 0;
  EXPECT_EQ(a.reflected.size(), b.reflected.size());
  EXPECT_EQ(a.transmitted.size(), b.transmitted.size());
  for (std::size_t k = 0; k < std::min(a.reflected.size(), b.reflected.size()); ++k)
    m = std::max(m, std::abs(a.reflected[k].efficiency - b.reflected[k].efficiency));
  for (std::size_t k = 0; k < std::min(a.transmitted.size(), b.transmitted.size()); ++k)
    m = std::max(m, std::abs(a.transmitted[k].efficiency - b.transmitted[k].efficiency));
  return m;
}

double qDiff(const QBlocks& a, const QBlocks& b) {
  double m = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (a(i, j).size()) m = std::max(m, (a(i, j) - b(i, j)).cwiseAbs().maxCoeff());
  return m;
}

double energy(const Result& r) {
  double s = 0;
  for (const auto& o : r.reflected) s += o.efficiency;
  for (const auto& o : r.transmitted) s += o.efficiency;
  return s;
}

}  // namespace

TEST(Sectioned1D, VerticalWallsEqualLamellar) {
  const std::vector<geometry::Point2> box{{0.1, 0.0}, {0.2, 0.0}, {0.2, h}, {0.1, h}};
  const Stack s = sectioned(box, Factorization::Auto);
  Stack l = base();
  l.layers.push_back({Lamellar1DLayer{h, {0.0, 0.1, 0.2, d}, {"air", "si", "air"}}, Factorization::Auto});
  const IndexMap map(8, 0);
  const auto ms = detail::buildLayerModel(s, 0, map, 40), ml = detail::buildLayerModel(l, 0, map, 40);
  for (double y : {0.01, 0.05, 0.09}) EXPECT_LT(qDiff(ms->Q(y), ml->Q(y)), 1e-13) << y;
  EXPECT_LT(maxDiff(run(s, SolverKind::RCWA, 8), run(l, SolverKind::RCWA, 8)), 1e-12);
  EXPECT_LT(maxDiff(run(s, SolverKind::RK, 6), run(l, SolverKind::RK, 6)), 1e-8);  // RK tol 1e-8 vs exact
}

TEST(Sectioned1D, TrapezoidLaurentEqualsProfile1D) {
  const Stack s = sectioned(trapezoid(), Factorization::Laurent), p = profile(Factorization::Laurent);
  const IndexMap map(8, 0);
  const auto ms = detail::buildLayerModel(s, 0, map, 40), mp = detail::buildLayerModel(p, 0, map, 40);
  for (double y : {0.003, 0.025, 0.05, 0.075, 0.097}) EXPECT_LT(qDiff(ms->Q(y), mp->Q(y)), 1e-13) << y;
  EXPECT_LT(maxDiff(run(s, SolverKind::RK, 6), run(p, SolverKind::RK, 6)), 1e-10);
  EXPECT_LT(maxDiff(run(s, SolverKind::RCWA, 8), run(p, SolverKind::RCWA, 8)), 1e-12);
}

TEST(Sectioned1D, FFFConvergesToTheProfileLimit) {
  // different normal fields inside the line (Plan section "M23: Normalenfeld"), same limit
  std::vector<double> diff;
  for (int N : {5, 20}) {
    const Result a = run(sectioned(trapezoid(), Factorization::FFF_sym), SolverKind::RCWA, N);
    const Result b = run(profile(Factorization::FFF_sym), SolverKind::RCWA, N);
    diff.push_back(maxDiff(a, b));
    std::printf("N = %d: |sectioned1d - profile1d| = %.3e\n", N, diff.back());
  }
  EXPECT_LT(diff.back(), diff.front());
  EXPECT_LT(diff.back(), 1e-4);
}

TEST(Sectioned1D, EnergyAndChebyshev) {
  const Stack s = sectioned(trapezoid(), Factorization::Auto);
  const Result r = run(s, SolverKind::RK, 6);
  EXPECT_NEAR(energy(r), 1.0, 1e-8);
  const Result c = run(sectioned(trapezoid(), Factorization::Auto, 16), SolverKind::RK, 6);
  EXPECT_LT(maxDiff(r, c), 1e-8);
  EXPECT_NEAR(energy(run(s, SolverKind::RCWA, 10)), 1.0, 1e-12);
}

TEST(Sectioned1D, AnisotropicLaurentEqualsProfile1D) {
  Mat3C e = Mat3C::Zero();
  e(0, 0) = 4.0;
  e(1, 1) = 3.0;
  e(2, 2) = 2.5;
  e(0, 1) = e(1, 0) = 0.3;
  Stack s = sectioned(trapezoid(), Factorization::Laurent), p = profile(Factorization::Laurent);
  s.materials["si"] = p.materials["si"] = Material::fromTensor("si", e);
  EXPECT_LT(maxDiff(run(s, SolverKind::RCWA, 6), run(p, SolverKind::RCWA, 6)), 1e-12);
  // default FFF_B (not exactly energy conserving, as for profile1d): the defect falls with N, and FFF_B
  // approaches Laurent's limit
  Stack b = sectioned(trapezoid(), Factorization::Auto);
  b.materials["si"] = Material::fromTensor("si", e);
  const double e6 = std::abs(energy(run(b, SolverKind::RCWA, 6)) - 1),
               e24 = std::abs(energy(run(b, SolverKind::RCWA, 24)) - 1);
  std::printf("FFF_B energy defect N = 6: %.3e, N = 24: %.3e\n", e6, e24);
  EXPECT_LT(e24, e6);
  EXPECT_LT(e6, 1e-4);
}

TEST(Sectioned1D, CoordinateTransformMatchesProfile1D) {
  Stack s = sectioned(trapezoid(), Factorization::Laurent), p = profile(Factorization::Laurent);
  for (Stack* st : {&s, &p}) {
    st->transform.x.type = AxisTransformSpec::Type::ASR;
    st->transform.x.eta = 0.6;
    st->transform.x.breakpoints_um = {0.0866, 0.2134};
  }
  EXPECT_LT(maxDiff(run(s, SolverKind::RCWA, 8), run(p, SolverKind::RCWA, 8)), 1e-12);
}

TEST(Sectioned1D, RejectsTwoDimensionalStacksAndSignedFactorizations) {
  Stack s = sectioned(trapezoid(), Factorization::Auto);
  s.dz_um = 0.3;
  EXPECT_THROW((void)detail::buildLayerModel(s, 0, IndexMap(3, 3), 10), Error);
  EXPECT_THROW((void)detail::buildLayerModel(sectioned(trapezoid(), Factorization::FFF_rig), 0, IndexMap(3, 0), 10),
               Error);
}
