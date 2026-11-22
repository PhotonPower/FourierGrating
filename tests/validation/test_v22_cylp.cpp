// V22 (M21, Plan section "Tests und Abnahme (M21)"): P-fold rotational symmetry with coupled harmonics.
// V22b: guided modes of an eccentric fiber core (P = 1) equal those of the centred core (Plan "Befund Prototyp").
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "fg/cyl/coupled.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/cyl/spectrum.hpp"

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

}  // namespace

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
