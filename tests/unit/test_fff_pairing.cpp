// +-lambda pairing of slice eigenvalues (Plan, criterion in "Varianten der isotropen FFF"):
// a y-invariant slice of an isotropic reciprocal material has eigenvalues in pairs (lambda, -lambda).
#include <gtest/gtest.h>

#include <cmath>

#include "fg/fourier/coefficients.hpp"
#include "fg/geometry/profile.hpp"
#include "fg/physics/qmatrix.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;

namespace {
double pairingDefect(int variant) {
  fg::Profile1DLayer p;
  p.depth_um = 0.25;
  const fg::IndexMap map(10, 0);
  const int Px = 2 * map.Nx(), n = map.size();
  const double d = 1.0, k0 = 2 * fg::PI / 0.6328, K = 2 * fg::PI / (k0 * d);
  VecC al(n), ga = VecC::Zero(n);
  for (int r = 0; r < n; ++r) al(r) = std::sin(20 * fg::PI / 180) + map.n(r) * K;
  const auto nf = fg::normalField(p, d, map);
  std::vector<double> br;
  std::vector<bool> below;
  fg::profileIntervals(p, 0.125, d, br, below);
  const cplx em = cplx(0.2, 3.4) * cplx(0.2, 3.4);
  std::vector<cplx> ve, vi;
  for (bool b : below) {
    ve.push_back(b ? em : cplx(1.0));
    vi.push_back(1.0 / ve.back());
  }
  const auto ce = fg::fourier::piecewiseConstant1D(br, ve, d, Px);
  const auto ci = fg::fourier::piecewiseConstant1D(br, vi, d, Px);
  fg::SliceCoefficients sc;
  sc.eps = fg::fourier::CoefficientGrid(Px, 0);
  sc.invEps = fg::fourier::CoefficientGrid(Px, 0);
  for (int q = -Px; q <= Px; ++q) {
    sc.eps(q, 0) = ce[static_cast<std::size_t>(q + Px)];
    sc.invEps(q, 0) = ci[static_cast<std::size_t>(q + Px)];
  }
  sc.normal = &nf;
  const fg::QBlocks Q = variant == 0   ? fg::qFFFIsotropic(sc, map, false)
                        : variant == 1 ? fg::qFFFIsotropic(sc, map, true)
                        : variant == 2 ? fg::qFFFIsotropicRigorous(sc, map)
                                       : fg::qLaurent(sc, map);
  const MatC M = fg::buildPropagationMatrix(Q, al, ga);
  const VecC l = Eigen::ComplexEigenSolver<MatC>(M, false).eigenvalues();
  double worst = 0;
  for (Eigen::Index i = 0; i < l.size(); ++i) {
    double best = 1e300;
    for (Eigen::Index j = 0; j < l.size(); ++j) best = std::min(best, std::abs(l(i) + l(j)));
    worst = std::max(worst, best / std::max(1.0, std::abs(l(i))));
  }
  return worst;
}
}  // namespace

TEST(FFFPairing, LaurentAndSymmetricPreservePairing) {
  EXPECT_LT(pairingDefect(3), 1e-8);  // Laurent
  EXPECT_LT(pairingDefect(1), 1e-8);  // FFF_sym
}

TEST(FFFPairing, BookAndRigorousVariantsDoNotPair) {
  // Structural property only, not a quality criterion: FFF_rig is the most accurate variant for
  // the gold sinus (V13) although it does not pair the eigenvalues (Plan, "Varianten der FFF").
  EXPECT_GT(pairingDefect(0), 1e-2);  // FFF_iso (book IV.32)
  EXPECT_GT(pairingDefect(2), 1e-2);  // FFF_rig
}
