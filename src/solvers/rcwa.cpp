#include "fg/solvers/rcwa.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "fg/core/error.hpp"

namespace fg {

bool isBlockAntiDiagonal(const MatC& M, double tol) {
  const Eigen::Index n = M.rows() / 2;
  const double ref = M.cwiseAbs().maxCoeff();
  return M.topLeftCorner(n, n).cwiseAbs().maxCoeff() <= tol * ref &&
         M.bottomRightCorner(n, n).cwiseAbs().maxCoeff() <= tol * ref;
}

LayerModes layerModesReduced(const MatC& M, double eta) {
  const Eigen::Index n = M.rows() / 2;  // n = 2 n_f
  const MatC P = M.topRightCorner(n, n), R = M.bottomLeftCorner(n, n);
  Eigen::ComplexEigenSolver<MatC> es(P * R, true);
  FG_REQUIRE(es.info() == Eigen::Success, "eigen decomposition failed (RCWA reduced)");
  LayerModes m;
  m.Wu.resize(2 * n, n);
  m.Wd.resize(2 * n, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    cplx kap = std::sqrt(-es.eigenvalues()(k));  // kappa with Im >= 0 (Plan eq:beta_branch)
    if (kap.imag() < 0 || (kap.imag() == 0 && kap.real() < 0)) kap = -kap;
    if (std::abs(kap) < eta) kap += cplx(0, eta);
    const VecC u = es.eigenvectors().col(k);
    const VecC Ru = R * u;
    const cplx lup = I * kap, ldn = -I * kap;
    m.lu(k) = lup;
    m.ld(k) = ldn;
    m.Wu.col(k) << u, Ru / lup;
    m.Wd.col(k) << u, Ru / ldn;
  }
  return m;
}

LayerModes layerModesGeneral(const MatC& M) {
  Eigen::ComplexEigenSolver<MatC> es(M, true);
  FG_REQUIRE(es.info() == Eigen::Success, "eigen decomposition failed (RCWA general)");
  const Eigen::Index N = M.rows(), n = N / 2, nf = N / 4;
  const VecC& lam = es.eigenvalues();
  const double tau = 1e-9 * lam.cwiseAbs().maxCoeff();
  std::vector<Eigen::Index> up, dn;
  for (Eigen::Index k = 0; k < N; ++k) {
    const double re = lam(k).real();
    bool isUp;
    if (re < -tau) isUp = true;
    else if (re > tau) isUp = false;
    else {
      const VecC w = es.eigenvectors().col(k);
      double s = 0;
      for (Eigen::Index r = 0; r < nf; ++r)
        s += 0.5 * std::real(w(nf + r) * std::conj(w(2 * nf + r)) - w(r) * std::conj(w(3 * nf + r)));
      isUp = s > 0;
    }
    (isUp ? up : dn).push_back(k);
  }
  FG_REQUIRE(static_cast<Eigen::Index>(up.size()) == n && static_cast<Eigen::Index>(dn.size()) == n,
             "mode classification failed: " + std::to_string(up.size()) + " up / " +
                 std::to_string(dn.size()) + " down");
  LayerModes m;
  m.Wu.resize(N, n);
  m.Wd.resize(N, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    m.Wu.col(k) = es.eigenvectors().col(up[static_cast<std::size_t>(k)]);
    m.lu(k) = lam(up[static_cast<std::size_t>(k)]);
    m.Wd.col(k) = es.eigenvectors().col(dn[static_cast<std::size_t>(k)]);
    m.ld(k) = lam(dn[static_cast<std::size_t>(k)]);
  }
  return m;
}

LayerModes layerModesHomogeneous(const ModeBasis& mb) {
  const Eigen::Index N = mb.Psi.rows(), n = N / 2;
  LayerModes m;
  m.Wd = mb.Psi.leftCols(n);   // columns (a_x, a_z): down
  m.Wu = mb.Psi.rightCols(n);  // columns (b_x, b_z): up
  m.ld = mb.lambda.head(n);
  m.lu = mb.lambda.tail(n);
  return m;
}

SMatrix stableLayerSMatrix(const LayerModes& m, double h, const ModeBasis& ref,
                           RCWASliceRecord* record) {
  const Eigen::Index n = m.lu.size();
  const MatC Cu = ref.PsiInv * m.Wu, Cd = ref.PsiInv * m.Wd;
  const MatC Au = Cu.topRows(n), Bu = Cu.bottomRows(n);
  const MatC Ad = Cd.topRows(n), Bd = Cd.bottomRows(n);
  const VecC Xu = (m.lu * h).array().exp().matrix();     // |.| <= 1
  const VecC Xd = (-m.ld * h).array().exp().matrix();    // |.| <= 1
  MatC G(2 * n, 2 * n), H(2 * n, 2 * n);
  G << Bu, Bd * Xd.asDiagonal(), Au * Xu.asDiagonal(), Ad;
  H << Bu * Xu.asDiagonal(), Bd, Au, Ad * Xd.asDiagonal();
  // S = H G^-1  ->  S^T = G^-T H^T
  const MatC S = G.transpose().partialPivLu().solve(H.transpose()).transpose();
  if (record) {
    record->modes = m;
    record->h = h;
    record->G = G;
  }
  SMatrix out;
  out.S11 = S.topLeftCorner(n, n);
  out.S12 = S.topRightCorner(n, n);
  out.S21 = S.bottomLeftCorner(n, n);
  out.S22 = S.bottomRightCorner(n, n);
  return out;
}

}  // namespace fg
