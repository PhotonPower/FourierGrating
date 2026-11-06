#include "fg/solvers/rcwa.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/linalg/backend.hpp"

namespace fg {

bool isBlockAntiDiagonal(const MatC& M, double tol) {
  const Eigen::Index n = M.rows() / 2;
  const double ref = M.cwiseAbs().maxCoeff();
  return M.topLeftCorner(n, n).cwiseAbs().maxCoeff() <= tol * ref &&
         M.bottomRightCorner(n, n).cwiseAbs().maxCoeff() <= tol * ref;
}

/// General eigenproblem: device backend for large matrices (M10), Eigen otherwise.
static linalg::EigResult eigen(const MatC& A, const char* what) {
  FG_PROF_SCOPE("linalg.eig");
  if (linalg::offload(A.rows())) return linalg::device().eig(A);
  Eigen::ComplexEigenSolver<MatC> es(A, true);
  FG_REQUIRE(es.info() == Eigen::Success, what);
  return {es.eigenvalues(), es.eigenvectors()};
}

LayerModes layerModesReduced(const MatC& M, double eta) {
  const Eigen::Index n = M.rows() / 2;  // n = 2 n_f
  const MatC P = M.topRightCorner(n, n), R = M.bottomLeftCorner(n, n);
  const linalg::EigResult es = eigen(P * R, "eigen decomposition failed (RCWA reduced)");
  LayerModes m;
  m.Wu.resize(2 * n, n);
  m.Wd.resize(2 * n, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    cplx kap = std::sqrt(-es.values(k));  // kappa with Im >= 0 (Plan eq:beta_branch)
    if (kap.imag() < 0 || (kap.imag() == 0 && kap.real() < 0)) kap = -kap;
    if (std::abs(kap) < eta) kap += cplx(0, eta);
    const VecC u = es.vectors.col(k);
    const VecC Ru = R * u;
    const cplx lup = I * kap, ldn = -I * kap;
    m.lu(k) = lup;
    m.ld(k) = ldn;
    m.Wu.col(k) << u, Ru / lup;
    m.Wd.col(k) << u, Ru / ldn;
  }
  return m;
}

double classificationTolerance() {
  // relative to max |lambda|; eigenvalues in single precision (M16) are accurate to ~eps32 max |lambda|
  return linalg::precision() == linalg::Precision::FP64 ? 1e-9 : 1e-4;
}

std::vector<Eigen::Index> modeOrder(const VecC& lam, const std::function<double(Eigen::Index)>& flux) {
  // Classification (Plan eq:classify): field ~ exp(lambda y); up-going modes decay towards +y.
  // Robust variant: sort by a score (Re lambda, or -sign(S_y) * tau for |Re lambda| <= tau)
  // and take the lower half as up-going. This always yields n_h modes per direction.
  const Eigen::Index N = lam.size();
  const double tau = classificationTolerance() * lam.cwiseAbs().maxCoeff();
  std::vector<double> score(static_cast<std::size_t>(N));
  for (Eigen::Index k = 0; k < N; ++k) {
    const double re = lam(k).real();
    if (std::abs(re) > tau) {
      score[static_cast<std::size_t>(k)] = re;
    } else {
      score[static_cast<std::size_t>(k)] = flux(k) > 0 ? -0.5 * tau : 0.5 * tau;
    }
  }
  std::vector<Eigen::Index> idx(static_cast<std::size_t>(N));
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](Eigen::Index i, Eigen::Index j) {
    return score[static_cast<std::size_t>(i)] < score[static_cast<std::size_t>(j)];
  });
  return idx;
}

double modeFlux(const VecC& w) {
  const Eigen::Index nf = w.size() / 4;
  double sy = 0;
  for (Eigen::Index r = 0; r < nf; ++r)
    sy += 0.5 * std::real(w(nf + r) * std::conj(w(2 * nf + r)) - w(r) * std::conj(w(3 * nf + r)));
  return sy;
}

double blockModeFlux(const VecC& w, double fluxSign) {
  const Eigen::Index n = w.size() / 2;
  double sy = 0;
  for (Eigen::Index r = 0; r < n; ++r) sy += fluxSign * 0.5 * std::real(w(r) * std::conj(w(n + r)));
  return sy;
}

LayerModes layerModesGeneral(const MatC& M) {
  const linalg::EigResult es = eigen(M, "eigen decomposition failed (RCWA general)");
  const Eigen::Index N = M.rows(), n = N / 2;
  const VecC& lam = es.values;
  const double tau = classificationTolerance() * lam.cwiseAbs().maxCoeff();
  const std::vector<Eigen::Index> idx = modeOrder(lam, [&](Eigen::Index k) { return modeFlux(es.vectors.col(k)); });
  LayerModes m;
  m.Wu.resize(N, n);
  m.Wd.resize(N, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    const Eigen::Index iu = idx[static_cast<std::size_t>(k)], id = idx[static_cast<std::size_t>(n + k)];
    m.Wu.col(k) = es.vectors.col(iu);
    m.lu(k) = lam(iu);
    m.Wd.col(k) = es.vectors.col(id);
    m.ld(k) = lam(id);
  }
  FG_REQUIRE(m.lu.real().maxCoeff() <= m.ld.real().minCoeff() + 2 * tau,  // flux-sorted modes: |Re| <= tau
             "mode classification inconsistent (up/down eigenvalues overlap)");
  return m;
}

/// Eigenmodes of one decoupled polarization block on (E, H) = 2 n rows, classified by Re lambda
/// and, for |Re lambda| <= tau, by the flux sign * 1/2 Re(E H*) (TE: +Ez Hx*, TM: -Ex Hz*).
static LayerModes blockModes(const MatC& Mb, double fluxSign) {
  const Eigen::Index n = Mb.rows() / 2;
  const linalg::EigResult es = eigen(Mb, "eigen decomposition failed (RCWA TE/TM block)");
  const VecC& lam = es.values;
  const std::vector<Eigen::Index> idx =
      modeOrder(lam, [&](Eigen::Index k) { return blockModeFlux(es.vectors.col(k), fluxSign); });
  LayerModes m;
  m.Wu.resize(2 * n, n);
  m.Wd.resize(2 * n, n);
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    const Eigen::Index iu = idx[static_cast<std::size_t>(k)], id = idx[static_cast<std::size_t>(n + k)];
    m.Wu.col(k) = es.vectors.col(iu);
    m.Wd.col(k) = es.vectors.col(id);
    m.lu(k) = lam(iu);
    m.ld(k) = lam(id);
  }
  return m;
}

LayerModes layerModesPolarization(const MatC& M, Polarization p) {
  if (isBlockAntiDiagonal(M)) return layerModesReduced(M);
  return blockModes(M, p == Polarization::TE ? 1.0 : -1.0);
}

LayerModes layerModesAuto(const MatC& M) {
  if (isBlockAntiDiagonal(M)) return layerModesReduced(M);
  const Eigen::Index n = M.rows() / 4;
  const int A[2] = {0, 3}, B[2] = {1, 2};  // (Ex, Hz) and (Ez, Hx)
  double cross = 0;
  for (int i : A)
    for (int j : B) {
      cross = std::max(cross, M.block(i * n, j * n, n, n).cwiseAbs().maxCoeff());
      cross = std::max(cross, M.block(j * n, i * n, n, n).cwiseAbs().maxCoeff());
    }
  if (cross > 0.0) return layerModesGeneral(M);
  // TE/TM-decoupled: eigenmodes per 2n block, embedded into the 4n state vector
  LayerModes out;
  out.Wu = MatC::Zero(4 * n, 2 * n);
  out.Wd = MatC::Zero(4 * n, 2 * n);
  out.lu.resize(2 * n);
  out.ld.resize(2 * n);
  Eigen::Index col = 0;
  for (const int* set : {A, B}) {
    MatC Mb(2 * n, 2 * n);
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) Mb.block(a * n, b * n, n, n) = M.block(set[a] * n, set[b] * n, n, n);
    const LayerModes bm = blockModes(Mb, set[0] == 0 ? -1.0 : 1.0);  // (Ex, Hz): -Ex Hz*; (Ez, Hx): +Ez Hx*
    for (int a = 0; a < 2; ++a) {
      out.Wu.block(set[a] * n, col, n, n) = bm.Wu.middleRows(a * n, n);
      out.Wd.block(set[a] * n, col, n, n) = bm.Wd.middleRows(a * n, n);
    }
    out.lu.segment(col, n) = bm.lu;
    out.ld.segment(col, n) = bm.ld;
    col += n;
  }
  return out;
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
  FG_PROF_SCOPE("rcwa.layer_smatrix");
  const Eigen::Index n = m.lu.size();
  const MatC Cu = linalg::product(ref.PsiInv, m.Wu), Cd = linalg::product(ref.PsiInv, m.Wd);
  const MatC Au = Cu.topRows(n), Bu = Cu.bottomRows(n);
  const MatC Ad = Cd.topRows(n), Bd = Cd.bottomRows(n);
  const VecC Xu = (m.lu * h).array().exp().matrix();     // |.| <= 1
  const VecC Xd = (-m.ld * h).array().exp().matrix();    // |.| <= 1
  MatC G(2 * n, 2 * n), H(2 * n, 2 * n);
  G << Bu, Bd * Xd.asDiagonal(), Au * Xu.asDiagonal(), Ad;
  H << Bu * Xu.asDiagonal(), Bd, Au, Ad * Xd.asDiagonal();
  // S = H G^-1  ->  S^T = G^-T H^T
  const MatC S = linalg::offload(2 * n) ? MatC(linalg::device().solve(G.transpose(), H.transpose()).transpose())
                                        : MatC(G.transpose().partialPivLu().solve(H.transpose()).transpose());
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
