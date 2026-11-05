#include "fg/solvers/resident.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"

namespace fg::resident {

using linalg::Backend;
using linalg::DeviceMatrix;

namespace {

/// [A | B] (same number of rows)
DeviceMatrix hcat(Backend& d, const DeviceMatrix& A, const DeviceMatrix& B) {
  DeviceMatrix C = d.zeros(A.rows(), A.cols() + B.cols());
  d.setBlock(C, 0, 0, A);
  d.setBlock(C, 0, A.cols(), B);
  return C;
}

/// [A; B] (same number of columns)
DeviceMatrix vcat(Backend& d, const DeviceMatrix& A, const DeviceMatrix& B) {
  DeviceMatrix C = d.zeros(A.rows() + B.rows(), A.cols());
  d.setBlock(C, 0, 0, A);
  d.setBlock(C, A.rows(), 0, B);
  return C;
}

/// Eigenvector columns needed by the classification are downloaded one by one (only modes with
/// |Re lambda| <= tau, i.e. propagating ones).
std::vector<Eigen::Index> order(Backend& d, const linalg::DeviceEig& es, double fluxSign,
                                bool block) {
  return modeOrder(es.values, [&](Eigen::Index k) {
    const VecC w = d.download(d.columns(es.vectors, {k}));
    return block ? blockModeFlux(w, fluxSign) : modeFlux(w);
  });
}

LayerModesD reduced(Backend& d, const MatC& M, double eta = 1e-10) {
  const Eigen::Index n = M.rows() / 2;  // n = 2 n_f
  const DeviceMatrix P = d.upload(M.topRightCorner(n, n)), R = d.upload(M.bottomLeftCorner(n, n));
  const linalg::DeviceEig es = d.eigResident(d.multiply(P, R));
  LayerModesD m;
  m.lu.resize(n);
  m.ld.resize(n);
  VecC iu(n), id(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    cplx kap = std::sqrt(-es.values(k));  // kappa with Im >= 0 (Plan eq:beta_branch)
    if (kap.imag() < 0 || (kap.imag() == 0 && kap.real() < 0)) kap = -kap;
    if (std::abs(kap) < eta) kap += cplx(0, eta);
    m.lu(k) = I * kap;
    m.ld(k) = -I * kap;
    iu(k) = 1.0 / m.lu(k);
    id(k) = 1.0 / m.ld(k);
  }
  const DeviceMatrix Ru = d.multiply(R, es.vectors);
  m.Wu = vcat(d, es.vectors, d.scaleCols(Ru, iu));
  m.Wd = vcat(d, es.vectors, d.scaleCols(Ru, id));
  return m;
}

/// Classified eigenmodes of M (full 4 n_f system, or one TE/TM block with block = true).
LayerModesD classified(Backend& d, const MatC& M, bool block, double fluxSign) {
  const Eigen::Index N = M.rows(), n = N / 2;
  const linalg::DeviceEig es = d.eigResident(d.upload(M));
  const std::vector<Eigen::Index> idx = order(d, es, fluxSign, block);
  LayerModesD m;
  m.Wu = d.columns(es.vectors, std::vector<Eigen::Index>(idx.begin(), idx.begin() + n));
  m.Wd = d.columns(es.vectors, std::vector<Eigen::Index>(idx.begin() + n, idx.end()));
  m.lu.resize(n);
  m.ld.resize(n);
  for (Eigen::Index k = 0; k < n; ++k) {
    m.lu(k) = es.values(idx[static_cast<std::size_t>(k)]);
    m.ld(k) = es.values(idx[static_cast<std::size_t>(n + k)]);
  }
  if (!block) {
    const double tau = 1e-9 * es.values.cwiseAbs().maxCoeff();
    FG_REQUIRE(m.lu.real().maxCoeff() <= m.ld.real().minCoeff() + tau,
               "mode classification inconsistent (up/down eigenvalues overlap)");
  }
  return m;
}

}  // namespace

SMatrixD upload(Backend& d, const SMatrix& S) {
  return {d.upload(S.S11), d.upload(S.S12), d.upload(S.S21), d.upload(S.S22)};
}

SMatrix download(Backend& d, const SMatrixD& S) {
  SMatrix out;
  out.S11 = d.download(S.S11);
  out.S12 = d.download(S.S12);
  out.S21 = d.download(S.S21);
  out.S22 = d.download(S.S22);
  return out;
}

SMatrixD star(Backend& d, const SMatrixD& A, const SMatrixD& B) {
  FG_PROF_SCOPE("resident.star");
  const Eigen::Index n = A.S11.rows();
  FG_REQUIRE(B.S11.rows() == n, "S-matrix size mismatch in star product");
  const DeviceMatrix Id = d.identity(n);
  // as fg::star with the device backend (M10): T11^-1 [A11 | A12 B22] and [B21 A11 | B22] in one
  // solve each
  const DeviceMatrix x1 =
      d.luSolve(d.multiply(A.S12, B.S21, -1.0, &Id, 1.0), hcat(d, A.S11, d.multiply(A.S12, B.S22)));
  const DeviceMatrix x2 =
      d.luSolve(d.multiply(B.S21, A.S12, -1.0, &Id, 1.0), hcat(d, d.multiply(B.S21, A.S11), B.S22));
  SMatrixD C;
  C.S11 = d.multiply(B.S11, d.block(x1, 0, 0, n, n));
  C.S12 = d.multiply(B.S11, d.block(x1, 0, n, n, n), 1.0, &B.S12, 1.0);
  C.S21 = d.multiply(A.S22, d.block(x2, 0, 0, n, n), 1.0, &A.S21, 1.0);
  C.S22 = d.multiply(A.S22, d.block(x2, 0, n, n, n));
  return C;
}

LayerModesD upload(Backend& d, const LayerModes& m) {
  return {d.upload(m.Wu), d.upload(m.Wd), m.lu, m.ld};
}

LayerModesD layerModes(Backend& d, const MatC& M, std::optional<Polarization> pol) {
  FG_PROF_SCOPE("resident.modes");
  if (isBlockAntiDiagonal(M)) return reduced(d, M);
  if (pol) return classified(d, M, true, *pol == Polarization::TE ? 1.0 : -1.0);
  const Eigen::Index n = M.rows() / 4;
  const int A[2] = {0, 3}, B[2] = {1, 2};  // (Ex, Hz) and (Ez, Hx)
  double cross = 0;
  for (int i : A)
    for (int j : B) {
      cross = std::max(cross, M.block(i * n, j * n, n, n).cwiseAbs().maxCoeff());
      cross = std::max(cross, M.block(j * n, i * n, n, n).cwiseAbs().maxCoeff());
    }
  if (cross > 0.0) return classified(d, M, false, 0.0);
  // TE/TM-decoupled: eigenmodes per 2n block, embedded into the 4n state vector (as layerModesAuto)
  LayerModesD out;
  out.Wu = d.zeros(4 * n, 2 * n);
  out.Wd = d.zeros(4 * n, 2 * n);
  out.lu.resize(2 * n);
  out.ld.resize(2 * n);
  Eigen::Index col = 0;
  for (const int* set : {A, B}) {
    MatC Mb(2 * n, 2 * n);
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
        Mb.block(a * n, b * n, n, n) = M.block(set[a] * n, set[b] * n, n, n);
    const LayerModesD bm = classified(d, Mb, true, set[0] == 0 ? -1.0 : 1.0);
    for (int a = 0; a < 2; ++a) {
      d.setBlock(out.Wu, set[a] * n, col, d.block(bm.Wu, a * n, 0, n, n));
      d.setBlock(out.Wd, set[a] * n, col, d.block(bm.Wd, a * n, 0, n, n));
    }
    out.lu.segment(col, n) = bm.lu;
    out.ld.segment(col, n) = bm.ld;
    col += n;
  }
  return out;
}

SMatrixD layerSMatrix(Backend& d, const LayerModesD& m, double h, const DeviceMatrix& refPsiInv) {
  FG_PROF_SCOPE("resident.layer_smatrix");
  const Eigen::Index n = m.lu.size();
  const DeviceMatrix Cu = d.multiply(refPsiInv, m.Wu), Cd = d.multiply(refPsiInv, m.Wd);
  const DeviceMatrix Au = d.block(Cu, 0, 0, n, n), Bu = d.block(Cu, n, 0, n, n);
  const DeviceMatrix Ad = d.block(Cd, 0, 0, n, n), Bd = d.block(Cd, n, 0, n, n);
  const VecC Xu = (m.lu * h).array().exp().matrix();   // |.| <= 1
  const VecC Xd = (-m.ld * h).array().exp().matrix();  // |.| <= 1
  DeviceMatrix G = d.zeros(2 * n, 2 * n), H = d.zeros(2 * n, 2 * n);
  d.setBlock(G, 0, 0, Bu);
  d.setBlock(G, 0, n, d.scaleCols(Bd, Xd));
  d.setBlock(G, n, 0, d.scaleCols(Au, Xu));
  d.setBlock(G, n, n, Ad);
  d.setBlock(H, 0, 0, d.scaleCols(Bu, Xu));
  d.setBlock(H, 0, n, Bd);
  d.setBlock(H, n, 0, Au);
  d.setBlock(H, n, n, d.scaleCols(Ad, Xd));
  // S = H G^-1  ->  S^T = G^-T H^T
  const DeviceMatrix GT = d.add(1.0, G, true, 0.0, {}, false),
                     HT = d.add(1.0, H, true, 0.0, {}, false);
  const DeviceMatrix S = d.add(1.0, d.luSolve(GT, HT), true, 0.0, {}, false);
  return {d.block(S, 0, 0, n, n), d.block(S, 0, n, n, n), d.block(S, n, 0, n, n),
          d.block(S, n, n, n, n)};
}

}  // namespace fg::resident
