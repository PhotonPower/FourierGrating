#include "fg/smatrix/smatrix.hpp"

#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"
#include "fg/linalg/backend.hpp"

namespace fg {

SMatrix SMatrix::identity(int nh) {
  SMatrix S;
  S.S11 = MatC::Identity(nh, nh);
  S.S22 = MatC::Identity(nh, nh);
  S.S12 = MatC::Zero(nh, nh);
  S.S21 = MatC::Zero(nh, nh);
  return S;
}

MatC SMatrix::full() const {
  const int n = nh();
  MatC F(2 * n, 2 * n);
  F << S11, S12, S21, S22;
  return F;
}

SMatrix fromTransfer(const MatC& T) {
  FG_PROF_SCOPE("smatrix.from_transfer");
  FG_REQUIRE(T.rows() == T.cols() && T.rows() % 2 == 0, "transfer matrix must be 2nh x 2nh");
  const Eigen::Index n = T.rows() / 2;
  const MatC T11 = T.topLeftCorner(n, n), T12 = T.topRightCorner(n, n);
  const MatC T21 = T.bottomLeftCorner(n, n), T22 = T.bottomRightCorner(n, n);
  SMatrix S;
  if (linalg::offload(n)) {  // device backend (M10): T11^-1 [1 | T12] in one solve
    linalg::Backend& d = linalg::device();
    MatC rhs(n, 2 * n);
    rhs << MatC::Identity(n, n), T12;
    const MatC X = d.solve(T11, rhs);
    S.S22 = X.leftCols(n);
    S.S21 = -X.rightCols(n);
    S.S12 = d.gemm(T21, S.S22);
    S.S11 = T22 + d.gemm(T21, S.S21);
    return S;
  }
  Eigen::PartialPivLU<MatC> lu(T11);
  S.S22 = lu.inverse();               // T11^-1 (needed explicitly as a block)
  S.S21 = -lu.solve(T12);             // -T11^-1 T12
  S.S12 = T21 * S.S22;                // T21 T11^-1
  S.S11 = T22 + T21 * S.S21;          // T22 - T21 T11^-1 T12
  return S;
}

namespace {

/// The 2 x 2 group block of an n_h = 2G block and its inverse (scatter).
MatC groupBlock(const MatC& M, Eigen::Index g, Eigen::Index G) {
  MatC b(2, 2);
  b << M(g, g), M(g, G + g), M(G + g, g), M(G + g, G + g);
  return b;
}
void setGroupBlock(MatC& M, const MatC& b, Eigen::Index g, Eigen::Index G) {
  M(g, g) = b(0, 0);
  M(g, G + g) = b(0, 1);
  M(G + g, g) = b(1, 0);
  M(G + g, G + g) = b(1, 1);
}
SMatrix zeroSMatrix(Eigen::Index n) {
  SMatrix S;
  S.S11 = S.S12 = S.S21 = S.S22 = MatC::Zero(n, n);
  return S;
}

/// X Y with X or Y grouped over G (rows / columns {g, G + g}): O(n^2 m) -> O(n m).
MatC groupedMul(const MatC& X, Eigen::Index gx, const MatC& Y, Eigen::Index gy) {
  if (gx > 0) {
    MatC C(X.rows(), Y.cols());
    for (Eigen::Index g = 0; g < gx; ++g) {
      const Eigen::Index i[2] = {g, gx + g};
      for (int r = 0; r < 2; ++r) C.row(i[r]) = X(i[r], i[0]) * Y.row(i[0]) + X(i[r], i[1]) * Y.row(i[1]);
    }
    return C;
  }
  if (gy > 0) {
    MatC C(X.rows(), Y.cols());
    for (Eigen::Index g = 0; g < gy; ++g) {
      const Eigen::Index i[2] = {g, gy + g};
      for (int c = 0; c < 2; ++c) C.col(i[c]) = X.col(i[0]) * Y(i[0], i[c]) + X.col(i[1]) * Y(i[1], i[c]);
    }
    return C;
  }
  return X * Y;
}

/// Star product with one grouped factor (the formulas of star; with srcC the sources too: srcC holds those of the
/// lower factor on entry and the result on return, srcB those of the upper one).
SMatrix mixedStar(const SMatrix& A, const SMatrix& B, const SourcedSMatrix* srcB, SourcedSMatrix* srcC);

}  // namespace

SMatrix groupedStar(const SMatrix& A, const SMatrix& B, Eigen::Index G) {
  FG_PROF_SCOPE("smatrix.star_grouped");
  FG_REQUIRE(A.nh() == 2 * G && B.nh() == 2 * G, "groupedStar: n_h = 2G");
  SMatrix C = zeroSMatrix(2 * G);
  C.groups = G;
  const MatC Id = MatC::Identity(2, 2);
  for (Eigen::Index g = 0; g < G; ++g) {
    const MatC a11 = groupBlock(A.S11, g, G), a12 = groupBlock(A.S12, g, G), a21 = groupBlock(A.S21, g, G),
               a22 = groupBlock(A.S22, g, G), b11 = groupBlock(B.S11, g, G), b12 = groupBlock(B.S12, g, G),
               b21 = groupBlock(B.S21, g, G), b22 = groupBlock(B.S22, g, G);
    const Eigen::PartialPivLU<MatC> D1(Id - a12 * b21), D2(Id - b21 * a12);  // as star
    setGroupBlock(C.S11, b11 * D1.solve(a11), g, G);
    setGroupBlock(C.S12, b12 + b11 * D1.solve(a12 * b22), g, G);
    setGroupBlock(C.S21, a21 + a22 * D2.solve(b21 * a11), g, G);
    setGroupBlock(C.S22, a22 * D2.solve(b22), g, G);
  }
  return C;
}

SMatrix groupedFromTransfer(const MatC& PsiInv, const MatC& Psi, Eigen::Index G) {
  FG_PROF_SCOPE("smatrix.from_transfer");
  FG_REQUIRE(PsiInv.rows() == 4 * G && Psi.rows() == 4 * G, "groupedFromTransfer: 4G states");
  SMatrix S = zeroSMatrix(2 * G);
  S.groups = G;
  for (Eigen::Index g = 0; g < G; ++g) {
    const Eigen::Index idx[4] = {g, G + g, 2 * G + g, 3 * G + g};
    MatC a(4, 4), b(4, 4);
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c) {
        a(r, c) = PsiInv(idx[r], idx[c]);
        b(r, c) = Psi(idx[r], idx[c]);
      }
    const MatC T = a * b;  // as fromTransfer
    const MatC T11 = T.topLeftCorner(2, 2), T12 = T.topRightCorner(2, 2), T21 = T.bottomLeftCorner(2, 2),
               T22 = T.bottomRightCorner(2, 2);
    const Eigen::PartialPivLU<MatC> lu(T11);
    const MatC S22 = lu.inverse(), S21 = -lu.solve(T12);
    setGroupBlock(S.S22, S22, g, G);
    setGroupBlock(S.S21, S21, g, G);
    setGroupBlock(S.S12, T21 * S22, g, G);
    setGroupBlock(S.S11, T22 + T21 * S21, g, G);
  }
  return S;
}

namespace {
SMatrix mixedStar(const SMatrix& A, const SMatrix& B, const SourcedSMatrix* srcB, SourcedSMatrix* srcC) {
  FG_PROF_SCOPE("smatrix.star_mixed");
  const Eigen::Index n = A.nh(), ga = A.groups, gb = B.groups;
  auto mul = [](const MatC& X, Eigen::Index gx, const MatC& Y, Eigen::Index gy) { return groupedMul(X, gx, Y, gy); };
  const MatC Id = MatC::Identity(n, n);
  const Eigen::PartialPivLU<MatC> D1(Id - mul(A.S12, ga, B.S21, gb));  // (1 - A12 B21)
  const Eigen::PartialPivLU<MatC> D2(Id - mul(B.S21, gb, A.S12, ga));  // (1 - B21 A12)
  SMatrix C;
  C.S11 = mul(B.S11, gb, D1.solve(A.S11), 0);
  C.S12 = B.S12 + mul(B.S11, gb, D1.solve(mul(A.S12, ga, B.S22, gb)), 0);
  C.S21 = A.S21 + mul(A.S22, ga, D2.solve(mul(B.S21, gb, A.S11, ga)), 0);
  C.S22 = mul(A.S22, ga, D2.solve(B.S22), 0);
  if (srcC) {  // srcC holds the sources of A (u, v), srcB those of B
    const MatC Au = srcC->u, Av = srcC->v, Bu = srcB->u, Bv = srcB->v;
    srcC->u = Bu + mul(B.S11, gb, D1.solve(Au + mul(A.S12, ga, Bv, 0)), 0);
    srcC->v = Av + mul(A.S22, ga, D2.solve(Bv + mul(B.S21, gb, Au, 0)), 0);
  }
  return C;
}
}  // namespace

SMatrix star(const SMatrix& A, const SMatrix& B) {
  FG_PROF_SCOPE("smatrix.star");
  FG_REQUIRE(A.nh() == B.nh(), "S-matrix size mismatch in star product");
  const int n = A.nh();
  const MatC Id = MatC::Identity(n, n);
  if (A.groups > 0 && A.groups == B.groups) return groupedStar(A, B, A.groups);
  if ((A.groups > 0 || B.groups > 0) && !linalg::offload(n)) return mixedStar(A, B, nullptr, nullptr);
  if (linalg::offload(n)) {  // device backend (M10)
    linalg::Backend& d = linalg::device();
    MatC r1(n, 2 * n), r2(n, 2 * n);
    r1 << A.S11, d.gemm(A.S12, B.S22);
    r2 << d.gemm(B.S21, A.S11), B.S22;
    const MatC x1 = d.solve(Id - d.gemm(A.S12, B.S21), r1);
    const MatC x2 = d.solve(Id - d.gemm(B.S21, A.S12), r2);
    SMatrix C;
    C.S11 = d.gemm(B.S11, x1.leftCols(n));
    C.S12 = B.S12 + d.gemm(B.S11, x1.rightCols(n));
    C.S21 = A.S21 + d.gemm(A.S22, x2.leftCols(n));
    C.S22 = d.gemm(A.S22, x2.rightCols(n));
    return C;
  }
  Eigen::PartialPivLU<MatC> D1(Id - A.S12 * B.S21);  // (1 - A12 B21)
  Eigen::PartialPivLU<MatC> D2(Id - B.S21 * A.S12);  // (1 - B21 A12)
  SMatrix C;
  C.S11 = B.S11 * D1.solve(A.S11);
  C.S12 = B.S12 + B.S11 * D1.solve(A.S12 * B.S22);
  C.S21 = A.S21 + A.S22 * D2.solve(B.S21 * A.S11);
  C.S22 = A.S22 * D2.solve(B.S22);
  return C;
}

VecC interfaceAmplitudes(const SMatrix& L, const SMatrix& U, const VecC& b_in, const VecC& a_in) {
  const int n = L.nh();
  const MatC Id = MatC::Identity(n, n);
  const VecC rhs = L.S11 * b_in + L.S12 * (U.S22 * a_in);
  const VecC bq = (Id - L.S12 * U.S21).partialPivLu().solve(rhs);
  const VecC aq = U.S21 * bq + U.S22 * a_in;
  VecC out(2 * n);
  out << aq, bq;
  return out;
}

SourcedSMatrix SourcedSMatrix::identity(int nh) {
  return {SMatrix::identity(nh), VecC::Zero(nh), VecC::Zero(nh)};
}

SourcedSMatrix star(const SourcedSMatrix& A, const SourcedSMatrix& B) {
  FG_PROF_SCOPE("smatrix.star_sourced");
  const int n = A.S.nh();
  const MatC Id = MatC::Identity(n, n);
  if (A.S.groups > 0 || B.S.groups > 0) {
    const bool noSources = A.u.isZero(0.0) && A.v.isZero(0.0) && B.u.isZero(0.0) && B.v.isZero(0.0);
    if (A.S.groups == B.S.groups && noSources) return {groupedStar(A.S, B.S, A.S.groups), A.u, A.v};
    SMatrix a = A.S, b = B.S;
    if (a.groups > 0 && b.groups > 0) b.groups = 0;  // sources with two grouped factors: one taken dense
    SourcedSMatrix C{SMatrix(), A.u, A.v};
    C.S = mixedStar(a, b, &B, &C);
    return C;
  }
  Eigen::PartialPivLU<MatC> D1(Id - A.S.S12 * B.S.S21);
  Eigen::PartialPivLU<MatC> D2(Id - B.S.S21 * A.S.S12);
  SourcedSMatrix C;
  C.S.S11 = B.S.S11 * D1.solve(A.S.S11);
  C.S.S12 = B.S.S12 + B.S.S11 * D1.solve(A.S.S12 * B.S.S22);
  C.S.S21 = A.S.S21 + A.S.S22 * D2.solve(B.S.S21 * A.S.S11);
  C.S.S22 = A.S.S22 * D2.solve(B.S.S22);
  C.u = B.u + B.S.S11 * D1.solve(A.u + A.S.S12 * B.v);
  C.v = A.v + A.S.S22 * D2.solve(B.v + B.S.S21 * A.u);
  return C;
}

VecC interfaceAmplitudes(const SourcedSMatrix& L, const SourcedSMatrix& U, const VecC& b_in, const VecC& a_in) {
  const int n = L.S.nh();
  const MatC Id = MatC::Identity(n, n);
  const VecC rhs = L.S.S11 * b_in + L.S.S12 * (U.S.S22 * a_in + U.v) + L.u;
  const VecC bq = (Id - L.S.S12 * U.S.S21).partialPivLu().solve(rhs);
  const VecC aq = U.S.S21 * bq + U.S.S22 * a_in + U.v;
  VecC out(2 * n);
  out << aq, bq;
  return out;
}

}  // namespace fg
