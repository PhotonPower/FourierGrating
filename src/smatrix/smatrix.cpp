#include "fg/smatrix/smatrix.hpp"

#include "fg/core/error.hpp"
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

SMatrix star(const SMatrix& A, const SMatrix& B) {
  FG_REQUIRE(A.nh() == B.nh(), "S-matrix size mismatch in star product");
  const int n = A.nh();
  const MatC Id = MatC::Identity(n, n);
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

}  // namespace fg
