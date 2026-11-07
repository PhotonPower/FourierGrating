// FourierGrating - adjoint rules of the building blocks (Plan section "M19: Konvention").
#include "fg/adjoint/primitives.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/fourier/toeplitz.hpp"

namespace fg::adjoint {

void accumulate(MatC& a, const MatC& b) {
  if (a.size() == 0) a = b;
  else a += b;
}
void accumulate(VecC& a, const VecC& b) {
  if (a.size() == 0) a = b;
  else a += b;
}

void product(cplx c, const MatC& P, const MatC& X, const MatC& Cbar, MatC* Pbar, MatC* Xbar) {
  const cplx cc = std::conj(c);
  if (Pbar) accumulate(*Pbar, cc * (Cbar * X.adjoint()));
  if (Xbar) accumulate(*Xbar, cc * (P.adjoint() * Cbar));
}

void solve(const MatC& A, const MatC& X, const MatC& Xbar, MatC* Abar, MatC* Bbar) {
  const MatC Bhat = A.adjoint().partialPivLu().solve(Xbar);
  if (Bbar) accumulate(*Bbar, Bhat);
  if (Abar) accumulate(*Abar, MatC(-Bhat * X.adjoint()));
}

void inverse(const MatC& C, const MatC& Cbar, MatC& Abar) { accumulate(Abar, MatC(-C.adjoint() * Cbar * C.adjoint())); }

MatC eig(const MatC& W, const VecC& lambda, const MatC& Wbar, const VecC& lambdabar, const std::vector<int>& cls) {
  const Eigen::Index n = W.rows();
  FG_REQUIRE(W.cols() == n && lambda.size() == n, "adjoint::eig: square W and n eigenvalues");
  const double lmax = lambda.size() ? lambda.cwiseAbs().maxCoeff() : 0.0;
  const double tau = 1e-10 * lmax, delta = 1e-12 * lmax;
  MatC K = MatC::Zero(n, n);
  if (Wbar.size()) {
    const MatC G = W.adjoint() * Wbar;
    for (Eigen::Index i = 0; i < n; ++i)
      for (Eigen::Index j = 0; j < n; ++j) {
        if (i == j) continue;
        const cplx dl = lambda(j) - lambda(i);
        const bool same = cls.empty() || cls[static_cast<std::size_t>(i)] == cls[static_cast<std::size_t>(j)];
        if (same && std::abs(dl) < tau) continue;              // degenerate cluster: subspace only
        const cplx cP = dl / (std::norm(dl) + delta * delta);   // conj(P_ij), Lorentz broadened
        K(i, j) = cP * G(i, j);
      }
  }
  if (lambdabar.size())
    for (Eigen::Index i = 0; i < n; ++i) K(i, i) += lambdabar(i);
  // Abar = W^-H K W^H
  return W.adjoint().partialPivLu().solve(MatC(K * W.adjoint()));
}

void toeplitz(const IndexMap& map, const MatC& Tbar, fourier::CoefficientGrid& fbar) {
  if (fbar.Px() < 2 * map.Nx() || fbar.Pz() < 2 * map.Nz()) fbar = fourier::CoefficientGrid(2 * map.Nx(), 2 * map.Nz());
  const int nf = map.size();
  for (int r = 0; r < nf; ++r) {
    const auto [n1, m1] = map.orders(r);
    for (int s = 0; s < nf; ++s) {
      const auto [n2, m2] = map.orders(s);
      fbar(n1 - n2, m1 - m2) += Tbar(r, s);
    }
  }
}

namespace {
void zeroIfEmpty(SMatrix& S, int n) {
  for (MatC* B : {&S.S11, &S.S12, &S.S21, &S.S22})
    if (B->size() == 0) *B = MatC::Zero(n, n);
}
}  // namespace

void star(const SMatrix& A, const SMatrix& B, const SMatrix& Cbar, SMatrix& Abar, SMatrix& Bbar) {
  const int n = A.nh();
  zeroIfEmpty(Abar, n);
  zeroIfEmpty(Bbar, n);
  const MatC Id = MatC::Identity(n, n);
  const MatC D1 = Id - A.S12 * B.S21, D2 = Id - B.S21 * A.S12;
  Eigen::PartialPivLU<MatC> L1(D1), L2(D2);
  const MatC X1 = L1.solve(A.S11), Y1 = L1.solve(MatC(A.S12 * B.S22));
  const MatC X2 = L2.solve(MatC(B.S21 * A.S11)), Y2 = L2.solve(B.S22);
  // C11 = B11 X1, C12 = B12 + B11 Y1, C21 = A21 + A22 X2, C22 = A22 Y2
  Bbar.S11 += Cbar.S11 * X1.adjoint() + Cbar.S12 * Y1.adjoint();
  Bbar.S12 += Cbar.S12;
  Abar.S21 += Cbar.S21;
  Abar.S22 += Cbar.S21 * X2.adjoint() + Cbar.S22 * Y2.adjoint();
  const MatC X1b = B.S11.adjoint() * Cbar.S11, Y1b = B.S11.adjoint() * Cbar.S12;
  const MatC X2b = A.S22.adjoint() * Cbar.S21, Y2b = A.S22.adjoint() * Cbar.S22;
  Eigen::PartialPivLU<MatC> L1h(MatC(D1.adjoint())), L2h(MatC(D2.adjoint()));
  MatC D1b = MatC::Zero(n, n), D2b = MatC::Zero(n, n);
  {  // X1 = D1^-1 A11
    const MatC R = L1h.solve(X1b);
    Abar.S11 += R;
    D1b -= R * X1.adjoint();
  }
  {  // Y1 = D1^-1 (A12 B22)
    const MatC R = L1h.solve(Y1b);
    Abar.S12 += R * B.S22.adjoint();
    Bbar.S22 += A.S12.adjoint() * R;
    D1b -= R * Y1.adjoint();
  }
  {  // X2 = D2^-1 (B21 A11)
    const MatC R = L2h.solve(X2b);
    Bbar.S21 += R * A.S11.adjoint();
    Abar.S11 += B.S21.adjoint() * R;
    D2b -= R * X2.adjoint();
  }
  {  // Y2 = D2^-1 B22
    const MatC R = L2h.solve(Y2b);
    Bbar.S22 += R;
    D2b -= R * Y2.adjoint();
  }
  // D1 = 1 - A12 B21, D2 = 1 - B21 A12
  Abar.S12 -= D1b * B.S21.adjoint();
  Bbar.S21 -= A.S12.adjoint() * D1b;
  Bbar.S21 -= D2b * A.S12.adjoint();
  Abar.S12 -= B.S21.adjoint() * D2b;
}

void fromTransfer(const MatC& T, const SMatrix& Sbar, MatC& Tbar) {
  const Eigen::Index n = T.rows() / 2;
  const MatC T11 = T.topLeftCorner(n, n), T12 = T.topRightCorner(n, n), T21 = T.bottomLeftCorner(n, n);
  Eigen::PartialPivLU<MatC> lu(T11);
  const MatC S22 = lu.inverse(), X = lu.solve(T12), S21 = -X;
  MatC T11b = MatC::Zero(n, n), T12b = MatC::Zero(n, n), T21b = MatC::Zero(n, n), T22b = MatC::Zero(n, n);
  MatC S21b = Sbar.S21.size() ? Sbar.S21 : MatC(MatC::Zero(n, n));
  MatC S22b = Sbar.S22.size() ? Sbar.S22 : MatC(MatC::Zero(n, n));
  if (Sbar.S11.size()) {  // S11 = T22 + T21 S21
    T22b += Sbar.S11;
    T21b += Sbar.S11 * S21.adjoint();
    S21b += T21.adjoint() * Sbar.S11;
  }
  if (Sbar.S12.size()) {  // S12 = T21 S22
    T21b += Sbar.S12 * S22.adjoint();
    S22b += T21.adjoint() * Sbar.S12;
  }
  {  // S21 = -T11^-1 T12
    const MatC R = MatC(T11.adjoint()).partialPivLu().solve(MatC(-S21b));
    T12b += R;
    T11b -= R * X.adjoint();
  }
  T11b -= S22.adjoint() * S22b * S22.adjoint();  // S22 = T11^-1
  MatC out(2 * n, 2 * n);
  out << T11b, T12b, T21b, T22b;
  accumulate(Tbar, out);
}

void layerSMatrix(const LayerModes& m, double h, const ModeBasis& ref, const SMatrix& Sbar, MatC& Wubar, MatC& Wdbar,
                  VecC& lubar, VecC& ldbar, double& hbar, MatC* PsiInvbar) {
  const Eigen::Index n = m.lu.size();
  const MatC Cu = ref.PsiInv * m.Wu, Cd = ref.PsiInv * m.Wd;
  const MatC Au = Cu.topRows(n), Bu = Cu.bottomRows(n), Ad = Cd.topRows(n), Bd = Cd.bottomRows(n);
  const VecC Xu = (m.lu * h).array().exp().matrix(), Xd = (-m.ld * h).array().exp().matrix();
  MatC G(2 * n, 2 * n), H(2 * n, 2 * n);
  G << Bu, Bd * Xd.asDiagonal(), Au * Xu.asDiagonal(), Ad;
  H << Bu * Xu.asDiagonal(), Bd, Au, Ad * Xd.asDiagonal();
  const MatC S = G.transpose().partialPivLu().solve(H.transpose()).transpose();  // H G^-1
  MatC Sb = MatC::Zero(2 * n, 2 * n);
  if (Sbar.S11.size()) Sb.topLeftCorner(n, n) = Sbar.S11;
  if (Sbar.S12.size()) Sb.topRightCorner(n, n) = Sbar.S12;
  if (Sbar.S21.size()) Sb.bottomLeftCorner(n, n) = Sbar.S21;
  if (Sbar.S22.size()) Sb.bottomRightCorner(n, n) = Sbar.S22;
  const MatC Hb = G.partialPivLu().solve(Sb.adjoint()).adjoint();  // Sbar G^-H = (G^-1 Sbar^H)^H
  const MatC Gb = -S.adjoint() * Hb;
  auto blk = [n](const MatC& M, int i, int j) { return MatC(M.block(i * n, j * n, n, n)); };
  const MatC G11 = blk(Gb, 0, 0), G12 = blk(Gb, 0, 1), G21 = blk(Gb, 1, 0), G22 = blk(Gb, 1, 1);
  const MatC H11 = blk(Hb, 0, 0), H12 = blk(Hb, 0, 1), H21 = blk(Hb, 1, 0), H22 = blk(Hb, 1, 1);
  const VecC cXu = Xu.conjugate(), cXd = Xd.conjugate();
  const MatC Bub = G11 + H11 * cXu.asDiagonal();
  const MatC Aub = G21 * cXu.asDiagonal() + H21;
  const MatC Bdb = H12 + G12 * cXd.asDiagonal();
  const MatC Adb = G22 + H22 * cXd.asDiagonal();
  VecC Xub(n), Xdb(n);
  for (Eigen::Index j = 0; j < n; ++j) {
    Xub(j) = Bu.col(j).dot(H11.col(j)) + Au.col(j).dot(G21.col(j));  // dot conjugates the first argument
    Xdb(j) = Bd.col(j).dot(G12.col(j)) + Ad.col(j).dot(H22.col(j));
  }
  MatC Cub(2 * n, n), Cdb(2 * n, n);
  Cub << Aub, Bub;
  Cdb << Adb, Bdb;
  accumulate(Wubar, MatC(ref.PsiInv.adjoint() * Cub));
  accumulate(Wdbar, MatC(ref.PsiInv.adjoint() * Cdb));
  if (PsiInvbar) accumulate(*PsiInvbar, MatC(Cub * m.Wu.adjoint() + Cdb * m.Wd.adjoint()));
  VecC lub(n), ldb(n);
  for (Eigen::Index j = 0; j < n; ++j) {
    lub(j) = std::conj(h * Xu(j)) * Xub(j);
    ldb(j) = std::conj(-h * Xd(j)) * Xdb(j);
    hbar += std::real(std::conj(Xub(j)) * m.lu(j) * Xu(j)) + std::real(std::conj(Xdb(j)) * (-m.ld(j)) * Xd(j));
  }
  accumulate(lubar, lub);
  accumulate(ldbar, ldb);
}

namespace {
/// Shared part of both propagation-matrix rules: Pa = alpha (or A), Pg = gamma (or G); Pab, Pgb: their
/// cotangents through the products and solves (the diagonal alpha gamma, alpha^2, gamma^2 terms are added by
/// the caller).
void propagationBackward(const QBlocks& Q, const MatC& Pa, const MatC& Pg, const MatC& Mbar, QBlocks& Qbar,
                         MatC* Pab = nullptr, MatC* Pgb = nullptr) {
  const int n = Q.nf();
  constexpr int X = 0, Y = 1, Z = 2;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (Qbar(i, j).size() == 0) Qbar(i, j) = MatC::Zero(n, n);
  Eigen::PartialPivLU<MatC> Eyy(Q(Y, Y));
  const MatC E_yx = Eyy.solve(Q(Y, X)), E_yz = Eyy.solve(Q(Y, Z)), E_g = Eyy.solve(Pg), E_a = Eyy.solve(Pa);
  auto Mb = [&](int i, int j) { return MatC(Mbar.block(i * n, j * n, n, n)); };
  MatC Eyxb = MatC::Zero(n, n), Eyzb = MatC::Zero(n, n), Egb = MatC::Zero(n, n), Eab = MatC::Zero(n, n);
  const cplx i1 = I;
  // rows 1, 2 (Plan eq:Mmu1a)
  product(-i1, Pa, E_yx, Mb(0, 0), Pab, &Eyxb);
  product(-i1, Pa, E_yz, Mb(0, 1), Pab, &Eyzb);
  product(-i1, Pa, E_g, Mb(0, 2), Pab, &Egb);
  product(i1, Pa, E_a, Mb(0, 3), Pab, &Eab);
  product(-i1, Pg, E_yx, Mb(1, 0), Pgb, &Eyxb);
  product(-i1, Pg, E_yz, Mb(1, 1), Pgb, &Eyzb);
  product(-i1, Pg, E_g, Mb(1, 2), Pgb, &Egb);
  product(i1, Pg, E_a, Mb(1, 3), Pgb, &Eab);
  // rows 3, 4 (Plan eq:Mmu1b)
  Qbar(Z, X) += std::conj(i1) * Mb(2, 0);
  product(-i1, Q(Z, Y), E_yx, Mb(2, 0), &Qbar(Z, Y), &Eyxb);
  Qbar(Z, Z) += std::conj(i1) * Mb(2, 1);
  product(-i1, Q(Z, Y), E_yz, Mb(2, 1), &Qbar(Z, Y), &Eyzb);
  product(-i1, Q(Z, Y), E_g, Mb(2, 2), &Qbar(Z, Y), &Egb);
  product(i1, Q(Z, Y), E_a, Mb(2, 3), &Qbar(Z, Y), &Eab);
  product(i1, Q(X, Y), E_yx, Mb(3, 0), &Qbar(X, Y), &Eyxb);
  Qbar(X, X) += std::conj(-i1) * Mb(3, 0);
  product(i1, Q(X, Y), E_yz, Mb(3, 1), &Qbar(X, Y), &Eyzb);
  Qbar(X, Z) += std::conj(-i1) * Mb(3, 1);
  product(i1, Q(X, Y), E_g, Mb(3, 2), &Qbar(X, Y), &Egb);
  product(-i1, Q(X, Y), E_a, Mb(3, 3), &Qbar(X, Y), &Eab);
  // E_* = Q_yy^-1 rhs
  solve(Q(Y, Y), E_yx, Eyxb, &Qbar(Y, Y), &Qbar(Y, X));
  solve(Q(Y, Y), E_yz, Eyzb, &Qbar(Y, Y), &Qbar(Y, Z));
  solve(Q(Y, Y), E_g, Egb, &Qbar(Y, Y), Pgb);
  solve(Q(Y, Y), E_a, Eab, &Qbar(Y, Y), Pab);
}
}  // namespace

void propagationMatrix(const QBlocks& Q, const VecC& alpha, const VecC& gamma, const MatC& Mbar, QBlocks& Qbar) {
  propagationBackward(Q, MatC(alpha.asDiagonal()), MatC(gamma.asDiagonal()), Mbar, Qbar);
}
void propagationMatrix(const QBlocks& Q, const MatC& A, const MatC& G, const MatC& Mbar, QBlocks& Qbar) {
  propagationBackward(Q, A, G, Mbar, Qbar);
}
void propagationMatrix(const QBlocks& Q, const MatC& A, const MatC& G, const MatC& Mbar, QBlocks& Qbar, MatC& Abar,
                       MatC& Gbar) {
  const int n = Q.nf();
  if (Abar.size() == 0) Abar = MatC::Zero(n, n);
  if (Gbar.size() == 0) Gbar = MatC::Zero(n, n);
  propagationBackward(Q, A, G, Mbar, Qbar, &Abar, &Gbar);
  // rows 3, 4: i A G, -i A A, i G G, -i G A (C = s X Y: Xbar += conj(s) Cbar Y^H, Ybar += conj(s) X^H Cbar)
  const cplx i1 = I;
  const MatC m20 = std::conj(i1) * Mbar.block(2 * n, 0, n, n), m21 = std::conj(-i1) * Mbar.block(2 * n, n, n, n);
  const MatC m30 = std::conj(i1) * Mbar.block(3 * n, 0, n, n), m31 = std::conj(-i1) * Mbar.block(3 * n, n, n, n);
  Abar += m20 * G.adjoint() + m21 * A.adjoint() + A.adjoint() * m21 + G.adjoint() * m31;
  Gbar += A.adjoint() * m20 + m30 * G.adjoint() + G.adjoint() * m30 + m31 * A.adjoint();
}
void propagationMatrix(const QBlocks& Q, const VecC& alpha, const VecC& gamma, const MatC& Mbar, QBlocks& Qbar,
                       VecC& alphabar, VecC& gammabar) {
  const int n = Q.nf();
  MatC Pab = MatC::Zero(n, n), Pgb = MatC::Zero(n, n);
  propagationBackward(Q, MatC(alpha.asDiagonal()), MatC(gamma.asDiagonal()), Mbar, Qbar, &Pab, &Pgb);
  if (alphabar.size() == 0) alphabar = VecC::Zero(n);
  if (gammabar.size() == 0) gammabar = VecC::Zero(n);
  const cplx i1 = I;
  for (int r = 0; r < n; ++r) {
    const cplx a = alpha(r), g = gamma(r);
    // diagonal terms of rows 3, 4: i alpha gamma, -i alpha^2, i gamma^2, -i alpha gamma
    const cplx m20 = Mbar(2 * n + r, r), m21 = Mbar(2 * n + r, n + r), m30 = Mbar(3 * n + r, r),
               m31 = Mbar(3 * n + r, n + r);
    alphabar(r) += Pab(r, r) + std::conj(i1 * g) * m20 + std::conj(-2.0 * i1 * a) * m21 + std::conj(-i1 * g) * m31;
    gammabar(r) += Pgb(r, r) + std::conj(i1 * a) * m20 + std::conj(2.0 * i1 * g) * m30 + std::conj(-i1 * a) * m31;
  }
}

void isotropicModes(cplx eps, const VecC& alpha, const VecC& gamma, const ModeBasis& mb, const MatC& Psibar,
                    const MatC& PsiInvbar, const VecC& betabar, cplx& epsbar, VecC& alphabar, VecC& gammabar) {
  // Psi = [1 1; -V V], PsiInv = 1/2 [1 -Vi; 1 Vi]; per order V = (p, q; -s, -p), Vi = (-p, -q; s, p) / eps with
  // p = alpha gamma / beta, q = (eps - alpha^2) / beta, s = (eps - gamma^2) / beta
  const Eigen::Index n = alpha.size();
  if (alphabar.size() == 0) alphabar = VecC::Zero(n);
  if (gammabar.size() == 0) gammabar = VecC::Zero(n);
  auto vb = [&](Eigen::Index i, Eigen::Index j) {  // Vbar
    if (Psibar.size() == 0) return cplx(0);
    return -Psibar(2 * n + i, j) + Psibar(2 * n + i, 2 * n + j);
  };
  auto vib = [&](Eigen::Index i, Eigen::Index j) {  // Vibar
    if (PsiInvbar.size() == 0) return cplx(0);
    return -0.5 * PsiInvbar(i, 2 * n + j) + 0.5 * PsiInvbar(2 * n + i, 2 * n + j);
  };
  const cplx ie = 1.0 / eps;
  for (Eigen::Index r = 0; r < n; ++r) {
    const cplx a = alpha(r), g = gamma(r), b = mb.beta(r);
    const cplx p = a * g / b, q = (eps - a * a) / b, s = (eps - g * g) / b;
    const cplx v00 = vb(r, r), v01 = vb(r, n + r), v10 = vb(n + r, r), v11 = vb(n + r, n + r);
    const cplx w00 = vib(r, r), w01 = vib(r, n + r), w10 = vib(n + r, r), w11 = vib(n + r, n + r);
    const cplx pb = v00 - v11 + std::conj(ie) * (w11 - w00), qb = v01 - std::conj(ie) * w01,
               sb = -v10 + std::conj(ie) * w10;
    // the factor 1 / eps of Vi
    const cplx de = -ie * ie;
    epsbar += std::conj(-de * p) * w00 + std::conj(-de * q) * w01 + std::conj(de * s) * w10 + std::conj(de * p) * w11;
    cplx bb = betabar.size() ? betabar(r) : cplx(0);
    bb += std::conj(-p / b) * pb + std::conj(-q / b) * qb + std::conj(-s / b) * sb;
    // beta^2 = eps - alpha^2 - gamma^2
    alphabar(r) += std::conj(g / b) * pb + std::conj(-2.0 * a / b) * qb + std::conj(-a / b) * bb;
    gammabar(r) += std::conj(a / b) * pb + std::conj(-2.0 * g / b) * sb + std::conj(-g / b) * bb;
    epsbar += std::conj(1.0 / b) * (qb + sb) + std::conj(0.5 / b) * bb;
  }
}

namespace {
MatC block(const QBlocks& Q, int i, int j, int n) { return Q(i, j).size() ? Q(i, j) : MatC(MatC::Zero(n, n)); }
}  // namespace

void qLaurent(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g) {
  FG_REQUIRE(c.isotropic, "adjoint::qLaurent: isotropic slice (anisotropic: qLaurentTensor)");
  const int n = map.size();
  toeplitz(map, MatC(block(Qbar, 0, 0, n) + block(Qbar, 1, 1, n) + block(Qbar, 2, 2, n)), g.eps);
}

void qLaurentTensor(const IndexMap& map, const QBlocks& Qbar, std::array<std::array<fourier::CoefficientGrid, 3>, 3>& g) {
  const int n = map.size();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      toeplitz(map, block(Qbar, i, j, n), g[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
}

void qLiLamellar(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g) {
  const int n = map.size();
  const MatC L = fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  MatC Tb;
  inverse(L, block(Qbar, 0, 0, n), Tb);
  toeplitz(map, Tb, g.invEps);
  toeplitz(map, MatC(block(Qbar, 1, 1, n) + block(Qbar, 2, 2, n)), g.eps);
}

void qFFFIsotropic(const SliceCoefficients& c, const IndexMap& map, bool symmetric, const QBlocks& Qbar, IsoGridBar& g) {
  const int n = map.size();
  const MatC L = fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  const MatC Delta = fourier::toeplitz(map, c.eps) - L;
  MatC Eb = MatC::Zero(n, n), Db = MatC::Zero(n, n);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const MatC Qb = block(Qbar, i, j, n);
      const MatC NN = fourier::toeplitz(map, c.normal->NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
      // Q_ij = delta_ij E - Delta NN (or -(Delta NN + NN Delta) / 2)
      if (symmetric) Db -= 0.5 * (Qb * NN.adjoint() + NN.adjoint() * Qb);
      else Db -= Qb * NN.adjoint();
      const MatC NNb = symmetric ? MatC(-0.5 * (Delta.adjoint() * Qb + Qb * Delta.adjoint())) : MatC(-Delta.adjoint() * Qb);
      toeplitz(map, NNb, g.NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
      if (i == j) Eb += Qb;
    }
  Eb += Db;  // Delta = E - L
  MatC Tb;
  inverse(L, MatC(-Db), Tb);
  toeplitz(map, Tb, g.invEps);
  toeplitz(map, Eb, g.eps);
}

void qFFFIsotropicRigorous(const SliceCoefficients& c, const IndexMap& map, const QBlocks& Qbar, IsoGridBar& g) {
  const int n = map.size();
  const MatC L = fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  std::array<MatC, 3> Nt;
  for (std::size_t i = 0; i < 3; ++i) Nt[i] = fourier::toeplitz(map, c.normal->N[i]);
  const MatC E = fourier::toeplitz(map, c.eps);
  MatC Eb = MatC::Zero(n, n), Lb = MatC::Zero(n, n);
  std::array<MatC, 3> Nb{MatC::Zero(n, n), MatC::Zero(n, n), MatC::Zero(n, n)};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) {
      const MatC Qb = block(Qbar, int(i), int(j), n);
      // Q_ij = delta_ij E - E (N_i N_j) + N_i L N_j
      Eb -= Qb * (Nt[i] * Nt[j]).adjoint();
      Lb += Nt[i].adjoint() * Qb * Nt[j].adjoint();
      if (i == j) Eb += Qb;
      Nb[i] += -E.adjoint() * Qb * Nt[j].adjoint() + Qb * (L * Nt[j]).adjoint();
      Nb[j] += -(E * Nt[i]).adjoint() * Qb + (Nt[i] * L).adjoint() * Qb;
    }
  for (std::size_t i = 0; i < 3; ++i) toeplitz(map, Nb[i], g.N[i]);
  MatC Tb;
  inverse(L, Lb, Tb);
  toeplitz(map, Tb, g.invEps);
  toeplitz(map, Eb, g.eps);
}

namespace {
MatC blockToeplitz(const std::vector<fourier::CoefficientGrid>& g, int off, const IndexMap& map) {
  const int n = map.size();
  MatC T(3 * n, 3 * n);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) T.block(i * n, j * n, n, n) = fourier::toeplitz(map, g[static_cast<std::size_t>(off + 3 * i + j)]);
  return T;
}
void blockToeplitzBar(const MatC& Tb, int off, const IndexMap& map, std::vector<fourier::CoefficientGrid>& gbar) {
  const int n = map.size();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) toeplitz(map, MatC(Tb.block(i * n, j * n, n, n)), gbar[static_cast<std::size_t>(off + 3 * i + j)]);
}
}  // namespace

void qAnisotropic(Factorization f, const std::vector<fourier::CoefficientGrid>& g, const IndexMap& map, const QBlocks& Qbar,
                  std::vector<fourier::CoefficientGrid>& gbar) {
  const int n = map.size();
  gbar.resize(g.size());
  for (std::size_t k = 0; k < g.size(); ++k)
    if (gbar[k].Px() < g[k].Px() || gbar[k].Pz() < g[k].Pz()) gbar[k] = fourier::CoefficientGrid(g[k].Px(), g[k].Pz());
  MatC Qb(3 * n, 3 * n);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Qb.block(i * n, j * n, n, n) = block(Qbar, i, j, n);
  if (f == Factorization::Laurent) {
    blockToeplitzBar(Qb, 0, map, gbar);
    return;
  }
  if (f == Factorization::FFF_A) {
    // Q = T_eC Y, Y = T_C^-1
    const MatC TC = blockToeplitz(g, 0, map), TeC = blockToeplitz(g, 9, map);
    const MatC Y = TC.partialPivLu().inverse();
    const MatC TeCb = Qb * Y.adjoint(), Yb = TeC.adjoint() * Qb;
    MatC TCb;
    inverse(Y, Yb, TCb);
    blockToeplitzBar(TCb, 0, map, gbar);
    blockToeplitzBar(TeCb, 9, map, gbar);
    return;
  }
  FG_REQUIRE(f == Factorization::FFF_B, "adjoint::qAnisotropic: factorization Laurent, FFF_A or FFF_B");
  // Q = G1 PT + G2 Z, Z = XiB^-1 R, R = 1 - H1 PT (Plan eq:QaniB)
  const MatC G1 = blockToeplitz(g, 0, map), G2 = blockToeplitz(g, 9, map), H1 = blockToeplitz(g, 18, map),
             PT = blockToeplitz(g, 27, map);
  const MatC Xi = fourier::toeplitz(map, g[36]);
  MatC XiB = MatC::Zero(3 * n, 3 * n);
  for (int i = 0; i < 3; ++i) XiB.block(i * n, i * n, n, n) = Xi;
  const MatC R = MatC::Identity(3 * n, 3 * n) - H1 * PT;
  const MatC Z = XiB.partialPivLu().solve(R);
  const MatC G1b = Qb * PT.adjoint(), G2b = Qb * Z.adjoint(), Zb = G2.adjoint() * Qb;
  MatC PTb = G1.adjoint() * Qb;
  MatC XiBb = MatC::Zero(3 * n, 3 * n), Rb = MatC::Zero(3 * n, 3 * n);
  solve(XiB, Z, Zb, &XiBb, &Rb);
  const MatC H1b = -Rb * PT.adjoint();
  PTb -= H1.adjoint() * Rb;
  MatC Xib = MatC::Zero(n, n);
  for (int i = 0; i < 3; ++i) Xib += XiBb.block(i * n, i * n, n, n);
  blockToeplitzBar(G1b, 0, map, gbar);
  blockToeplitzBar(G2b, 9, map, gbar);
  blockToeplitzBar(H1b, 18, map, gbar);
  blockToeplitzBar(PTb, 27, map, gbar);
  toeplitz(map, Xib, gbar[36]);
}

void rkStep(const std::vector<MatC>& nodes, double h, const ButcherTableau& tab, const MatC& Pbar, std::vector<MatC>& nodesBar,
            double* hbar) {
  const int s = tab.stages();
  const Eigen::Index n = nodes[0].rows();
  MatC Abig = MatC::Identity(s * n, s * n), Mstack(s * n, n);
  for (int i = 0; i < s; ++i) {
    for (int j = 0; j < s; ++j) Abig.block(i * n, j * n, n, n) -= (h * tab.A(i, j)) * nodes[static_cast<std::size_t>(i)];
    Mstack.block(i * n, 0, n, n) = nodes[static_cast<std::size_t>(i)];
  }
  Eigen::PartialPivLU<MatC> lu(Abig);
  const MatC X = lu.solve(Mstack);
  // P = 1 + h sum_i b_i X_i
  MatC Xb(s * n, n);
  for (int i = 0; i < s; ++i) Xb.block(i * n, 0, n, n) = (h * tab.b(i)) * Pbar;
  const MatC R = MatC(Abig.adjoint()).partialPivLu().solve(Xb);  // Mstack bar
  const MatC Ab = -R * X.adjoint();
  if (hbar) {  // P = 1 + h sum b_i X_i, Abig = 1 - h (a_ij M_i)
    for (int i = 0; i < s; ++i) {
      *hbar += tab.b(i) * std::real((Pbar.adjoint() * X.block(i * n, 0, n, n)).trace());
      for (int j = 0; j < s; ++j)
        *hbar -= tab.A(i, j) * std::real((Ab.block(i * n, j * n, n, n).adjoint() * nodes[static_cast<std::size_t>(i)]).trace());
    }
  }
  if (nodesBar.size() != static_cast<std::size_t>(s)) nodesBar.assign(static_cast<std::size_t>(s), MatC());
  for (int i = 0; i < s; ++i) {
    MatC Mb = R.block(i * n, 0, n, n);
    for (int j = 0; j < s; ++j) Mb -= (h * tab.A(i, j)) * Ab.block(i * n, j * n, n, n);
    accumulate(nodesBar[static_cast<std::size_t>(i)], Mb);
  }
}

}  // namespace fg::adjoint
