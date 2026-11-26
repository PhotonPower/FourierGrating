#include "fg/cyl/propagation.hpp"

#include "fg/core/error.hpp"

namespace fg::cyl {

MatC propagationMatrix(const VecC& k, const CylQ& Q) {
  const Eigen::Index n = k.size();
  FG_REQUIRE(Q.Qyy.rows() == n && Q.Qpp.rows() == n, "cyl::propagationMatrix: sizes");
  const MatC Id = MatC::Identity(n, n);
  const MatC K = k.asDiagonal();
  const MatC KEK = 0.5 * (K * Q.Qyy.partialPivLu().solve(K));  // (1/2) kappa E kappa
  const MatC K2 = 0.5 * k.cwiseProduct(k).asDiagonal().toDenseMatrix();
  MatC M = MatC::Zero(4 * n, 4 * n);
  M.block(0, 2 * n, n, n) = KEK - Id;
  M.block(0, 3 * n, n, n) = KEK;
  M.block(n, 2 * n, n, n) = -KEK;
  M.block(n, 3 * n, n, n) = Id - KEK;
  M.block(2 * n, 0, n, n) = Q.Qpp - K2;
  M.block(2 * n, n, n, n) = Q.Qpm - K2;
  M.block(3 * n, 0, n, n) = K2 - Q.Qmp;
  M.block(3 * n, n, n, n) = K2 - Q.Qmm;
  return M;
}

ModeBasis homogeneousModes(cplx eps, const VecC& k, double eta) {
  const Eigen::Index n = k.size();
  ModeBasis mb;
  mb.beta.resize(n);
  mb.Psi = MatC::Zero(4 * n, 4 * n);
  mb.PsiInv = MatC::Zero(4 * n, 4 * n);
  mb.lambda.resize(4 * n);
  for (Eigen::Index j = 0; j < n; ++j) {
    const cplx b = betaBranch(eps, k(j), 0.0, eta), k2 = k(j) * k(j);
    mb.beta(j) = b;
    // calB = [[eps - k^2/2, -k^2/2], [k^2/2, -eps + k^2/2]], V = -i/beta calB, V^-1 = i/(eps beta)
    // calB
    const cplx B[2][2] = {{eps - 0.5 * k2, -0.5 * k2}, {0.5 * k2, -eps + 0.5 * k2}};
    const Eigen::Index idx[2] = {j, n + j};
    for (int r = 0; r < 2; ++r) {
      mb.Psi(idx[r], idx[r]) = 1.0;          // e from a
      mb.Psi(idx[r], 2 * n + idx[r]) = 1.0;  // e from b
      for (int c = 0; c < 2; ++c) {
        const cplx V = -I / b * B[r][c], Vi = I / (eps * b) * B[r][c];
        mb.Psi(2 * n + idx[r], idx[c]) = -V;
        mb.Psi(2 * n + idx[r], 2 * n + idx[c]) = V;
        // Psi^-1 = 1/2 [[1, -V^-1], [1, V^-1]]
        mb.PsiInv(idx[r], 2 * n + idx[c]) = -0.5 * Vi;
        mb.PsiInv(2 * n + idx[r], 2 * n + idx[c]) = 0.5 * Vi;
      }
      mb.PsiInv(idx[r], idx[r]) = 0.5;
      mb.PsiInv(2 * n + idx[r], idx[r]) = 0.5;
      mb.lambda(idx[r]) = -I * b;
      mb.lambda(2 * n + idx[r]) = I * b;
    }
  }
  return mb;
}

cplx homogeneousModesBar(cplx eps, const VecC& k, const MatC& Psibar, const MatC& PsiInvbar, double eta) {
  const Eigen::Index n = k.size();
  FG_REQUIRE((Psibar.size() == 0 || Psibar.rows() == 4 * n) && (PsiInvbar.size() == 0 || PsiInvbar.rows() == 4 * n),
             "homogeneousModesBar: sizes");
  // dX = X' d eps for the holomorphic entries X: d L = Re sum conj(Xbar) X' d eps, so epsbar = sum Xbar conj(X')
  cplx bar = 0.0;
  for (Eigen::Index j = 0; j < n; ++j) {
    const cplx b = betaBranch(eps, k(j), 0.0, eta), k2 = k(j) * k(j), db = 0.5 / b;
    const cplx B[2][2] = {{eps - 0.5 * k2, -0.5 * k2}, {0.5 * k2, -eps + 0.5 * k2}};
    const double dB[2][2] = {{1.0, 0.0}, {0.0, -1.0}};
    const Eigen::Index idx[2] = {j, n + j};
    for (int r = 0; r < 2; ++r)
      for (int c = 0; c < 2; ++c) {
        // V = -i/beta B, Vi = i/(eps beta) B
        const cplx dV = -I * (dB[r][c] / b - B[r][c] * db / (b * b));
        const cplx dVi = I * (dB[r][c] / (eps * b) - B[r][c] * (b + eps * db) / (eps * eps * b * b));
        if (Psibar.size())
          bar += (-Psibar(2 * n + idx[r], idx[c]) + Psibar(2 * n + idx[r], 2 * n + idx[c])) * std::conj(dV);
        if (PsiInvbar.size())
          bar += 0.5 * (-PsiInvbar(idx[r], 2 * n + idx[c]) + PsiInvbar(2 * n + idx[r], 2 * n + idx[c])) * std::conj(dVi);
      }
  }
  return bar;
}

double flux(const SpectralGrid& grid, const VecC& F) {
  const Eigen::Index n = grid.k.size();
  FG_REQUIRE(F.size() == 4 * n, "cyl::flux: F has 4 N_k entries");
  cplx s = 0.0;
  for (Eigen::Index j = 0; j < n; ++j) {
    if (grid.k(j).imag() != 0.0) continue;
    s += grid.w(j) * grid.k(j) *
         (F(j) * std::conj(F(2 * n + j)) - F(n + j) * std::conj(F(3 * n + j)));
  }
  return 0.5 * PI * s.imag();  // Plan eq:cyl_flux (one harmonic)
}

VecC fluxBar(const SpectralGrid& grid, const VecC& F, double scale) {
  const Eigen::Index n = grid.k.size();
  FG_REQUIRE(F.size() == 4 * n, "cyl::fluxBar: F has 4 N_k entries");
  // L = s Im(x conj(y)) per term: xbar = i s y, ybar = -i s x
  VecC Fb = VecC::Zero(4 * n);
  for (Eigen::Index j = 0; j < n; ++j) {
    if (grid.k(j).imag() != 0.0) continue;
    const double s = 0.5 * PI * scale * (grid.w(j) * grid.k(j)).real();
    Fb(j) += I * s * F(2 * n + j);
    Fb(2 * n + j) += -I * s * F(j);
    Fb(n + j) += -I * s * F(3 * n + j);
    Fb(3 * n + j) += I * s * F(n + j);
  }
  return Fb;
}

CylQ propagationMatrixBar(const VecC& k, const CylQ& Q, const MatC& Mbar) {
  const Eigen::Index n = k.size();
  FG_REQUIRE(Mbar.rows() == 4 * n && Mbar.cols() == 4 * n, "cyl::propagationMatrixBar: sizes");
  auto B = [&](int r, int c) { return Mbar.block(r * n, c * n, n, n); };
  CylQ Qb;
  Qb.Qpp = B(2, 0);
  Qb.Qpm = B(2, 1);
  Qb.Qmp = -B(3, 0);
  Qb.Qmm = -B(3, 1);
  // KEK = (1/2) K X K, X = Q_yy^{-1}: Xbar = (1/2) K^H KEKbar K^H, Q_yybar = -X^H Xbar X^H
  const MatC KEKb = B(0, 2) + B(0, 3) - B(1, 2) - B(1, 3);
  const VecC kc = k.conjugate();
  const MatC Xb = 0.5 * (kc.asDiagonal() * KEKb * kc.asDiagonal());
  const MatC X = Q.Qyy.partialPivLu().solve(MatC::Identity(n, n));
  Qb.Qyy = -(X.adjoint() * Xb * X.adjoint());
  return Qb;
}

}  // namespace fg::cyl
