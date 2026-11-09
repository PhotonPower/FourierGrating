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

}  // namespace fg::cyl
