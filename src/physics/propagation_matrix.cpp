#include "fg/physics/propagation_matrix.hpp"

#include "fg/core/error.hpp"
#include "fg/core/prof.hpp"

namespace fg {

QBlocks homogeneousIsotropicQ(cplx eps, int nf) {
  QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      Q(i, j) = (i == j) ? MatC(eps * MatC::Identity(nf, nf)) : MatC(MatC::Zero(nf, nf));
  return Q;
}

QBlocks homogeneousTensorQ(const Mat3C& eps, int nf) {
  QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Q(i, j) = eps(i, j) * MatC::Identity(nf, nf);
  return Q;
}

MatC buildPropagationMatrix(const QBlocks& Q, const VecC& alpha, const VecC& gamma) {
  FG_PROF_SCOPE("physics.M");
  const int n = Q.nf();
  FG_REQUIRE(alpha.size() == n && gamma.size() == n, "alpha/gamma size mismatch");
  constexpr int X = 0, Y = 1, Z = 2;
  const auto al = alpha.asDiagonal();
  const auto ga = gamma.asDiagonal();
  const MatC Id = MatC::Identity(n, n);
  Eigen::PartialPivLU<MatC> Eyy(Q(Y, Y));  // E = Q_yy^-1, applied via solves
  const MatC E_yx = Eyy.solve(Q(Y, X));
  const MatC E_yz = Eyy.solve(Q(Y, Z));
  const MatC E_g = Eyy.solve(MatC(ga));
  const MatC E_a = Eyy.solve(MatC(al));

  MatC M = MatC::Zero(4 * n, 4 * n);
  auto B = [&](int i, int j) { return M.block(i * n, j * n, n, n); };
  // rows 1, 2 (Plan eq:Mmu1a)
  B(0, 0) = -I * (al * E_yx);
  B(0, 1) = -I * (al * E_yz);
  B(0, 2) = -I * (al * E_g);
  B(0, 3) = I * (al * E_a) - I * Id;
  B(1, 0) = -I * (ga * E_yx);
  B(1, 1) = -I * (ga * E_yz);
  B(1, 2) = -I * (ga * E_g) + I * Id;
  B(1, 3) = I * (ga * E_a);
  // rows 3, 4 (Plan eq:Mmu1b)
  const MatC alga = (alpha.array() * gamma.array()).matrix().asDiagonal();
  const MatC al2 = (alpha.array() * alpha.array()).matrix().asDiagonal();
  const MatC ga2 = (gamma.array() * gamma.array()).matrix().asDiagonal();
  B(2, 0) = I * (Q(Z, X) - Q(Z, Y) * E_yx + alga);
  B(2, 1) = I * (Q(Z, Z) - Q(Z, Y) * E_yz - al2);
  B(2, 2) = -I * (Q(Z, Y) * E_g);
  B(2, 3) = I * (Q(Z, Y) * E_a);
  B(3, 0) = I * (Q(X, Y) * E_yx - Q(X, X) + ga2);
  B(3, 1) = I * (Q(X, Y) * E_yz - Q(X, Z) - alga);
  B(3, 2) = I * (Q(X, Y) * E_g);
  B(3, 3) = -I * (Q(X, Y) * E_a);
  return M;
}

MatC buildPropagationMatrix(const QBlocks& Q, const VecC& alpha, Polarization p) {
  FG_PROF_SCOPE("physics.M");
  const int n = Q.nf();
  FG_REQUIRE(alpha.size() == n, "alpha size mismatch");
  constexpr int X = 0, Y = 1, Z = 2;
  const auto al = alpha.asDiagonal();
  const MatC Id = MatC::Identity(n, n);
  MatC M = MatC::Zero(2 * n, 2 * n);
  if (p == Polarization::TE) {  // Plan eq:TE
    const MatC al2 = (alpha.array() * alpha.array()).matrix().asDiagonal();
    M.topRightCorner(n, n) = I * Id;
    M.bottomLeftCorner(n, n) = I * (Q(Z, Z) - al2);
    return M;
  }
  // Plan eq:TM
  Eigen::PartialPivLU<MatC> Eyy(Q(Y, Y));
  const MatC E_yx = Eyy.solve(Q(Y, X));
  const MatC E_a = Eyy.solve(MatC(al));
  M.topLeftCorner(n, n) = -I * (al * E_yx);
  M.topRightCorner(n, n) = I * (al * E_a) - I * Id;
  M.bottomLeftCorner(n, n) = I * (Q(X, Y) * E_yx - Q(X, X));
  M.bottomRightCorner(n, n) = -I * (Q(X, Y) * E_a);
  return M;
}

VecC computeEy(const QBlocks& Q, const VecC& alpha, const VecC& gamma, const VecC& F) {
  const int n = Q.nf();
  FG_REQUIRE(F.size() == 4 * n, "state vector size mismatch");
  const VecC Ex = F.segment(0, n), Ez = F.segment(n, n), Hx = F.segment(2 * n, n),
             Hz = F.segment(3 * n, n);
  const VecC rhs = alpha.cwiseProduct(Hz) - gamma.cwiseProduct(Hx) - Q(1, 0) * Ex - Q(1, 2) * Ez;
  return Q(1, 1).partialPivLu().solve(rhs);
}

VecC computeHy(const VecC& alpha, const VecC& gamma, const VecC& F) {
  const Eigen::Index n = alpha.size();
  return gamma.cwiseProduct(F.segment(0, n)) - alpha.cwiseProduct(F.segment(n, n));
}

}  // namespace fg
