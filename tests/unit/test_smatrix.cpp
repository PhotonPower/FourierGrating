#include <gtest/gtest.h>

#include <random>

#include "fg/smatrix/smatrix.hpp"

using fg::MatC;

namespace {
MatC randomTransfer(int n, std::mt19937& rng) {
  std::normal_distribution<double> nd;
  MatC T = MatC::Identity(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) T(i, j) += 0.3 * fg::cplx(nd(rng), nd(rng));
  return T;
}
}  // namespace

TEST(SMatrix, StarEqualsTransferProduct) {
  std::mt19937 rng(42);
  const int nh = 4;
  MatC Tacc = randomTransfer(2 * nh, rng);
  fg::SMatrix Sacc = fg::fromTransfer(Tacc);
  for (int k = 0; k < 5; ++k) {
    const MatC T = randomTransfer(2 * nh, rng);
    Sacc = fg::star(Sacc, fg::fromTransfer(T));
    Tacc = T * Tacc;
  }
  EXPECT_LT((Sacc.full() - fg::fromTransfer(Tacc).full()).cwiseAbs().maxCoeff(), 1e-10);
}

TEST(SMatrix, IdentityIsNeutral) {
  std::mt19937 rng(7);
  const fg::SMatrix S = fg::fromTransfer(randomTransfer(6, rng));
  const fg::SMatrix E = fg::SMatrix::identity(3);
  EXPECT_LT((fg::star(E, S).full() - S.full()).cwiseAbs().maxCoeff(), 1e-13);
  EXPECT_LT((fg::star(S, E).full() - S.full()).cwiseAbs().maxCoeff(), 1e-13);
}

TEST(SMatrix, InterfaceAmplitudesMatchPropagation) {
  std::mt19937 rng(3);
  const int nh = 3;
  std::vector<MatC> T;
  for (int k = 0; k < 4; ++k) T.push_back(randomTransfer(2 * nh, rng));
  const fg::SMatrix L = fg::star(fg::fromTransfer(T[0]), fg::fromTransfer(T[1]));
  const fg::SMatrix U = fg::star(fg::fromTransfer(T[2]), fg::fromTransfer(T[3]));
  const fg::VecC b_in = fg::VecC::Random(nh), a_in = fg::VecC::Random(nh);
  const fg::SMatrix tot = fg::star(L, U);
  fg::VecC inp(2 * nh);
  inp << b_in, a_in;
  const fg::VecC a_out = tot.S21 * b_in + tot.S22 * a_in;
  fg::VecC v(2 * nh);
  v << a_out, b_in;  // (a; b) at the bottom
  v = T[1] * (T[0] * v);
  EXPECT_LT((fg::interfaceAmplitudes(L, U, b_in, a_in) - v).cwiseAbs().maxCoeff(), 1e-11);
}
