#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

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

// ---- M15: star-product chains (Plan section "M15", step 3) ----
namespace {
bool bitEqual(const fg::SMatrix& a, const fg::SMatrix& b) {
  const MatC fa = a.full(), fb = b.full();
  return fa.size() == fb.size() && std::equal(fa.data(), fa.data() + fa.size(), fb.data());
}
std::vector<fg::SMatrix> randomChain(int K, int nh, unsigned seed) {
  std::mt19937 rng(seed);
  std::vector<fg::SMatrix> v;
  for (int k = 0; k < K; ++k) v.push_back(fg::fromTransfer(randomTransfer(2 * nh, rng)));
  return v;
}
}  // namespace

TEST(StarChain, SequentialIsTheFormerLoop) {
  const auto v = randomChain(9, 4, 11);
  fg::SMatrix ref = v[0];
  for (std::size_t k = 1; k < v.size(); ++k) ref = fg::star(ref, v[k]);
  for (int threads : {1, 3, 8}) {
    fg::StarChain c(false, threads);
    for (const auto& S : v) c.push(S);
    EXPECT_TRUE(bitEqual(c.result(), ref)) << threads;
  }
}

TEST(StarChain, TreeHasFixedShapeForEveryThreadCount) {
  // bit-identical to the level-wise reference for every chain length and block size
  for (int K : {1, 2, 3, 5, 6, 7, 8, 11, 13, 16, 17, 31, 40}) {
    const auto v = randomChain(K, 3, 100u + static_cast<unsigned>(K));
    const fg::SMatrix ref = fg::starTree(v, 1);
    for (int threads : {1, 2, 3, 4, 5, 8, 16}) {
      fg::StarChain c(true, threads);
      for (const auto& S : v) c.push(S);
      EXPECT_TRUE(bitEqual(c.result(), ref)) << "K = " << K << ", threads = " << threads;
      EXPECT_TRUE(bitEqual(fg::starTree(v, threads), ref)) << "K = " << K << ", threads = " << threads;
    }
  }
}

TEST(StarChain, TreeEqualsSequentialTo1e12) {
  const auto v = randomChain(40, 6, 5);
  fg::StarChain seq(false), tree(true, 4);
  for (const auto& S : v) {
    seq.push(S);
    tree.push(S);
  }
  const MatC a = seq.result().full(), b = tree.result().full();
  EXPECT_LT((a - b).cwiseAbs().maxCoeff(), 1e-12 * a.cwiseAbs().maxCoeff());
}
