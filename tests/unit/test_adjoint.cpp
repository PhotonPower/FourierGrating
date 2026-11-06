// Unit tests of the adjoint rules (Plan section "M19: Tests und Abnahme": every backward function against
// a central difference quotient with random matrices, h = 1e-7, real and imaginary parts, 1e-6).
// Directional form: L = Re tr(K^H out(x)) with random K, then dL/dt along a random complex direction D
// equals Re tr(xbar^H D).
#include <gtest/gtest.h>

#include <functional>
#include <random>

#include "fg/adjoint/primitives.hpp"
#include "fg/fourier/toeplitz.hpp"
#include "fg/physics/modes.hpp"

using fg::cplx;
using fg::MatC;
using fg::VecC;

namespace {

std::mt19937 rng(11);

MatC random(Eigen::Index r, Eigen::Index c, double scale = 1.0) {
  std::uniform_real_distribution<double> u(-1, 1);
  MatC M(r, c);
  for (Eigen::Index i = 0; i < r; ++i)
    for (Eigen::Index j = 0; j < c; ++j) M(i, j) = scale * cplx(u(rng), u(rng));
  return M;
}

double re(const MatC& K, const MatC& X) { return std::real((K.adjoint() * X).trace()); }

/// Central difference of L(x + t D) at t = 0.
double fd(const std::function<double(double)>& L, double h = 1e-7) { return (L(h) - L(-h)) / (2 * h); }

void expectClose(double adj, double num, double rel = 1e-6) {
  EXPECT_NEAR(adj, num, rel * std::max(1.0, std::abs(num))) << "adjoint " << adj << " fd " << num;
}

fg::SMatrix randomS(int n, double scale) {
  return {random(n, n, scale), random(n, n, scale), random(n, n, scale), random(n, n, scale)};
}
double reS(const fg::SMatrix& K, const fg::SMatrix& S) {
  return re(K.S11, S.S11) + re(K.S12, S.S12) + re(K.S21, S.S21) + re(K.S22, S.S22);
}
fg::SMatrix axpy(const fg::SMatrix& A, double t, const fg::SMatrix& D) {
  return {A.S11 + t * D.S11, A.S12 + t * D.S12, A.S21 + t * D.S21, A.S22 + t * D.S22};
}

}  // namespace

TEST(Adjoint, SolveAndInverse) {
  const int n = 6;
  const MatC A = random(n, n) + 3.0 * MatC::Identity(n, n), B = random(n, 3), K = random(n, 3), Kc = random(n, n);
  const MatC DA = random(n, n), DB = random(n, 3);
  const MatC X = A.partialPivLu().solve(B);
  MatC Ab, Bb;
  fg::adjoint::solve(A, X, K, &Ab, &Bb);
  expectClose(re(Ab, DA) + re(Bb, DB),
              fd([&](double t) { return re(K, MatC((A + t * DA).partialPivLu().solve(MatC(B + t * DB)))); }));
  const MatC C = A.inverse();
  MatC Ab2;
  fg::adjoint::inverse(C, Kc, Ab2);
  expectClose(re(Ab2, DA), fd([&](double t) { return re(Kc, MatC((A + t * DA).inverse())); }));
}

TEST(Adjoint, StarProduct) {
  const int n = 5;
  const auto A = randomS(n, 0.4), B = randomS(n, 0.4), K = randomS(n, 1.0), DA = randomS(n, 1.0), DB = randomS(n, 1.0);
  fg::SMatrix Ab, Bb;
  fg::adjoint::star(A, B, K, Ab, Bb);
  expectClose(reS(Ab, DA) + reS(Bb, DB),
              fd([&](double t) { return reS(K, fg::star(axpy(A, t, DA), axpy(B, t, DB))); }));
}

TEST(Adjoint, FromTransfer) {
  const int n = 4;
  const MatC T = random(2 * n, 2 * n) + 3.0 * MatC::Identity(2 * n, 2 * n), D = random(2 * n, 2 * n);
  const auto K = randomS(n, 1.0);
  MatC Tb;
  fg::adjoint::fromTransfer(T, K, Tb);
  expectClose(re(Tb, D), fd([&](double t) { return reS(K, fg::fromTransfer(MatC(T + t * D))); }));
}

TEST(Adjoint, PropagationMatrix) {
  const int n = 4;
  fg::QBlocks Q, D;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Q(i, j) = random(n, n, 0.3) + (i == j ? MatC(2.5 * MatC::Identity(n, n)) : MatC(MatC::Zero(n, n)));
      D(i, j) = random(n, n);
    }
  VecC al(n), ga(n);
  for (int r = 0; r < n; ++r) {
    al(r) = 0.3 + 0.4 * r;
    ga(r) = 0.1 - 0.2 * r;
  }
  const MatC K = random(4 * n, 4 * n);
  auto shifted = [&](double t) {
    fg::QBlocks q;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) q(i, j) = Q(i, j) + t * D(i, j);
    return q;
  };
  fg::QBlocks Qb;
  fg::adjoint::propagationMatrix(Q, al, ga, K, Qb);
  double adj = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) adj += re(Qb(i, j), D(i, j));
  expectClose(adj, fd([&](double t) { return re(K, fg::buildPropagationMatrix(shifted(t), al, ga)); }));
  // dense A, G (M17)
  const MatC A = random(n, n), G = random(n, n);
  fg::QBlocks Qb2;
  fg::adjoint::propagationMatrix(Q, A, G, K, Qb2);
  adj = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) adj += re(Qb2(i, j), D(i, j));
  expectClose(adj, fd([&](double t) { return re(K, fg::buildPropagationMatrix(shifted(t), A, G)); }));
}

namespace {
/// S-matrix of a slice with matrix M (eigen-decomposition + stable slice S) and its adjoint w.r.t. M.
fg::SMatrix sliceS(const MatC& M, double h, const fg::ModeBasis& ref) {
  return fg::stableLayerSMatrix(fg::layerModesAuto(M), h, ref);
}
MatC sliceBackward(const MatC& M, double h, const fg::ModeBasis& ref, const fg::SMatrix& K, double* hbar = nullptr) {
  const fg::LayerModes m = fg::layerModesAuto(M);
  MatC Wub, Wdb;
  VecC lub, ldb;
  double hb = 0;
  fg::adjoint::layerSMatrix(m, h, ref, K, Wub, Wdb, lub, ldb, hb);
  if (hbar) *hbar = hb;
  const Eigen::Index n = m.lu.size();
  MatC W(M.rows(), 2 * n), Wb(M.rows(), 2 * n);
  W << m.Wu, m.Wd;
  Wb << Wub, Wdb;
  VecC lam(2 * n), lamb(2 * n);
  lam << m.lu, m.ld;
  lamb << lub, ldb;
  std::vector<int> cls(static_cast<std::size_t>(2 * n), 0);
  for (Eigen::Index k = n; k < 2 * n; ++k) cls[static_cast<std::size_t>(k)] = 1;
  return fg::adjoint::eig(W, lam, Wb, lamb, cls);
}
}  // namespace

// A general (lossy, anisotropic-like) slice: eigen-decomposition and stable S-matrix, w.r.t. M and h.
TEST(Adjoint, SliceSMatrixGeneral) {
  const int nf = 3;
  fg::QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      Q(i, j) = random(nf, nf, 0.4) + (i == j ? MatC(cplx(3.0, 0.2) * MatC::Identity(nf, nf)) : MatC(MatC::Zero(nf, nf)));
  VecC al(nf), ga(nf);
  for (int r = 0; r < nf; ++r) {
    al(r) = 0.2 + 0.9 * (r - 1);
    ga(r) = 0.15;
  }
  const MatC M = fg::buildPropagationMatrix(Q, al, ga), D = random(4 * nf, 4 * nf, 0.2);
  const fg::ModeBasis ref = fg::isotropicModes(1.0, al, ga);
  const auto K = randomS(2 * nf, 1.0);
  const double h = 0.7;
  double hb = 0;
  const MatC Mb = sliceBackward(M, h, ref, K, &hb);
  expectClose(re(Mb, D), fd([&](double t) { return reS(K, sliceS(MatC(M + t * D), h, ref)); }), 1e-5);
  expectClose(hb, fd([&](double t) { return reS(K, sliceS(M, h + t, ref)); }), 1e-6);
}

// Exact degeneracy (homogeneous isotropic layer: TE/TM pairs, normal incidence: +-n orders): the
// gradient is finite and equals the difference quotient along a symmetric perturbation (Plan "Entartung").
TEST(Adjoint, SliceSMatrixDegenerate) {
  const int N = 2, nf = 2 * N + 1;
  VecC al(nf), ga = VecC::Zero(nf);
  for (int r = 0; r < nf; ++r) al(r) = 0.6 * (r - N);  // normal incidence: alpha_{-n} = -alpha_n
  // symmetric lamellar-like layer: Toeplitz of an even real profile
  std::vector<cplx> c(4 * N + 1);
  for (int p = -2 * N; p <= 2 * N; ++p) c[static_cast<std::size_t>(p + 2 * N)] = p == 0 ? cplx(3.0) : cplx(0.4 / (1 + p * p));
  const MatC E = fg::fourier::toeplitz1D(c, N);
  fg::QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Q(i, j) = i == j ? E : MatC(MatC::Zero(nf, nf));
  const MatC M = fg::buildPropagationMatrix(Q, al, ga);
  // symmetric perturbation (even in the order index) keeps the degeneracy structure
  std::vector<cplx> dc(4 * N + 1);
  for (int p = -2 * N; p <= 2 * N; ++p) dc[static_cast<std::size_t>(p + 2 * N)] = cplx(0.1 / (1 + std::abs(p)), 0.02);
  const MatC dE = fg::fourier::toeplitz1D(dc, N);
  fg::QBlocks dQ;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) dQ(i, j) = i == j ? dE : MatC(MatC::Zero(nf, nf));
  auto Mt = [&](double t) {
    fg::QBlocks q;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) q(i, j) = Q(i, j) + t * dQ(i, j);
    return fg::buildPropagationMatrix(q, al, ga);
  };
  const MatC D = Mt(1e-3) - Mt(-1e-3);  // direction of M (linear in t for the diagonal Q path)
  const fg::ModeBasis ref = fg::isotropicModes(1.0, al, ga);
  const auto K = randomS(2 * nf, 1.0);
  const MatC Mb = sliceBackward(M, 0.5, ref, K);
  ASSERT_TRUE(Mb.allFinite());
  expectClose(re(Mb, D) / 2e-3, fd([&](double t) { return reS(K, sliceS(Mt(t), 0.5, ref)); }), 1e-5);
}

TEST(Adjoint, QFactorizations) {
  const int N = 2;
  const fg::IndexMap map(N, 0);
  auto grid = [&](double scale, double mean) {
    fg::fourier::CoefficientGrid g(2 * N, 0);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int p = -2 * N; p <= 2 * N; ++p) g(p, 0) = scale * cplx(u(rng), u(rng)) + (p == 0 ? mean : 0.0);
    return g;
  };
  fg::NormalFieldCoefficients nn;
  for (auto& row : nn.NN)
    for (auto& g : row) g = grid(0.2, 0.3);
  for (auto& g : nn.N) g = grid(0.2, 0.5);
  nn.hasSigned = true;
  fg::SliceCoefficients c;
  c.eps = grid(0.3, 3.0);
  c.invEps = grid(0.03, 0.35);
  c.normal = &nn;
  const auto dEps = grid(1.0, 0.0), dInv = grid(0.1, 0.0);
  fg::NormalFieldCoefficients dnn;
  for (auto& row : dnn.NN)
    for (auto& g : row) g = grid(0.5, 0.0);
  for (auto& g : dnn.N) g = grid(0.5, 0.0);
  fg::QBlocks K;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) K(i, j) = random(map.size(), map.size());
  auto reQ = [&](const fg::QBlocks& Q) {
    double s = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s += re(K(i, j), Q(i, j));
    return s;
  };
  fg::NormalFieldCoefficients nnT;  // normal field at t (kept alive for the returned coefficients)
  auto shifted = [&](double t) {
    fg::SliceCoefficients d = c;
    for (int p = -2 * N; p <= 2 * N; ++p) {
      d.eps(p, 0) += t * dEps(p, 0);
      d.invEps(p, 0) += t * dInv(p, 0);
    }
    nnT = nn;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        for (int p = -2 * N; p <= 2 * N; ++p) nnT.NN[i][j](p, 0) += t * dnn.NN[i][j](p, 0);
      for (int p = -2 * N; p <= 2 * N; ++p) nnT.N[i](p, 0) += t * dnn.N[i](p, 0);
    }
    d.normal = &nnT;
    return d;
  };
  auto dot = [&](const fg::adjoint::IsoGridBar& g) {
    double s = 0;
    for (int p = -2 * N; p <= 2 * N; ++p) {
      if (g.eps.Px() >= 2 * N) s += std::real(std::conj(g.eps(p, 0)) * dEps(p, 0));
      if (g.invEps.Px() >= 2 * N) s += std::real(std::conj(g.invEps(p, 0)) * dInv(p, 0));
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          if (g.NN[i][j].Px() >= 2 * N) s += std::real(std::conj(g.NN[i][j](p, 0)) * dnn.NN[i][j](p, 0));
        if (g.N[i].Px() >= 2 * N) s += std::real(std::conj(g.N[i](p, 0)) * dnn.N[i](p, 0));
      }
    }
    return s;
  };
  {
    fg::adjoint::IsoGridBar g;
    fg::adjoint::qLaurent(c, map, K, g);
    expectClose(dot(g), fd([&](double t) { return reQ(fg::qLaurent(shifted(t), map)); }));
  }
  {
    fg::adjoint::IsoGridBar g;
    fg::adjoint::qLiLamellar(c, map, K, g);
    expectClose(dot(g), fd([&](double t) { return reQ(fg::qLiLamellar(shifted(t), map)); }));
  }
  for (bool sym : {false, true}) {
    fg::adjoint::IsoGridBar g;
    fg::adjoint::qFFFIsotropic(c, map, sym, K, g);
    expectClose(dot(g), fd([&](double t) { return reQ(fg::qFFFIsotropic(shifted(t), map, sym)); }));
  }
  {
    fg::adjoint::IsoGridBar g;
    fg::adjoint::qFFFIsotropicRigorous(c, map, K, g);
    expectClose(dot(g), fd([&](double t) { return reQ(fg::qFFFIsotropicRigorous(shifted(t), map)); }));
  }
}

TEST(Adjoint, RKStep) {
  const int n = 4;
  for (const auto& name : {"gauss2", "gauss3", "radau3"}) {
    SCOPED_TRACE(name);
    const auto tab = fg::ButcherTableau::byName(name);
    std::vector<MatC> nodes, D;
    for (int i = 0; i < tab.stages(); ++i) {
      nodes.push_back(random(n, n));
      D.push_back(random(n, n));
    }
    const MatC K = random(n, n);
    std::vector<MatC> nb;
    fg::adjoint::rkStep(nodes, 0.3, tab, K, nb);
    double adj = 0;
    for (int i = 0; i < tab.stages(); ++i) adj += re(nb[static_cast<std::size_t>(i)], D[static_cast<std::size_t>(i)]);
    expectClose(adj, fd([&](double t) {
                  std::vector<MatC> m;
                  for (int i = 0; i < tab.stages(); ++i)
                    m.push_back(nodes[static_cast<std::size_t>(i)] + t * D[static_cast<std::size_t>(i)]);
                  return re(K, fg::rkStepPropagator(m, 0.3, tab));
                }));
  }
}
