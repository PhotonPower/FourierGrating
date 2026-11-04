#include "fg/physics/qmatrix.hpp"

#include "fg/core/error.hpp"
#include "fg/fourier/toeplitz.hpp"

namespace fg {

namespace {
QBlocks zeros(int nf) {
  QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Q(i, j) = MatC::Zero(nf, nf);
  return Q;
}
}  // namespace

QBlocks qLaurent(const SliceCoefficients& c, const IndexMap& map) {
  QBlocks Q = zeros(map.size());
  if (c.isotropic) {
    const MatC E = fourier::toeplitz(map, c.eps);
    for (int i = 0; i < 3; ++i) Q(i, i) = E;
  } else {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        Q(i, j) = fourier::toeplitz(map, c.epsTensor[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
  }
  return Q;
}

QBlocks qLiLamellar(const SliceCoefficients& c, const IndexMap& map) {
  FG_REQUIRE(c.isotropic, "qLiLamellar requires an isotropic layer");
  QBlocks Q = zeros(map.size());
  const MatC E = fourier::toeplitz(map, c.eps);
  Q(0, 0) = fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  Q(1, 1) = E;
  Q(2, 2) = E;
  return Q;
}

QBlocks qFFFIsotropic(const SliceCoefficients& c, const IndexMap& map, bool symmetric) {
  FG_REQUIRE(c.isotropic, "qFFFIsotropic requires an isotropic layer");
  FG_REQUIRE(c.normal != nullptr, "qFFFIsotropic requires a normal vector field");
  const MatC E = fourier::toeplitz(map, c.eps);
  const MatC Delta = E - fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  QBlocks Q;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const MatC NN = fourier::toeplitz(
          map, c.normal->NN[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
      Q(i, j) = symmetric ? MatC(-0.5 * (Delta * NN + NN * Delta)) : MatC(-Delta * NN);
      if (i == j) Q(i, j) += E;
    }
  return Q;
}

QBlocks qFFFIsotropicRigorous(const SliceCoefficients& c, const IndexMap& map) {
  FG_REQUIRE(c.isotropic && c.normal && c.normal->hasSigned, "qFFFIsotropicRigorous needs signed N");
  const MatC E = fourier::toeplitz(map, c.eps);
  const MatC L = fourier::toeplitz(map, c.invEps).partialPivLu().inverse();
  std::array<MatC, 3> Nt;
  for (std::size_t i = 0; i < 3; ++i) Nt[i] = fourier::toeplitz(map, c.normal->N[i]);
  QBlocks Q;
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) {
      Q(int(i), int(j)) = -E * (Nt[i] * Nt[j]) + Nt[i] * (L * Nt[j]);
      if (i == j) Q(int(i), int(j)) += E;
    }
  return Q;
}

namespace {
MatC blockToeplitz(const std::array<std::array<fourier::CoefficientGrid, 3>, 3>& g, const IndexMap& map) {
  const int n = map.size();
  MatC T(3 * n, 3 * n);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) T.block(int(i) * n, int(j) * n, n, n) = fourier::toeplitz(map, g[i][j]);
  return T;
}
QBlocks split(const MatC& Q, int n) {
  QBlocks out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) out(i, j) = Q.block(i * n, j * n, n, n);
  return out;
}
}  // namespace

QBlocks qFFFAnisotropicA(const SliceCoefficients& c, const IndexMap& map) {
  // Q = [[eps C]] [[C]]^-1  ->  Q^T = [[C]]^-T [[eps C]]^T
  const MatC TC = blockToeplitz(c.C, map), TeC = blockToeplitz(c.epsC, map);
  const MatC Q = TC.transpose().partialPivLu().solve(TeC.transpose()).transpose();
  return split(Q, map.size());
}

QBlocks qFFFAnisotropicB(const SliceCoefficients& c, const IndexMap& map) {
  // Q = [[G1]][[P_T]] + [[G2]] [[1/xi]]^-1 (1 - [[H1]][[P_T]])  (Plan eq:QaniB)
  const int n = map.size();
  const MatC G1 = blockToeplitz(c.G1, map), G2 = blockToeplitz(c.G2, map);
  const MatC H1 = blockToeplitz(c.H1, map), PT = blockToeplitz(c.PT, map);
  const MatC Xi = fourier::toeplitz(map, c.invXi);
  MatC XiB = MatC::Zero(3 * n, 3 * n);
  for (int i = 0; i < 3; ++i) XiB.block(i * n, i * n, n, n) = Xi;
  const MatC rhs = MatC::Identity(3 * n, 3 * n) - H1 * PT;
  const MatC Q = G1 * PT + G2 * XiB.partialPivLu().solve(rhs);
  return split(Q, n);
}

}  // namespace fg
