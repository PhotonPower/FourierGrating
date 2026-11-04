// FourierGrating - propagation matrix M of dF/dy = M F, F = ([Ex],[Ez],[Hx],[Hz])
// (Plan eq:Mmu1a/eq:Mmu1b, normalized units, mu = 1).
#pragma once

#include <array>

#include "fg/core/types.hpp"

namespace fg {

/// 3x3 blocks Q_ij (i,j in x=0, y=1, z=2) with [D] = Q [E]; each block n_f x n_f.
struct QBlocks {
  std::array<std::array<MatC, 3>, 3> q;
  const MatC& operator()(int i, int j) const { return q[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]; }
  MatC& operator()(int i, int j) { return q[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]; }
  int nf() const { return static_cast<int>(q[0][0].rows()); }
};

/// Q blocks of a homogeneous isotropic medium: Q_ij = eps * delta_ij * 1.
QBlocks homogeneousIsotropicQ(cplx eps, int nf);
/// Q blocks of a homogeneous anisotropic medium: Q_ij = eps_ij * 1.
QBlocks homogeneousTensorQ(const Mat3C& eps, int nf);

/// Full 4 n_f x 4 n_f matrix M for Q_mu = 1 (Plan eq:Mmu1a, eq:Mmu1b).
/// alpha, gamma: normalized wave numbers of the retained orders (single-index order).
MatC buildPropagationMatrix(const QBlocks& Qe, const VecC& alpha, const VecC& gamma);

/// [E_y] from F (Plan eq:Ey) and [H_y] from F (Plan eq:Hy, Q_mu = 1).
VecC computeEy(const QBlocks& Qe, const VecC& alpha, const VecC& gamma, const VecC& F);
VecC computeHy(const VecC& alpha, const VecC& gamma, const VecC& F);

}  // namespace fg
