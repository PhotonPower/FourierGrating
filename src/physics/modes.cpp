#include "fg/physics/modes.hpp"

#include <cmath>

#include <algorithm>
#include <array>

#include "fg/core/error.hpp"
#include "fg/physics/propagation_matrix.hpp"

namespace fg {

cplx betaBranch(cplx eps, cplx alpha, cplx gamma, double eta) {
  cplx b = std::sqrt(eps - alpha * alpha - gamma * gamma);
  if (b.imag() < 0.0 || (b.imag() == 0.0 && b.real() < 0.0)) b = -b;
  if (std::abs(b) < eta) b += cplx(0.0, eta);
  return b;
}

ModeBasis isotropicModes(cplx eps, const VecC& alpha, const VecC& gamma, double eta) {
  const Eigen::Index n = alpha.size();
  FG_REQUIRE(gamma.size() == n, "alpha/gamma size mismatch");
  ModeBasis mb;
  mb.beta.resize(n);
  // V maps (Ex;Ez) -> (Hx;Hz) for up-going waves: V = K / beta (Plan eq:Kmat)
  MatC V = MatC::Zero(2 * n, 2 * n), Vi = MatC::Zero(2 * n, 2 * n);
  for (Eigen::Index r = 0; r < n; ++r) {
    const cplx a = alpha(r), g = gamma(r);
    const cplx b = betaBranch(eps, a, g, eta);
    mb.beta(r) = b;
    V(r, r) = a * g / b;
    V(r, n + r) = (eps - a * a) / b;
    V(n + r, r) = -(eps - g * g) / b;
    V(n + r, n + r) = -a * g / b;
    const cplx f = 1.0 / (eps * b);
    Vi(r, r) = -a * g * f;
    Vi(r, n + r) = -(eps - a * a) * f;
    Vi(n + r, r) = (eps - g * g) * f;
    Vi(n + r, n + r) = a * g * f;
  }
  const MatC Id = MatC::Identity(2 * n, 2 * n);
  mb.Psi.resize(4 * n, 4 * n);
  mb.Psi << Id, Id, -V, V;
  mb.PsiInv.resize(4 * n, 4 * n);
  mb.PsiInv << 0.5 * Id, -0.5 * Vi, 0.5 * Id, 0.5 * Vi;
  mb.lambda.resize(4 * n);
  for (Eigen::Index r = 0; r < n; ++r) {
    mb.lambda(r) = mb.lambda(n + r) = -I * mb.beta(r);              // down (a_x, a_z)
    mb.lambda(2 * n + r) = mb.lambda(3 * n + r) = I * mb.beta(r);   // up   (b_x, b_z)
  }
  return mb;
}

ModeBasis anisotropicModes(const Mat3C& eps, const VecC& alpha, const VecC& gamma) {
  // Per-order 4x4 eigenproblem (Berreman), classification Plan eq:classify (sort-and-split).
  const Eigen::Index n = alpha.size();
  ModeBasis mb;
  mb.Psi = MatC::Zero(4 * n, 4 * n);
  mb.PsiInv = MatC::Zero(4 * n, 4 * n);
  mb.lambda.resize(4 * n);
  const QBlocks Q = homogeneousTensorQ(eps, 1);
  for (Eigen::Index r = 0; r < n; ++r) {
    VecC a(1), g(1);
    a(0) = alpha(r);
    g(0) = gamma(r);
    const MatC M = buildPropagationMatrix(Q, a, g);  // 4x4: rows/cols (Ex, Ez, Hx, Hz)
    Eigen::ComplexEigenSolver<MatC> es(M, true);
    FG_REQUIRE(es.info() == Eigen::Success, "anisotropic mode eigenproblem failed");
    const VecC& lam = es.eigenvalues();
    const double tau = 1e-9 * std::max(1.0, lam.cwiseAbs().maxCoeff());
    std::array<double, 4> score{};
    for (int k = 0; k < 4; ++k) {
      const double re = lam(k).real();
      if (std::abs(re) > tau) {
        score[static_cast<std::size_t>(k)] = re;
      } else {
        const VecC w = es.eigenvectors().col(k);
        const double sy = 0.5 * std::real(w(1) * std::conj(w(2)) - w(0) * std::conj(w(3)));
        score[static_cast<std::size_t>(k)] = sy > 0 ? -0.5 * tau : 0.5 * tau;
      }
    }
    std::array<int, 4> idx{0, 1, 2, 3};
    std::sort(idx.begin(), idx.end(), [&](int i, int j) { return score[static_cast<std::size_t>(i)] < score[static_cast<std::size_t>(j)]; });
    // up: idx[0], idx[1]; down: idx[2], idx[3]. Columns: down mode 1/2 at r, n+r; up at 2n+r, 3n+r
    const int order[4] = {idx[2], idx[3], idx[0], idx[1]};
    MatC P(4, 4);
    for (int c = 0; c < 4; ++c) {
      VecC w = es.eigenvectors().col(order[c]);
      w /= w.cwiseAbs().maxCoeff();
      P.col(c) = w;
      mb.lambda(c * n + r) = lam(order[c]);
    }
    const MatC Pi = P.partialPivLu().inverse();
    for (int i = 0; i < 4; ++i)
      for (int c = 0; c < 4; ++c) {
        mb.Psi(i * n + r, c * n + r) = P(i, c);
        mb.PsiInv(c * n + r, i * n + r) = Pi(c, i);
      }
  }
  return mb;
}

double poyntingY(const VecC& F) {
  const Eigen::Index n = F.size() / 4;
  double s = 0.0;
  for (Eigen::Index r = 0; r < n; ++r) {
    const cplx Ex = F(r), Ez = F(n + r), Hx = F(2 * n + r), Hz = F(3 * n + r);
    s += 0.5 * std::real(Ez * std::conj(Hx) - Ex * std::conj(Hz));
  }
  return s;
}

}  // namespace fg
