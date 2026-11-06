#include "fg/solvers/rk.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/core/prof.hpp"

namespace fg {

namespace {
double stepFactor(const ButcherTableau& tab) { return tab.stages() >= 3 ? 1.0 : 0.5; }

MatC matPow(MatC P, long k) {
  MatC R = MatC::Identity(P.rows(), P.cols());
  while (k > 0) {
    if (k & 1) R = P * R;
    P = P * P;
    k >>= 1;
  }
  return R;
}
}  // namespace

MatC rkTransferConstant(const MatC& M, double H, const ButcherTableau& tab, double tol,
                        double rho, RKStats& stats, RKSteps* record) {
  FG_PROF_SCOPE("rk.transfer");
  long k = std::max(1L, static_cast<long>(std::ceil(H * rho / stepFactor(tab))));
  const std::size_t s = static_cast<std::size_t>(tab.stages());
  MatC T1 = matPow(rkStepPropagator(std::vector<MatC>(s, M), H / static_cast<double>(k), tab), k);
  for (int iter = 0; iter < 30; ++iter) {
    const MatC T2 =
        matPow(rkStepPropagator(std::vector<MatC>(s, M), H / static_cast<double>(2 * k), tab), 2 * k);
    stats.steps += 3 * k;
    const double err = (T1 - T2).norm() / T2.norm();
    if (err <= tol) {
      if (record)
        for (long i = 0; i < 2 * k; ++i) {
          record->y.push_back(H * static_cast<double>(i) / static_cast<double>(2 * k));
          record->h.push_back(H / static_cast<double>(2 * k));
        }
      return T2;
    }
    ++stats.rejected;
    T1 = T2;
    k *= 2;
  }
  throw Error("rkTransferConstant: no convergence");
}

MatC rkTransfer(const std::function<MatC(double)>& M, double y0, double y1,
                const ButcherTableau& tab, double tol, double rho, RKStats& stats, RKSteps* record) {
  FG_PROF_SCOPE("rk.transfer");
  const Eigen::Index n = M(y0).rows();
  MatC T = MatC::Identity(n, n);
  const int s = tab.stages();
  auto step = [&](double y, double h) {
    std::vector<MatC> nodes;
    for (int i = 0; i < s; ++i) nodes.push_back(M(y + tab.c(i) * h));
    return rkStepPropagator(nodes, h, tab);
  };
  double y = y0;
  double h = std::min(y1 - y0, stepFactor(tab) / std::max(rho, 1e-12));
  const double hmin = 1e-12 * (y1 - y0);
  while (y < y1 - 1e-15 * (y1 - y0)) {
    h = std::min(h, y1 - y);
    const MatC Pfull = step(y, h);
    const MatC Phalf = step(y + 0.5 * h, 0.5 * h) * step(y, 0.5 * h);
    stats.steps += 3;
    const double err = (Pfull - Phalf).norm() / Phalf.norm();
    if (err <= tol || h <= hmin) {
      T = Phalf * T;
      if (record) {
        record->y.insert(record->y.end(), {y, y + 0.5 * h});
        record->h.insert(record->h.end(), {0.5 * h, 0.5 * h});
      }
      y += h;
      const double fac = err > 0 ? 0.9 * std::pow(tol / err, 1.0 / (tab.order + 1)) : 2.0;
      h *= std::clamp(fac, 0.2, 2.0);
    } else {
      ++stats.rejected;
      h *= std::clamp(0.9 * std::pow(tol / err, 1.0 / (tab.order + 1)), 0.1, 0.5);
    }
  }
  return T;
}

MatC rkTransferAuto(const std::function<MatC(double)>& M, double y0, double y1,
                    const ButcherTableau& tab, double tol, double rho, RKStats& stats) {
  const MatC Mmid = M(0.5 * (y0 + y1));
  const Eigen::Index n = Mmid.rows() / 4;
  // index sets: TM-like (Ex, Hz) = blocks 0, 3; TE-like (Ez, Hx) = blocks 1, 2
  const int A[2] = {0, 3}, B[2] = {1, 2};
  double cross = 0;
  for (int i : A)
    for (int j : B) {
      cross = std::max(cross, Mmid.block(i * n, j * n, n, n).cwiseAbs().maxCoeff());
      cross = std::max(cross, Mmid.block(j * n, i * n, n, n).cwiseAbs().maxCoeff());
    }
  if (n == 0 || cross > 0.0) return rkTransfer(M, y0, y1, tab, tol, rho, stats);
  auto extract = [n](const MatC& Mf, const int* set) {
    MatC R(2 * n, 2 * n);
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) R.block(a * n, b * n, n, n) = Mf.block(set[a] * n, set[b] * n, n, n);
    return R;
  };
  MatC T = MatC::Zero(4 * n, 4 * n);
  for (const int* set : {A, B}) {
    const MatC Tb = rkTransfer([&](double y) { return extract(M(y), set); }, y0, y1, tab, tol, rho, stats);
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) T.block(set[a] * n, set[b] * n, n, n) = Tb.block(a * n, b * n, n, n);
  }
  return T;
}

SMatrix rkSliceSMatrix(const MatC& TF, const ModeBasis& ref) {
  FG_PROF_SCOPE("rk.slice_smatrix");
  return fromTransfer(linalg::product(linalg::product(ref.PsiInv, TF), ref.Psi));
}

}  // namespace fg
