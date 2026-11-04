#include "fg/solvers/rk.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"

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
                        double rho, RKStats& stats) {
  long k = std::max(1L, static_cast<long>(std::ceil(H * rho / stepFactor(tab))));
  const std::size_t s = static_cast<std::size_t>(tab.stages());
  MatC T1 = matPow(rkStepPropagator(std::vector<MatC>(s, M), H / static_cast<double>(k), tab), k);
  for (int iter = 0; iter < 30; ++iter) {
    const MatC T2 =
        matPow(rkStepPropagator(std::vector<MatC>(s, M), H / static_cast<double>(2 * k), tab), 2 * k);
    stats.steps += 3 * k;
    const double err = (T1 - T2).norm() / T2.norm();
    if (err <= tol) return T2;
    ++stats.rejected;
    T1 = T2;
    k *= 2;
  }
  throw Error("rkTransferConstant: no convergence");
}

MatC rkTransfer(const std::function<MatC(double)>& M, double y0, double y1,
                const ButcherTableau& tab, double tol, double rho, RKStats& stats) {
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

SMatrix rkSliceSMatrix(const MatC& TF, const ModeBasis& ref) {
  return fromTransfer(ref.PsiInv * TF * ref.Psi);
}

}  // namespace fg
