#include "fg/solvers/butcher.hpp"

#include <cmath>

#include "fg/core/error.hpp"

namespace fg {

ButcherTableau ButcherTableau::gauss2() {
  const double s3 = std::sqrt(3.0);
  ButcherTableau t;
  t.name = "gauss2";
  t.order = 4;
  t.A.resize(2, 2);
  t.A << 0.25, 0.25 - s3 / 6.0, 0.25 + s3 / 6.0, 0.25;
  t.b.resize(2);
  t.b << 0.5, 0.5;
  t.c.resize(2);
  t.c << 0.5 - s3 / 6.0, 0.5 + s3 / 6.0;
  return t;
}

ButcherTableau ButcherTableau::gauss3() {
  const double s15 = std::sqrt(15.0);
  ButcherTableau t;
  t.name = "gauss3";
  t.order = 6;
  t.A.resize(3, 3);
  t.A << 5.0 / 36, 2.0 / 9 - s15 / 15, 5.0 / 36 - s15 / 30,  //
      5.0 / 36 + s15 / 24, 2.0 / 9, 5.0 / 36 - s15 / 24,      //
      5.0 / 36 + s15 / 30, 2.0 / 9 + s15 / 15, 5.0 / 36;
  t.b.resize(3);
  t.b << 5.0 / 18, 4.0 / 9, 5.0 / 18;
  t.c.resize(3);
  t.c << 0.5 - s15 / 10, 0.5, 0.5 + s15 / 10;
  return t;
}

ButcherTableau ButcherTableau::radau3() {
  const double s6 = std::sqrt(6.0);
  ButcherTableau t;
  t.name = "radau3";
  t.order = 5;
  t.A.resize(3, 3);
  t.A << (88 - 7 * s6) / 360, (296 - 169 * s6) / 1800, (-2 + 3 * s6) / 225,  //
      (296 + 169 * s6) / 1800, (88 + 7 * s6) / 360, (-2 - 3 * s6) / 225,     //
      (16 - s6) / 36, (16 + s6) / 36, 1.0 / 9;
  t.b = t.A.row(2).transpose();
  t.c.resize(3);
  t.c << (4 - s6) / 10, (4 + s6) / 10, 1.0;
  return t;
}

ButcherTableau ButcherTableau::byName(const std::string& name) {
  if (name == "gauss2") return gauss2();
  if (name == "gauss3") return gauss3();
  if (name == "radau3") return radau3();
  throw Error("unknown Butcher tableau: " + name);
}

MatC rkStepPropagator(const std::vector<MatC>& M_nodes, double h, const ButcherTableau& tab) {
  const int s = tab.stages();
  FG_REQUIRE(static_cast<int>(M_nodes.size()) == s, "need M at every stage node");
  const Eigen::Index n = M_nodes[0].rows();
  // (1 - h Abig) X = Mstack,  Abig_ij = a_ij M_i   (Plan eq:rkprop)
  MatC Abig = MatC::Identity(s * n, s * n);
  MatC Mstack(s * n, n);
  for (int i = 0; i < s; ++i) {
    for (int j = 0; j < s; ++j)
      Abig.block(i * n, j * n, n, n) -= (h * tab.A(i, j)) * M_nodes[static_cast<std::size_t>(i)];
    Mstack.block(i * n, 0, n, n) = M_nodes[static_cast<std::size_t>(i)];
  }
  const MatC X = Abig.partialPivLu().solve(Mstack);
  MatC P = MatC::Identity(n, n);
  for (int i = 0; i < s; ++i) P += (h * tab.b(i)) * X.block(i * n, 0, n, n);
  return P;
}

}  // namespace fg
