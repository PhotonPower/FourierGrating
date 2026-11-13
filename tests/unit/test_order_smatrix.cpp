// S-matrix of the stack in the s/p basis of the orders (Plan eq:Ssp): the central column reproduces
// solve() for any polarization, light from below mirrors light from above for a mirror-symmetric
// stack, and Simulation gives the solver's matrices.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <string>

#include "fg/physics/efficiency.hpp"
#include "fg/physics/incidence.hpp"
#include "fg/solvers/simulation.hpp"
#include "fg/solvers/solver.hpp"

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

std::string config(double psi_deg, const char* sub, const char* lamellar_materials) {
  return std::string(R"({
  "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 6},
  "period_um": {"x": 0.7},
  "incidence": {"lambda_um": 0.6, "phi_deg": 17, "psi_deg": )") +
         std::to_string(psi_deg) + R"(, "delta_deg": 0},
  "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "hi": {"n": [2.1, 0.02]}},
  "superstrate": "air", "substrate": ")" +
         sub + R"(",
  "layers": [
    {"type": "lamellar1d", "thickness_um": 0.25, "breakpoints_um": [0.0, 0.3, 0.7], "materials": )" +
         lamellar_materials + R"(}
  ]
})";
}

/// Order index r of (n, m) in the matrix, -1 if absent.
int indexOf(const fg::OrderSMatrix& S, int n, int m) {
  for (std::size_t r = 0; r < S.orders.size(); ++r)
    if (S.orders[r].first == n && S.orders[r].second == m) return static_cast<int>(r);
  return -1;
}

void expectCentralColumnEqualsSolve(const std::string& json, double delta_deg) {
  fg::Simulation sim(json);
  fg::Wavelength w = sim.config().scan.front();
  w.delta = delta_deg * kDeg;
  const fg::Result res = sim.solve(w);
  const fg::OrderSMatrix S = sim.orderSMatrix(w);
  const int nf = static_cast<int>(S.orders.size());
  const int c = indexOf(S, 0, 0);
  ASSERT_GE(c, 0);
  // the incident wave as the solver makes it, in the s/p basis of the incident order
  const fg::IncidentWave inc = fg::makeIncidence(1.0, w.phi, w.psi, w.delta);
  const auto o = fg::planeWaveOrder(1.0, inc.alpha0, inc.gamma0, inc.beta0, -1, inc.A(0), inc.A(2));
  fg::VecC in = fg::VecC::Zero(2 * nf);
  in(c) = o.amp_s;
  in(nf + c) = o.amp_p;
  const fg::VecC r = S.R_top * in, t = S.T_down * in;
  double scale = 0;
  for (const auto& e : res.reflected)
    scale = std::max(scale, std::abs(e.amp_s) + std::abs(e.amp_p));
  ASSERT_GT(scale, 0);
  for (const auto& e : res.reflected) {
    const int k = indexOf(S, e.n, e.m);
    ASSERT_GE(k, 0);
    EXPECT_NEAR(std::abs(r(k) - e.amp_s), 0.0, 1e-12 * scale) << "R s order " << e.n;
    EXPECT_NEAR(std::abs(r(nf + k) - e.amp_p), 0.0, 1e-12 * scale) << "R p order " << e.n;
  }
  for (const auto& e : res.transmitted) {
    const int k = indexOf(S, e.n, e.m);
    ASSERT_GE(k, 0);
    EXPECT_NEAR(std::abs(t(k) - e.amp_s), 0.0, 1e-12 * scale) << "T s order " << e.n;
    EXPECT_NEAR(std::abs(t(nf + k) - e.amp_p), 0.0, 1e-12 * scale) << "T p order " << e.n;
  }
}

}  // namespace

TEST(OrderSMatrix, CentralColumnEqualsSolveConical) {
  const std::string json = config(10, "glass", R"(["hi", "air"])");
  for (double delta : {0.0, 35.0, 90.0}) expectCentralColumnEqualsSolve(json, delta);
}

TEST(OrderSMatrix, CentralColumnEqualsSolveInPlane) {
  // 1D in-plane: solve() computes only the excited polarization; the matrix holds both
  const std::string json = config(0, "glass", R"(["hi", "air"])");
  for (double delta : {0.0, 90.0, 60.0}) expectCentralColumnEqualsSolve(json, delta);
}

TEST(OrderSMatrix, LightFromBelowMirrorsLightFromAbove) {
  // a lamellar layer between equal media is symmetric under y -> -y; the mirror keeps s and turns
  // p into -p (p = k x s, M(a x b) = -(Ma x Mb)), so R_bottom = D R_top D and T_up = D T_down D
  fg::Simulation sim(config(10, "air", R"(["hi", "air"])"));
  const fg::OrderSMatrix S = sim.orderSMatrix(sim.config().scan.front());
  const int nf = static_cast<int>(S.orders.size());
  fg::VecC d = fg::VecC::Ones(2 * nf);
  d.tail(nf).setConstant(-1.0);
  const fg::MatC D = d.asDiagonal();
  const double scale = S.R_top.cwiseAbs().maxCoeff() + S.T_down.cwiseAbs().maxCoeff();
  EXPECT_LT((S.R_bottom - D * S.R_top * D).cwiseAbs().maxCoeff(), 1e-11 * scale);
  EXPECT_LT((S.T_up - D * S.T_down * D).cwiseAbs().maxCoeff(), 1e-11 * scale);
  EXPECT_GT(S.R_top.cwiseAbs().maxCoeff(), 1e-3);  // content guard: it does reflect
}

TEST(OrderSMatrix, SimulationEqualsSolver) {
  const std::string json = config(10, "glass", R"(["hi", "air"])");
  fg::Simulation sim(json);
  const auto& c = sim.config();
  const fg::OrderSMatrix a = sim.orderSMatrix(c.scan.front());
  const fg::OrderSMatrix b = fg::makeSolver(c.solver)->orderSMatrix(c.stack, c.scan.front());
  EXPECT_EQ(a.orders, b.orders);
  EXPECT_TRUE(a.R_top == b.R_top && a.T_down == b.T_down && a.T_up == b.T_up &&
              a.R_bottom == b.R_bottom);
}
