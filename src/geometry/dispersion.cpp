// FourierGrating - dispersive material models (Plan section "Dispersive Materialien (M22)").
#include "fg/geometry/dispersion.hpp"

#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "fg/core/error.hpp"

namespace fg {

namespace {


std::string fmt(double v) {
  std::ostringstream s;
  s.precision(6);
  s << v;
  return s.str();
}

const char* modelKey(const std::variant<Dispersion::Sellmeier, Dispersion::Cauchy, Dispersion::DrudeLorentz,
                                        Dispersion::Table>& m) {
  switch (m.index()) {
    case 0: return "sellmeier";
    case 1: return "cauchy";
    case 2: return "drude_lorentz";
    default: return "table";
  }
}

std::size_t index(const std::string& tok, std::size_t size, const std::string& what) {
  std::size_t pos = 0;
  unsigned long v = 0;
  try {
    v = std::stoul(tok, &pos);
  } catch (const std::exception&) {
    pos = 0;
  }
  if (pos != tok.size() || v >= size) throw Error(what + ": index '" + tok + "' out of range");
  return static_cast<std::size_t>(v);
}

/// Interval j with x_j <= x <= x_{j+1} (x inside the nodes; the last node uses the last interval).
std::size_t interval(const std::vector<double>& x, double v) {
  const auto it = std::upper_bound(x.begin(), x.end(), v);
  const std::size_t j = static_cast<std::size_t>(std::max<std::ptrdiff_t>(it - x.begin(), 1)) - 1;
  return std::min(j, x.size() - 2);
}

/// Value and derivative of one table column at lambda inside the nodes (Plan "Tabelle").
std::pair<double, double> column(const Dispersion::Table& t, const std::vector<double>& y,
                                 const std::vector<double>& M, double lam) {
  const auto& x = t.lambda_um;
  const std::size_t j = interval(x, lam);
  const double h = x[j + 1] - x[j], s = lam - x[j];
  if (!t.cubic) {
    const double slope = (y[j + 1] - y[j]) / h;
    return {lam == x[j] ? y[j] : y[j] + slope * s, slope};
  }
  // s(x) = y_j + b t + M_j t^2 / 2 + (M_{j+1} - M_j) t^3 / (6 h); constant data stay exact
  const double b = (y[j + 1] - y[j]) / h - h * (2 * M[j] + M[j + 1]) / 6;
  const double c3 = (M[j + 1] - M[j]) / (6 * h);
  const double v = lam == x[j] ? y[j] : y[j] + s * (b + s * (0.5 * M[j] + s * c3));
  return {v, b + s * (M[j] + 3 * c3 * s)};
}

bool outside(const Dispersion::Table& t, double lam) {
  return lam < t.lambda_um.front() || lam > t.lambda_um.back();
}

}  // namespace

std::vector<double> splineSecondDerivatives(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t K = x.size();
  FG_REQUIRE(K >= 4 && y.size() == K, "cubic spline needs at least 4 nodes");
  // interior rows: h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1} = 6 (d_i - d_{i-1});
  // not-a-knot at x_1 and x_{K-2}: the third derivative is continuous there
  using Sp = Eigen::SparseMatrix<double>;
  std::vector<Eigen::Triplet<double>> tr;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(K));
  auto h = [&](std::size_t i) { return x[i + 1] - x[i]; };
  const auto e = [](std::size_t i) { return static_cast<int>(i); };
  tr.emplace_back(0, 0, h(1));
  tr.emplace_back(0, 1, -(h(0) + h(1)));
  tr.emplace_back(0, 2, h(0));
  for (std::size_t i = 1; i + 1 < K; ++i) {
    tr.emplace_back(e(i), e(i - 1), h(i - 1));
    tr.emplace_back(e(i), e(i), 2 * (h(i - 1) + h(i)));
    tr.emplace_back(e(i), e(i + 1), h(i));
    rhs(e(i)) = 6 * ((y[i + 1] - y[i]) / h(i) - (y[i] - y[i - 1]) / h(i - 1));
  }
  // the last row replaces nothing: rows 1..K-2 are the interior equations, row K-1 the second condition
  tr.emplace_back(e(K - 1), e(K - 3), h(K - 2));
  tr.emplace_back(e(K - 1), e(K - 2), -(h(K - 3) + h(K - 2)));
  tr.emplace_back(e(K - 1), e(K - 1), h(K - 3));
  Sp A(static_cast<Eigen::Index>(K), static_cast<Eigen::Index>(K));
  A.setFromTriplets(tr.begin(), tr.end());
  Eigen::SparseLU<Sp, Eigen::COLAMDOrdering<int>> lu;
  lu.compute(A);
  FG_REQUIRE(lu.info() == Eigen::Success, "cubic spline: singular system (nodes not distinct?)");
  const Eigen::VectorXd M = lu.solve(rhs);
  return {M.data(), M.data() + M.size()};
}

Dispersion::Dispersion(std::string pointer, Quantity q, std::variant<Sellmeier, Cauchy, DrudeLorentz, Table> model)
    : pointer_(std::move(pointer)), quantity_(q), model_(std::move(model)) {
  if (auto* t = std::get_if<Table>(&model_)) {
    const std::string where = pointer_ + "/table";
    FG_REQUIRE(t->lambda_um.size() >= 2, where + ": at least two rows needed");
    FG_REQUIRE(t->a.size() == t->lambda_um.size() && t->b.size() == t->lambda_um.size(), where + ": column sizes");
    for (std::size_t i = 0; i + 1 < t->lambda_um.size(); ++i)
      FG_REQUIRE(t->lambda_um[i] < t->lambda_um[i + 1], where + ": wavelengths must be strictly monotonic");
    FG_REQUIRE(t->lambda_um.front() > 0, where + ": wavelengths must be positive");
    if (t->row.empty())
      for (std::size_t i = 0; i < t->lambda_um.size(); ++i) t->row.push_back(static_cast<int>(i));
    if (t->cubic) {
      FG_REQUIRE(t->lambda_um.size() >= 4, where + ": \"cubic\" needs at least 4 rows");
      t->Ma = splineSecondDerivatives(t->lambda_um, t->a);
      t->Mb = splineSecondDerivatives(t->lambda_um, t->b);
    }
  }
  if (const auto* s = std::get_if<Sellmeier>(&model_))
    FG_REQUIRE(s->B.size() == s->C_um2.size(), pointer_ + "/sellmeier: B and C_um2 need the same length");
}

std::pair<double, double> Dispersion::validRange() const {
  if (const auto* t = std::get_if<Table>(&model_)) return {t->lambda_um.front(), t->lambda_um.back()};
  if (const auto* s = std::get_if<Sellmeier>(&model_)) return {s->rangeLo_um, s->rangeHi_um};
  return {0.0, std::numeric_limits<double>::infinity()};
}

void Dispersion::check(double lam) const {
  const std::string where = pointer_ + "/" + modelKey(model_);
  FG_REQUIRE(std::isfinite(lam) && lam > 0, where + ": lambda " + fmt(lam) + " um is not a positive wavelength");
  if (const auto* t = std::get_if<Table>(&model_)) {
    if (outside(*t, lam) && !t->extrapolateConstant)
      throw Error(where + ": lambda " + fmt(lam) + " um outside [" + fmt(t->lambda_um.front()) + ", " +
                  fmt(t->lambda_um.back()) + "] um");
  } else if (const auto* s = std::get_if<Sellmeier>(&model_)) {
    if (lam < s->rangeLo_um || lam > s->rangeHi_um)
      throw Error(where + ": lambda " + fmt(lam) + " um outside [" + fmt(s->rangeLo_um) + ", " + fmt(s->rangeHi_um) +
                  "] um (range_um)");
    for (std::size_t k = 0; k < s->C_um2.size(); ++k)
      if (std::abs(lam * lam - s->C_um2[k]) <= 1e-9 * lam * lam)
        throw Error(where + ": lambda " + fmt(lam) + " um at the pole sqrt(C_um2[" + std::to_string(k) + "]) = " +
                    fmt(std::sqrt(std::abs(s->C_um2[k]))) + " um");
  }
}

bool Dispersion::extrapolated(double lam) const {
  const auto* t = std::get_if<Table>(&model_);
  return t && t->extrapolateConstant && outside(*t, lam);
}

bool Dispersion::atLinearNode(double lam) const {
  const auto* t = std::get_if<Table>(&model_);
  return t && !t->cubic && std::binary_search(t->lambda_um.begin(), t->lambda_um.end(), lam);
}

cplx Dispersion::value(double lam) const {
  check(lam);
  return std::visit(
      [&](const auto& m) -> cplx {
        using T = std::decay_t<decltype(m)>;
        if constexpr (std::is_same_v<T, Sellmeier>) {  // Plan eq:disp_sellmeier
          const double l2 = lam * lam;
          double e = m.A;
          for (std::size_t k = 0; k < m.B.size(); ++k) e += m.B[k] * l2 / (l2 - m.C_um2[k]);
          return {e, 0.0};
        } else if constexpr (std::is_same_v<T, Cauchy>) {  // Plan eq:disp_cauchy
          const double l2 = lam * lam;
          return {m.A + m.B_um2 / l2 + m.C_um4 / (l2 * l2), m.k};
        } else if constexpr (std::is_same_v<T, DrudeLorentz>) {  // Plan eq:disp_dl
          const double w = HC_EV_UM / lam, wp2 = m.omegaP_eV * m.omegaP_eV;
          cplx e = m.epsInf - m.f0 * wp2 / (w * (w + I * m.gamma0_eV));
          for (const auto& o : m.oscillators) e += o.f * wp2 / (o.omega_eV * o.omega_eV - w * w - I * w * o.gamma_eV);
          return e;
        } else {
          if (outside(m, lam)) {
            const std::size_t j = lam < m.lambda_um.front() ? 0 : m.lambda_um.size() - 1;
            return {m.a[j], m.b[j]};
          }
          return {column(m, m.a, m.Ma, lam).first, column(m, m.b, m.Mb, lam).first};
        }
      },
      model_);
}

cplx Dispersion::eps(double lam) const {
  const cplx v = value(lam);
  return quantity_ == Quantity::N ? v * v : v;
}

cplx Dispersion::dValueDLambda(double lam) const {
  check(lam);
  return std::visit(
      [&](const auto& m) -> cplx {
        using T = std::decay_t<decltype(m)>;
        if constexpr (std::is_same_v<T, Sellmeier>) {
          // d/dl [B l^2 / (l^2 - C)] = -2 B C l / (l^2 - C)^2
          const double l2 = lam * lam;
          double d = 0;
          for (std::size_t k = 0; k < m.B.size(); ++k) {
            const double q = l2 - m.C_um2[k];
            d += -2 * m.B[k] * m.C_um2[k] * lam / (q * q);
          }
          return {d, 0.0};
        } else if constexpr (std::is_same_v<T, Cauchy>) {
          const double l2 = lam * lam;
          return {-2 * m.B_um2 / (l2 * lam) - 4 * m.C_um4 / (l2 * l2 * lam), 0.0};
        } else if constexpr (std::is_same_v<T, DrudeLorentz>) {
          const double w = HC_EV_UM / lam, dw = -HC_EV_UM / (lam * lam), wp2 = m.omegaP_eV * m.omegaP_eV;
          const cplx dd = w * w + I * m.gamma0_eV * w;
          cplx de = m.f0 * wp2 * (2 * w + I * m.gamma0_eV) / (dd * dd);
          for (const auto& o : m.oscillators) {
            const cplx q = o.omega_eV * o.omega_eV - w * w - I * w * o.gamma_eV;
            de += o.f * wp2 * (2 * w + I * o.gamma_eV) / (q * q);
          }
          return de * dw;
        } else {
          if (outside(m, lam)) return 0.0;
          return {column(m, m.a, m.Ma, lam).second, column(m, m.b, m.Mb, lam).second};
        }
      },
      model_);
}

cplx Dispersion::dEpsDLambda(double lam) const {
  const cplx d = dValueDLambda(lam);
  return quantity_ == Quantity::N ? 2.0 * value(lam) * d : d;
}

cplx Dispersion::dValueDParam(const std::vector<std::string>& t, double lam) const {
  check(lam);
  const std::string key = modelKey(model_);
  std::string path = pointer_;
  for (const auto& s : t) path += "/" + s;
  FG_REQUIRE(!t.empty() && t[0] == key, "gradient: '" + path + "' is not a value of the model '" + key + "'");
  const std::string bad = "gradient: '" + path + "' is not a parameter of the model '" + key + "'";
  const double l2 = lam * lam;
  if (const auto* s = std::get_if<Sellmeier>(&model_)) {
    if (t.size() == 2 && t[1] == "A") return 1.0;
    FG_REQUIRE(t.size() == 3, bad);
    const std::size_t k = index(t[2], s->B.size(), path);
    const double q = l2 - s->C_um2[k];
    if (t[1] == "B") return l2 / q;
    if (t[1] == "C_um2") return s->B[k] * l2 / (q * q);
    throw Error(bad);
  }
  if (std::holds_alternative<Cauchy>(model_)) {
    FG_REQUIRE(t.size() == 2, bad);
    if (t[1] == "A") return 1.0;
    if (t[1] == "B_um2") return 1.0 / l2;
    if (t[1] == "C_um4") return 1.0 / (l2 * l2);
    if (t[1] == "k") return I;
    throw Error(bad);
  }
  if (const auto* d = std::get_if<DrudeLorentz>(&model_)) {
    const double w = HC_EV_UM / lam, wp = d->omegaP_eV, wp2 = wp * wp;
    const cplx dd = w * (w + I * d->gamma0_eV);
    if (t.size() == 2) {
      if (t[1] == "eps_inf") return 1.0;
      if (t[1] == "f0") return -wp2 / dd;
      if (t[1] == "gamma0_eV") return d->f0 * wp2 * I * w / (dd * dd);
      if (t[1] == "omega_p_eV") {
        cplx s = -d->f0 / dd;
        for (const auto& o : d->oscillators) s += o.f / (o.omega_eV * o.omega_eV - w * w - I * w * o.gamma_eV);
        return 2 * wp * s;
      }
      throw Error(bad);
    }
    FG_REQUIRE(t.size() == 4 && t[1] == "oscillators", bad);
    const auto& o = d->oscillators[index(t[2], d->oscillators.size(), path)];
    const cplx q = o.omega_eV * o.omega_eV - w * w - I * w * o.gamma_eV;
    if (t[3] == "f") return wp2 / q;
    if (t[3] == "omega_eV") return -o.f * wp2 * 2 * o.omega_eV / (q * q);
    if (t[3] == "gamma_eV") return o.f * wp2 * I * w / (q * q);
    throw Error(bad);
  }
  const auto& tb = std::get<Table>(model_);
  FG_REQUIRE(tb.inlineData && t.size() == 4 && t[1] == "data",
             bad + " (table values are parameters only with inline \"data\")");
  const std::size_t r = index(t[2], tb.row.size(), path);
  const std::size_t c = index(t[3], 3, path);
  FG_REQUIRE(c != 0, "gradient: '" + path + "': table wavelengths are not parameters");
  const std::size_t node =
      static_cast<std::size_t>(std::find(tb.row.begin(), tb.row.end(), static_cast<int>(r)) - tb.row.begin());
  // the value is linear in the column: derivative = interpolant of the unit vector
  double w = 0;
  if (outside(tb, lam)) {
    w = (lam < tb.lambda_um.front() ? 0 : tb.lambda_um.size() - 1) == node ? 1.0 : 0.0;
  } else {
    std::vector<double> e(tb.lambda_um.size(), 0.0);
    e[node] = 1.0;
    const std::vector<double> M = tb.cubic ? splineSecondDerivatives(tb.lambda_um, e) : std::vector<double>{};
    w = column(tb, e, M, lam).first;
  }
  return c == 1 ? cplx(w) : I * w;
}

}  // namespace fg
