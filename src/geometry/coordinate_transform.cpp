// FourierGrating - coordinate transformations (Plan section "M17", eq:ct_map, eq:asr).
#include "fg/geometry/coordinate_transform.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/geometry/stack.hpp"

namespace fg {

namespace {

double wrap(double x, double d) {
  double r = std::fmod(x, d);
  if (r < 0) r += d;
  return r >= d ? 0.0 : r;
}

}  // namespace

CoordinateMap::CoordinateMap(double period, std::vector<double> breakpoints, double eta) : d_(period), eta_(eta) {
  FG_REQUIRE(period > 0, "coordinate transform: period must be positive");
  FG_REQUIRE(eta >= 0 && eta < 1, "coordinate transform: eta must satisfy 0 <= eta < 1");
  FG_REQUIRE(!breakpoints.empty(), "coordinate transform: at least one breakpoint");
  for (double& b : breakpoints) b = wrap(b, period);
  std::sort(breakpoints.begin(), breakpoints.end());
  const double tol = 1e-12 * period;
  for (double b : breakpoints)
    if (br_.empty() || b - br_.back() > tol) br_.push_back(b);
  if (br_.size() > 1 && br_.front() + period - br_.back() <= tol) br_.pop_back();
}

void CoordinateMap::locate(double X, std::size_t& j, double& shift, double& u) const {
  const double b0 = br_.front();
  const double k = std::floor((X - b0) / d_);
  shift = k * d_;
  double Xr = X - shift;
  if (Xr >= b0 + d_) {  // rounding
    Xr -= d_;
    shift += d_;
  }
  j = static_cast<std::size_t>(std::upper_bound(br_.begin(), br_.end(), Xr) - br_.begin());
  j = j == 0 ? 0 : j - 1;
  u = Xr - br_[j];
}

double CoordinateMap::f(double X) const {
  std::size_t j;
  double shift, u;
  locate(X, j, shift, u);
  const double D = (j + 1 < br_.size() ? br_[j + 1] : br_.front() + d_) - br_[j];
  return shift + br_[j] + u - eta_ * D / (2 * PI) * std::sin(2 * PI * u / D);  // Plan eq:asr
}

double CoordinateMap::df(double X) const {
  std::size_t j;
  double shift, u;
  locate(X, j, shift, u);
  const double D = (j + 1 < br_.size() ? br_[j + 1] : br_.front() + d_) - br_[j];
  return 1.0 - eta_ * std::cos(2 * PI * u / D);
}

double CoordinateMap::inverse(double x) const {
  // f maps every interval [x_{j-1}, x_j] onto itself: solve u - c sin(w u) = t on [0, Delta]
  std::size_t j;
  double shift, t;
  locate(x, j, shift, t);
  const double D = (j + 1 < br_.size() ? br_[j + 1] : br_.front() + d_) - br_[j];
  const double c = eta_ * D / (2 * PI), w = 2 * PI / D;
  double lo = 0, hi = D, u = t;
  for (int it = 0; it < 100; ++it) {
    const double g = u - c * std::sin(w * u) - t;
    (g > 0 ? hi : lo) = u;
    const double dg = 1.0 - eta_ * std::cos(w * u);
    double un = u - g / dg;
    if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);  // bisection step
    const double du = un - u;
    u = un;
    if (std::abs(du) <= 1e-16 * D) break;
  }
  return shift + br_[j] + u;
}

CoordinateMap::Nodes CoordinateMap::nodes(int pmax, const std::vector<double>& extraBreaks) const {
  const double b0 = br_.front();
  std::vector<double> br(br_.begin(), br_.end());
  for (double e : extraBreaks) br.push_back(b0 + wrap(e - b0, d_));
  std::sort(br.begin(), br.end());
  br.push_back(b0 + d_);
  const auto& g = quadrature::gaussLegendre(16);
  // panel width: phase <= 1 rad at |p| = pmax, and at most half the distance of the nearest pole of
  // 1/f' (cos(2 pi u / Delta) = 1/eta) from the real axis
  const double hp = d_ / (2 * PI * std::max(pmax, 1));
  Nodes out;
  for (std::size_t k = 0; k + 1 < br.size(); ++k) {
    const double a = br[k], b = br[k + 1];
    if (!(b > a)) continue;
    std::size_t j;
    double shift, u;
    locate(0.5 * (a + b), j, shift, u);
    const double D = (j + 1 < br_.size() ? br_[j + 1] : br_.front() + d_) - br_[j];
    const double hmax = eta_ > 0 ? std::min(hp, 0.5 * D * std::acosh(1.0 / eta_) / (2 * PI)) : hp;
    const int m = std::max(1, static_cast<int>(std::ceil((b - a) / hmax)));
    const double h = (b - a) / m;
    for (int s = 0; s < m; ++s)
      for (std::size_t q = 0; q < g.x.size(); ++q) {
        out.X.push_back(a + (s + 0.5 * (1.0 + g.x[q])) * h);
        out.w.push_back(0.5 * h * g.w[q] / d_);
      }
  }
  return out;
}

std::vector<cplx> CoordinateMap::inverseDerivativeCoefficients(int pmax) const {
  const Nodes n = nodes(pmax);
  const double K = 2 * PI / d_;
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1), 0.0);
  for (std::size_t l = 0; l < n.X.size(); ++l) {
    const double v = n.w[l] / df(n.X[l]);
    for (int p = -pmax; p <= pmax; ++p) c[static_cast<std::size_t>(p + pmax)] += v * std::exp(-I * (p * K * n.X[l]));
  }
  return c;
}

MatC CoordinateMap::toPhysical(double alpha0, int N) const {
  const double K = 2 * PI / d_;
  const Nodes q = nodes(2 * N + 2 + static_cast<int>(std::ceil(std::abs(alpha0) / K)));
  const auto L = static_cast<Eigen::Index>(q.X.size());
  MatC Ex(L, 2 * N + 1), Ef(L, 2 * N + 1);
  for (Eigen::Index l = 0; l < L; ++l) {
    const double X = q.X[static_cast<std::size_t>(l)], fx = f(X);
    const double wl = q.w[static_cast<std::size_t>(l)] * df(X);
    for (int n = -N; n <= N; ++n) {
      const double a = alpha0 + n * K;
      Ex(l, n + N) = std::exp(I * (a * X));
      Ef(l, n + N) = wl * std::exp(-I * (a * fx));
    }
  }
  return Ef.transpose() * Ex;  // Plan eq:ct_rayleigh
}

MatC CoordinateMap::fromPhysical(double alpha0, int N) const {
  const double K = 2 * PI / d_;
  const Nodes q = nodes(2 * N + 2 + static_cast<int>(std::ceil(std::abs(alpha0) / K)));
  const auto L = static_cast<Eigen::Index>(q.X.size());
  MatC Ex(L, 2 * N + 1), Ef(L, 2 * N + 1);
  for (Eigen::Index l = 0; l < L; ++l) {
    const double X = q.X[static_cast<std::size_t>(l)], fx = f(X);
    const double wl = q.w[static_cast<std::size_t>(l)];
    for (int n = -N; n <= N; ++n) {
      const double a = alpha0 + n * K;
      Ex(l, n + N) = wl * std::exp(-I * (a * X));
      Ef(l, n + N) = std::exp(I * (a * fx));
    }
  }
  return Ex.transpose() * Ef;  // Plan eq:ct_inc
}

namespace detail {

void mappedStarts(const CoordinateMap& map, const std::vector<double>& breaks, std::vector<double>& starts,
                  std::vector<std::size_t>& index) {
  const double d = map.period();
  FG_REQUIRE(breaks.size() >= 2, "mapIntervals: need at least one interval");
  // (start in [0, d), physical interval) of every interval with positive width
  std::vector<std::pair<double, std::size_t>> s;
  for (std::size_t k = 0; k + 1 < breaks.size(); ++k)
    if (breaks[k + 1] > breaks[k]) s.emplace_back(wrap(map.inverse(breaks[k]), d), k);
  std::stable_sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  // value at X = 0: interval of the last start (cyclic) unless one starts at 0
  auto at = [&](double X) {
    std::size_t k = s.back().second;
    for (const auto& e : s)
      if (e.first <= X) k = e.second;
    return k;
  };
  std::vector<std::pair<double, std::size_t>> all;
  if (s.front().first > 0) all.emplace_back(0.0, at(0.0));
  for (const auto& e : s) all.push_back(e);
  for (double b : map.breakpoints())
    if (b > 0) all.emplace_back(b, at(b));
  std::stable_sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  starts.clear();
  index.clear();
  for (std::size_t i = 0; i < all.size(); ++i) {
    // a later start at the same X wins (zero-width interval dropped)
    if (i + 1 < all.size() && all[i + 1].first <= all[i].first) continue;
    starts.push_back(all[i].first);
    index.push_back(all[i].second);
  }
  if (starts.front() > 0) {  // only possible if the first start was dropped
    starts.insert(starts.begin(), 0.0);
    index.insert(index.begin(), at(0.0));
  }
}

}  // namespace detail

namespace {

/// Material jumps along x (axis 0) or z (axis 1) of all layers (Plan "Standard: alle Materialspruenge").
std::vector<double> autoBreakpoints(const Stack& st, int axis) {
  std::vector<double> out;
  for (const auto& le : st.layers) {
    if (const auto* l = std::get_if<Lamellar1DLayer>(&le.spec)) {
      if (axis != 0) continue;
      const std::size_t K = l->materials.size();
      for (std::size_t k = 0; k < K; ++k)
        if (l->materials[k] != l->materials[(k + K - 1) % K]) out.push_back(l->breakpoints_um[k]);
    } else if (const auto* p = std::get_if<Polygon2DLayer>(&le.spec)) {
      for (const auto& inc : p->inclusions)
        for (const auto& v : inc.vertices) out.push_back(v[static_cast<std::size_t>(axis)]);
    } else if (const auto* c = std::get_if<Contour2DLayer>(&le.spec)) {
      auto add = [&](const geometry::Contour2D& con) {
        for (const auto& s : con.segments)
          if (s.degree() == 1) out.push_back(s.P.front()[static_cast<std::size_t>(axis)]);
      };
      for (const auto& r : c->regions) {
        add(r.outer);
        for (const auto& h : r.holes) add(h);
      }
    }
  }
  return out;
}

std::shared_ptr<const CoordinateMap> resolveAxis(const Stack& st, const AxisTransformSpec& a, int axis, bool force) {
  if (a.type == AxisTransformSpec::Type::None) return nullptr;
  FG_REQUIRE(a.eta >= 0 && a.eta < 1, "coordinate_transform: eta must satisfy 0 <= eta < 1");
  if (a.eta == 0 && !force) return nullptr;  // identity: the diagonal path (bit-identical)
  const char* name = axis == 0 ? "coordinate_transform.x" : "coordinate_transform.z";
  FG_REQUIRE(axis == 0 || st.is2D(), std::string(name) + " needs a 2D stack");
  std::vector<double> br = a.breakpoints_um.empty() ? autoBreakpoints(st, axis) : a.breakpoints_um;
  FG_REQUIRE(!br.empty(), std::string(name) + ": no material jumps found along the axis, give breakpoints_um");
  return std::make_shared<const CoordinateMap>(axis == 0 ? st.dx_um : st.dz_um, std::move(br), a.eta);
}

}  // namespace

ResolvedTransform resolveTransform(const Stack& stack) {
  ResolvedTransform t;
  t.force = stack.transform.forceMatrixPath;
  t.x = resolveAxis(stack, stack.transform.x, 0, t.force);
  t.z = resolveAxis(stack, stack.transform.z, 1, t.force);
  return t;
}

}  // namespace fg
