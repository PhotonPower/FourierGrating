// FourierGrating - polynomials in Bernstein form (see the header).
#include "fg/geometry/bernstein.hpp"

#include <algorithm>

namespace fg::geometry {

namespace {

/// Coefficients on [0, 1/2] and [1/2, 1] (de Casteljau at 1/2).
void split(const std::vector<double>& c, std::vector<double>& l, std::vector<double>& r) {
  std::vector<double> w = c;
  const std::size_t n = c.size();
  l.assign(n, 0.0);
  r.assign(n, 0.0);
  l[0] = w[0];
  r[n - 1] = w[n - 1];
  for (std::size_t k = 1; k < n; ++k) {
    for (std::size_t j = 0; j + k < n; ++j) w[j] = 0.5 * (w[j] + w[j + 1]);
    l[k] = w[0];
    r[n - 1 - k] = w[n - 1 - k];
  }
}

int signChanges(const std::vector<double>& c) {
  int v = 0, last = 0;
  for (double x : c) {
    const int s = x > 0 ? 1 : x < 0 ? -1 : 0;
    if (s == 0) continue;
    if (last != 0 && s != last) ++v;
    last = s;
  }
  return v;
}

/// Simple roots of `full` inside (a, b), whose coefficients on [a, b] are c.
void isolate(const std::vector<double>& full, const std::vector<double>& c, double a, double b, int depth,
             std::vector<double>& roots) {
  const int v = signChanges(c);
  if (v == 0) return;
  if (v == 1 && c.front() != 0 && c.back() != 0) {
    double lo = a, hi = b, flo = c.front();
    for (int it = 0; it < 200 && hi - lo > 1e-17; ++it) {
      const double m = 0.5 * (lo + hi), fm = bernsteinValue(full, m);
      if (fm == 0) {
        lo = hi = m;
        break;
      }
      if ((fm > 0) == (flo > 0)) {
        lo = m;
        flo = fm;
      } else {
        hi = m;
      }
    }
    roots.push_back(0.5 * (lo + hi));
    return;
  }
  if (depth > 60) return;  // a multiple root without a sign change, not separable
  std::vector<double> l, r;
  split(c, l, r);
  const double m = 0.5 * (a + b);
  isolate(full, l, a, m, depth + 1, roots);
  if (l.back() == 0) roots.push_back(m);  // root exactly at the split point
  isolate(full, r, m, b, depth + 1, roots);
}

}  // namespace

double binom(int n, int k) {
  double r = 1;
  for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
  return r;
}

double bernsteinValue(std::vector<double> c, double t) {
  for (std::size_t r = 1; r < c.size(); ++r)
    for (std::size_t j = 0; j + r < c.size(); ++j) c[j] = (1 - t) * c[j] + t * c[j + 1];
  return c.front();
}

std::vector<double> bernsteinProduct(const std::vector<double>& a, const std::vector<double>& b) {
  const int n = static_cast<int>(a.size()) - 1, m = static_cast<int>(b.size()) - 1;
  std::vector<double> c(static_cast<std::size_t>(n + m + 1), 0.0);
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= m; ++j)
      c[static_cast<std::size_t>(i + j)] += binom(n, i) * binom(m, j) / binom(n + m, i + j) *
                                            a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(j)];
  return c;
}

std::vector<double> bernsteinDerivative(const std::vector<double>& c) {
  const std::size_t n = c.size() - 1;
  if (n == 0) return {0.0};
  std::vector<double> d(n);
  for (std::size_t j = 0; j < n; ++j) d[j] = static_cast<double>(n) * (c[j + 1] - c[j]);
  return d;
}

std::vector<double> bernsteinRoots(const std::vector<double>& c) {
  std::vector<double> t;
  isolate(c, c, 0.0, 1.0, 0, t);
  std::sort(t.begin(), t.end());
  return t;
}

}  // namespace fg::geometry
