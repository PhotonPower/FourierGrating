// FourierGrating - PML zones and background stacks (Plan section "M18", eq:pml).
#include "fg/geometry/aperiodic.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <variant>

#include "fg/core/error.hpp"
#include "fg/core/quadrature.hpp"
#include "fg/geometry/profile.hpp"
#include "fg/geometry/section.hpp"
#include "fg/geometry/stack.hpp"

namespace fg {

PmlMap::PmlMap(double period, double width, double beta) : d_(period), w_(width), beta_(beta) {
  FG_REQUIRE(period > 0, "PML: period must be positive");
  FG_REQUIRE(width > 0 && 2 * width < period, "PML: width must satisfy 0 < width < period / 2");
  FG_REQUIRE(beta >= 0, "PML: beta must be non-negative");
}

cplx PmlMap::inverseDerivative(double X) const {
  X -= d_ * std::floor(X / d_);
  const double xi = X < w_ ? (w_ - X) / w_ : X > d_ - w_ ? (X - d_ + w_) / w_ : 0.0;
  if (xi <= 0) return 1.0;
  const double s = std::sin(0.5 * PI * xi), u = s * s;
  return (1.0 - u) / cplx(1.0, beta_ * u);  // Plan eq:pml
}

std::vector<cplx> PmlMap::inverseDerivativeCoefficients(int pmax) const {
  const double br[4] = {0.0, w_, d_ - w_, d_};
  const auto& g = quadrature::gaussLegendre(16);
  const double K = 2 * PI / d_;
  // panels: phase <= 1 rad at |p| = pmax; in the zones also <= w / 8 (u varies on the scale w)
  const double hp = d_ / (2 * PI * std::max(pmax, 1));
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1), 0.0);
  for (int k = 0; k < 3; ++k) {
    const double a = br[k], b = br[k + 1];
    const double hmax = k == 1 ? hp : std::min(hp, w_ / 8);
    const int m = std::max(1, static_cast<int>(std::ceil((b - a) / hmax)));
    const double h = (b - a) / m;
    for (int s = 0; s < m; ++s)
      for (std::size_t q = 0; q < g.x.size(); ++q) {
        const double X = a + (s + 0.5 * (1.0 + g.x[q])) * h;
        const cplx v = 0.5 * h * g.w[q] / d_ * inverseDerivative(X);
        for (int p = -pmax; p <= pmax; ++p) c[static_cast<std::size_t>(p + pmax)] += v * std::exp(-I * (p * K * X));
      }
  }
  return c;
}

namespace {

/// Zone test of one axis: [lo, hi] (micrometres) must lie in the physical region [w, d - w].
void requireInside(double lo, double hi, double w, double d, const std::string& what, char axis) {
  const double tol = 1e-12 * d;
  FG_REQUIRE(lo >= w - tol && hi <= d - w + tol,
             what + " reaches into the PML zone along " + std::string(1, axis) + " (structure in [" + std::to_string(lo) +
                 ", " + std::to_string(hi) + "] um, physical region [" + std::to_string(w) + ", " +
                 std::to_string(d - w) + "] um)");
}

void checkRegions(const std::vector<geometry::Region>& regions, const Stack& st, const std::string& what) {
  for (const auto& r : regions) {
    const auto b = r.bounds();  // {xmin, zmin, xmax, zmax} of the control points
    if (st.aperiodic.x.active()) requireInside(b[0], b[2], st.aperiodic.x.width_um, st.dx_um, what, 'x');
    if (st.aperiodic.z.active()) requireInside(b[1], b[3], st.aperiodic.z.width_um, st.dz_um, what, 'z');
  }
}

}  // namespace

void checkPmlZones(const Stack& st) {
  const auto& ap = st.aperiodic;
  if (!ap.enabled) return;
  const double wx = ap.x.width_um, d = st.dx_um;
  std::vector<std::pair<std::string, const LayerSpec*>> all;
  for (std::size_t k = 0; k < st.layers.size(); ++k) all.push_back({"layer " + std::to_string(k), &st.layers[k].spec});
  if (st.superstrateStructure) all.push_back({"structured superstrate", &st.superstrateStructure->spec});
  if (st.substrateStructure) all.push_back({"structured substrate", &st.substrateStructure->spec});
  for (const auto& [what, sp] : all) {
    const auto& spec = *sp;
    if (const auto* l = std::get_if<Lamellar1DLayer>(&spec)) {
      if (!ap.x.active()) continue;
      const std::size_t K = l->materials.size();
      for (std::size_t j = 0; j < K; ++j)
        if (l->materials[j] != l->materials[(j + K - 1) % K]) {
          const double x = l->breakpoints_um[j];
          requireInside(x, x, wx, d, what + " (lamellar material jump)", 'x');
        }
    } else if (const auto* p = std::get_if<Profile1DLayer>(&spec)) {
      if (!ap.x.active()) continue;
      const double g0 = profileValue(*p, 0.0, d);
      for (int s = 0; s <= 64; ++s)
        for (double x : {wx * s / 64.0, d - wx * s / 64.0})
          FG_REQUIRE(std::abs(profileValue(*p, x, d) - g0) <= 1e-12 * std::max(p->depth_um, 1e-300),
                     what + " (profile) is not flat in the PML zone along x");
    } else if (const auto* q = std::get_if<Polygon2DLayer>(&spec)) {
      for (const auto& inc : q->inclusions) {
        double x0 = 1e300, x1 = -1e300, z0 = 1e300, z1 = -1e300;
        for (const auto& v : inc.vertices) {
          x0 = std::min(x0, v[0]);
          x1 = std::max(x1, v[0]);
          z0 = std::min(z0, v[1]);
          z1 = std::max(z1, v[1]);
        }
        if (ap.x.active()) requireInside(x0, x1, wx, d, what + " (polygon)", 'x');
        if (ap.z.active()) requireInside(z0, z1, ap.z.width_um, st.dz_um, what + " (polygon)", 'z');
      }
    } else if (const auto* c = std::get_if<Contour2DLayer>(&spec)) {
      checkRegions(c->regions, st, what + " (contour2d)");
    } else if (const auto* s = std::get_if<SectionedLayer>(&spec)) {
      const double t = s->provider->thickness_um();
      for (int i = 0; i < 9; ++i) checkRegions(s->provider->section(t * (i + 0.5) / 9).regions, st, what + " (cad)");
    } else if (std::holds_alternative<Profile2DLayer>(spec)) {
      throw Error(what + ": profile2d layers are not supported with PML (use contour2d or cad)");
    }
  }
}

Stack backgroundStack(const Stack& st) {
  Stack bg = st;
  bg.layers.clear();
  auto homogeneous = [&](double t, const std::string& mat) {
    if (t > 0) bg.layers.push_back({HomogeneousLayer{t, mat}, Factorization::Auto});
  };
  for (const auto& le : st.layers) {
    const auto& spec = le.spec;
    if (std::holds_alternative<HomogeneousLayer>(spec)) {
      bg.layers.push_back(le);
    } else if (const auto* l = std::get_if<Lamellar1DLayer>(&spec)) {
      std::size_t j = 0;  // interval containing x = 0 (the zone at the cell boundary)
      while (j + 1 < l->materials.size() && l->breakpoints_um[j + 1] <= 0) ++j;
      homogeneous(l->thickness_um, l->materials[j]);
    } else if (const auto* p = std::get_if<Profile1DLayer>(&spec)) {
      const double g0 = profileValue(*p, 0.0, st.dx_um);
      homogeneous(g0, p->below);
      homogeneous(p->depth_um - g0, p->above);
    } else if (const auto* q = std::get_if<Polygon2DLayer>(&spec)) {
      homogeneous(q->thickness_um, q->background);
    } else if (const auto* c = std::get_if<Contour2DLayer>(&spec)) {
      homogeneous(c->thickness_um, c->background);
    } else if (const auto* s = std::get_if<SectionedLayer>(&spec)) {
      homogeneous(s->provider->thickness_um(), s->background);
    } else if (const auto* r = std::get_if<RingsLayer>(&spec)) {
      homogeneous(r->thickness_um, r->materials.back());  // M20: outside the rings
    } else {
      throw Error("background stack: profile2d layers are not supported with PML");
    }
  }
  return bg;
}

}  // namespace fg
