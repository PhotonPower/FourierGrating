// FourierGrating - rings_profile layers (Plan section "M20: Schichttypen"): radii as functions of y and their
// RCWA staircase of RingsLayers.
#include <algorithm>
#include <cmath>
#include <sstream>

#include "fg/core/error.hpp"
#include "fg/geometry/layer.hpp"

namespace fg {

double RingsProfileLayer::Radius::at(double y, double thickness) const {
  double r = 0;
  switch (shape) {
    case Shape::Constant: r = r0_um; break;
    case Shape::Cone: r = r0_um + (r1_um - r0_um) * y / thickness; break;
    case Shape::Sphere: {
      const double d = y - centerY_um;
      r = std::abs(d) < sphereR_um ? std::sqrt(sphereR_um * sphereR_um - d * d) : 0.0;
      break;
    }
    case Shape::Table: {
      FG_REQUIRE(ys_um.size() >= 2 && ys_um.size() == rs_um.size(), "rings_profile: table needs >= 2 points (y_um, r_um)");
      if (y <= ys_um.front()) return rs_um.front();
      if (y >= ys_um.back()) return rs_um.back();
      const auto k = static_cast<std::size_t>(std::upper_bound(ys_um.begin(), ys_um.end(), y) - ys_um.begin());
      const double t = (y - ys_um[k - 1]) / (ys_um[k] - ys_um[k - 1]);
      r = rs_um[k - 1] + t * (rs_um[k] - rs_um[k - 1]);
      break;
    }
    case Shape::Function:
      FG_REQUIRE(static_cast<bool>(f), "rings_profile: radius function missing");
      r = f(y);
      break;
  }
  FG_REQUIRE(std::isfinite(r) && r >= 0, "rings_profile: radius must be finite and >= 0");
  return r;
}

std::vector<LayerEntry> ringsStaircase(const RingsProfileLayer& l, int defaultSlices) {
  FG_REQUIRE(l.thickness_um > 0, "rings_profile: thickness_um > 0");
  FG_REQUIRE(!l.radii.empty() && l.materials.size() == l.radii.size() + 1,
             "rings_profile: n radii and n + 1 materials (the last one outside)");
  for (const auto& r : l.radii) {
    if (r.shape == RingsProfileLayer::Radius::Shape::Sphere)
      FG_REQUIRE(r.sphereR_um > 0, "rings_profile: sphere radius_um > 0");
    if (r.shape == RingsProfileLayer::Radius::Shape::Table)
      FG_REQUIRE(r.ys_um.size() >= 2 && r.ys_um.size() == r.rs_um.size() &&
                     std::adjacent_find(r.ys_um.begin(), r.ys_um.end(), std::greater_equal<>()) == r.ys_um.end(),
                 "rings_profile: table y_um strictly increasing, as many r_um, at least 2 points");
  }
  const int K = l.slices > 0 ? l.slices : defaultSlices;
  FG_REQUIRE(K >= 1, "rings_profile: slices >= 1");
  const double h = l.thickness_um / K;
  std::vector<LayerEntry> out;
  for (int k = 0; k < K; ++k) {
    const double y = (k + 0.5) * h;
    RingsLayer rl;
    rl.thickness_um = h;
    double prev = 0;
    for (std::size_t i = 0; i < l.radii.size(); ++i) {
      const double r = l.radii[i].at(y, l.thickness_um), tol = 1e-12 * std::max(1.0, r);
      if (r < prev - tol) {
        std::ostringstream m;
        m << "rings_profile: radius " << i + 1 << " (" << r << " um) is below radius " << i << " (" << prev
          << " um) at y = " << y << " um (radii must not cross)";
        throw Error(m.str());
      }
      if (r <= prev + tol) continue;  // ring i has no width here
      rl.radii_um.push_back(r);
      rl.materials.push_back(l.materials[i]);
      prev = r;
    }
    LayerEntry e;
    if (rl.radii_um.empty()) {
      e.spec = HomogeneousLayer{h, l.materials.back()};
    } else {
      rl.materials.push_back(l.materials.back());
      e.spec = std::move(rl);
    }
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace fg
