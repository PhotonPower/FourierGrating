// FourierGrating - cross-sections that change with the height (Plan section "M13: Hoehenabhaengige
// Schichtmodelle"); implemented by the CAD slicer (M14) and by analytic test providers.
#pragma once

#include <vector>

#include "fg/geometry/boundary.hpp"
#include "fg/geometry/contour.hpp"

namespace fg {

/// Cross-section of a layer at one height: regions (materials) in the background of the layer and
/// boundary samples with 3D normals (tilted walls: n_y != 0). Empty chains: sampled from the region
/// contours with n_y = 0.
struct CrossSection {
  std::vector<geometry::Region> regions;
  std::vector<BoundaryChain> chains;
};

class SectionProvider {
 public:
  virtual ~SectionProvider() = default;
  /// Cross-section at y (micrometres, relative to the layer bottom, 0 < y < thickness).
  virtual CrossSection section(double y_um) const = 0;
  virtual double thickness_um() const = 0;
};

}  // namespace fg
