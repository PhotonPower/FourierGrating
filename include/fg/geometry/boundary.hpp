// FourierGrating - boundary samples of material interfaces (Plan section "M13: Normalenfeld mit
// Randstichproben"): input of Algorithm NF, delivered by contours (Contour2D) or section providers.
#pragma once

#include <vector>

#include "fg/core/types.hpp"

namespace fg {

/// Boundary sample of Algorithm NF (Plan section "M13: Normalenfeld mit Randstichproben"): point
/// (x, z) on a material boundary and unit normal (3D; y-invariant layers: n_y = 0), oriented out of
/// the region it belongs to. exactChord: the chord to the next sample is the boundary itself
/// (polygon edge, normal of this sample along the whole chord); otherwise it approximates a curve
/// with sagitta <= 1e-3 raster spacings (normals interpolated between the two samples).
struct BoundarySample {
  double x = 0, z = 0;
  Vec3D n = Vec3D(1, 0, 0);
  bool exactChord = false;
};
/// Consecutive samples of one boundary curve (closed: chord from the last to the first sample).
struct BoundaryChain {
  std::vector<BoundarySample> samples;
  bool closed = true;
};

}  // namespace fg
