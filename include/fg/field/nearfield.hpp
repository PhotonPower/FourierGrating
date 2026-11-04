// FourierGrating - near-field reconstruction and cuts (Plan section "Nahfeld").
#pragma once

#include <string>
#include <vector>

#include "fg/core/types.hpp"

namespace fg {

enum class Reconstruction { Direct, FFF };

struct NearFieldRequest {
  double horizontal_y_um = 0;   ///< horizontal cut: y = const (1D: line E(x), 2D: map (x, z))
  double vertical_z_um = 0;     ///< vertical cut: z = const, map (x, y)
  double y_min_um = -0.2, y_max_um = 0.5;
  int nx = 400, ny = 300, nz = 200;
  int periods = 1;
  Reconstruction reconstruction = Reconstruction::FFF;
  std::vector<std::string> quantities{"E", "H", "absE2"};
};

/// Real-space field on a cut: complex components E(3), H(3) per grid point, row-major.
struct FieldCut {
  std::vector<double> axis1_um, axis2_um;   ///< axis2 empty for a 1D line
  std::vector<cplx> E[3], H[3];
  std::vector<double> materialIndex;        ///< for contour overlays
};

/// TODO(M6/M9): computed from a solved stack with storeInternal = true.
class NearFieldWriter {
 public:
  /// Writes one .npy per quantity plus meta.json into dir (Plan "Schnitte und Ausgabe").
  static void write(const FieldCut& cut, const std::string& dir, const std::string& name);
};

}  // namespace fg
