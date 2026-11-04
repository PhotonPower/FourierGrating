// FourierGrating - near-field reconstruction and cuts (Plan section "Nahfeld").
#pragma once

#include <string>
#include <vector>

#include "fg/core/types.hpp"

namespace fg {

enum class Reconstruction { Direct, FFF };

struct NearFieldRequest {
  double horizontal_y_um = 0;   ///< horizontal cut: y = const (1D: line u(x), 2D: map (x, z))
  double vertical_z_um = 0;     ///< vertical cut: z = const, map (x, y)
  double y_min_um = -0.2, y_max_um = 0.5;
  int nx = 400, ny = 300, nz = 200;
  int periods = 1;
  Reconstruction reconstruction = Reconstruction::FFF;
  std::vector<std::string> quantities{"E", "H", "absE2"};
};

/// Real-space field on a cut. Arrays are row-major with axis2 outer (size n2 * n1);
/// for a line (horizontal cut of a 1D grating) axis2 is empty and the arrays have n1 entries.
/// E in units of the incident amplitude (|E_inc| = 1), H normalized as Z0 * H.
struct FieldCut {
  std::string cut;                          ///< "horizontal" | "vertical"
  std::string axis1Name = "x", axis2Name = "y";
  std::vector<double> axis1_um, axis2_um;
  std::vector<cplx> E[3], H[3];
  std::vector<double> material;             ///< Re(eps_xx) for contour overlays
  double lambda_um = 0;
};

struct NearFieldResult {
  FieldCut horizontal, vertical;
};

class NearFieldWriter {
 public:
  /// Writes Ex..Hz (complex128), absE2, absH2, Sy (float64), material.npy and meta.json into dir
  /// (format read by scripts/plot_nearfield.py).
  static void write(const FieldCut& cut, const std::string& dir);
};

}  // namespace fg
