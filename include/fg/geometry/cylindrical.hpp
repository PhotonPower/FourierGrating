// FourierGrating - rotationally symmetric structures about the y axis (M20, Plan section
// "M20: Schichttypen, CAD und Konfiguration"). Lengths in micrometres.
#pragma once

#include <string>
#include <vector>

namespace fg {

struct CylindricalSpec {
  bool enabled = false;
  /// Azimuthal harmonics m; empty: those of the incidence (Plan sec:cylexc).
  std::vector<int> harmonics;
  // spectrum (Plan sec:cylgrid)
  double kappaMax = 6.0;  ///< k_max = kappa_max max Re sqrt(eps)
  double dkFactor = 1.0;  ///< panel width pi / (2 R_s dk_factor)
  int nodesPerPanel = 8;
  enum class Deform { Auto, Off, On } deform = Deform::Auto;
  double deltaRel = 0.05;  ///< depth delta = delta_rel (k_b - k_a)
  bool liR = true;         ///< factorization: Li in r (true) or Laurent
  // incidence: beam on the axis, waist in the superstrate (Plan sec:cylexc)
  double waist_um = 0.0;
  double focus_y_um = 0.0;  ///< height of the waist (y = 0 at the substrate interface)
  /// "x", "z", "circular_plus" (E_z = +i E_x), "circular_minus" (E_z = -i E_x), "radial",
  /// "azimuthal"
  std::string polarization = "x";
};

}  // namespace fg
