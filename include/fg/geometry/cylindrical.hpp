// FourierGrating - rotationally symmetric structures about the y axis (M20, Plan section
// "M20: Schichttypen, CAD und Konfiguration"). Lengths in micrometres.
#pragma once

#include <string>
#include <vector>

namespace fg {

struct CylindricalSpec {
  bool enabled = false;
  /// Azimuthal harmonics m; empty: those of the incidence (Plan sec:cylexc). With rotation > 0: the classes of
  /// these harmonics are computed.
  std::vector<int> harmonics;
  /// P-fold rotational symmetry (M21, Plan sec:cylP): 0 = axial (M20, one harmonic at a time); P >= 1: classes
  /// m0 + l P, |l| <= classL, with polygon2d / contour2d layers (Laurent there).
  int rotation = 0;
  /// Truncation L of a class; < 0: ceil(m_cut / P) + 2 (Plan sec:cylP "Anregung und Symmetrieklassen").
  int classL = -1;
  // spectrum (Plan sec:cylgrid)
  double kappaMax = 6.0;  ///< k_max = kappa_max max Re sqrt(eps)
  double dkFactor = 1.0;  ///< panel width pi / (2 R_s dk_factor)
  int nodesPerPanel = 8;
  enum class Deform { Auto, Off, On } deform = Deform::Auto;
  double deltaRel = 0.05;  ///< depth delta = delta_rel (k_b - k_a)
  bool liR = true;         ///< factorization: Li in r (true) or Laurent
  /// Plane wave: evaluate the optical theorem as a control of sigma_ext (volume integrals over the rings,
  /// Plan v2.25); with lossy rings they are computed anyway (sigma_abs).
  bool opticalTheorem = false;
  /// Far field (Plan "M20: Ausgaben"): number of azimuths phi_k in [0, 360) deg; 0: none (the angular spectrum
  /// per harmonic is always computed).
  int farFieldPhi = 0;
  /// Plane wave (scattered field, angles of the scan) instead of the beam (Plan sec:cylexc).
  bool planeWave = false;
  // incidence: beam on the axis, waist in the superstrate (Plan sec:cylexc)
  double waist_um = 0.0;
  double focus_y_um = 0.0;  ///< height of the waist (y = 0 at the substrate interface)
  /// "x", "z", "circular_plus" (E_z = +i E_x), "circular_minus" (E_z = -i E_x), "radial",
  /// "azimuthal"
  std::string polarization = "x";
};

}  // namespace fg
