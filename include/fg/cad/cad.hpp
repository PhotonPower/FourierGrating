// FourierGrating - CAD slicer: STEP bodies -> layers (Plan section "M14: CAD-Slicer").
// The STEP geometry is read with Nurbforge (FG_ENABLE_NURBFORGE); only src/cad includes its
// headers.
#pragma once

#include <array>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "fg/geometry/layer.hpp"
#include "fg/geometry/stack.hpp"

namespace fg::cad {

/// Material of a STEP body, by solid name or by index (order of the bodies in the file).
struct BodyAssignment {
  std::string name;
  int index = -1;
  std::string material;
};

/// The "cad" layer block of the configuration (docs/cad.md).
struct CadOptions {
  std::string file;      ///< STEP file (absolute or relative to the working directory)
  std::string stepText;  ///< STEP contents; used instead of `file` if not empty
  /// Unit of the coordinates in the file: "auto" (as declared), "mm", "um" (overrides a wrong
  /// declaration)
  std::string unit = "auto";
  /// CAD axes ("+X", "-Y", ...) mapped to FG y (stacking direction) and x; FG z = x cross y. A
  /// given z must equal that (mirrored, left-handed axes are rejected: they change chiral
  /// structures).
  std::string axisY = "+Z", axisX = "+X", axisZ;
  std::optional<std::array<double, 2>> yRange_um;  ///< default: bounding box of the bodies
  std::optional<std::array<double, 2>>
      cellOrigin_um;                       ///< (x0, z0); default: minimum of the bounding box
  std::string background = "superstrate";  ///< material name, "superstrate" or "substrate"
  std::vector<BodyAssignment> bodies;
  std::string overlap = "error";  ///< "error" | "priority" (earlier bodies win, Boolean difference)
  bool periodicUnion = false;     ///< unite bodies with their periodic images before slicing
  double criticalTol_um = 1e-6;   ///< critical heights closer than this are merged
  double minBand_um = 1e-5;       ///< thinner bands are merged into the neighbour (warning)
  int rcwaSlices = 40;            ///< RCWA staircase slices of y-variable bands
  std::string yInterpolation =
      "none";  ///< "none" | "chebyshev" (Q blocks at Chebyshev nodes per band)
  int chebyshevNodes = 16;
  Factorization factorization = Factorization::Auto;
  NormalFieldOptions normalField;
  double lambda_um = 0;  ///< smallest wavelength of the run (UNIT_SUSPECT check; 0: unknown)
};

/// One finding of the geometry check (Plan section "Geometriepruefung fuer das Sprachmodell"):
/// fixed code, text, position in FG coordinates (micrometres; NaN if not applicable).
struct Finding {
  std::string code, message;
  std::array<double, 3> where_um{std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN()};
};

/// One band [y0, y1] of the decomposition and the layer type it becomes.
struct Band {
  double y0_um = 0, y1_um = 0;
  std::string kind;      ///< "homogeneous" | "extrusion" (Contour2D) | "sectioned" (SectionedLayer)
  std::string material;  ///< homogeneous bands
  std::vector<std::pair<std::string, double>>
      areaFractions;  ///< per material at mid height (incl. background)
};

struct Report {
  bool ok = true;
  std::vector<Finding> errors, warnings;
  std::vector<Band> bands;
  std::array<double, 2> cellOrigin_um{0, 0};
  /// JSON (fgsolve --cad-report, Python fg.cad.check).
  std::string json() const;
};

/// True if this build contains the slicer (FG_ENABLE_NURBFORGE).
bool available();

/// Checks the geometry and decomposes it into bands; geometry problems are findings, not
/// exceptions. context: periods, materials and the names of super- and substrate.
Report check(const CadOptions& opt, const Stack& context);

/// The layers (bottom to top) replacing the cad block; throws fg::Error listing the error codes.
std::vector<LayerEntry> slice(const CadOptions& opt, const Stack& context,
                              Report* report = nullptr);

/// SVG of the section at height y (FG coordinates, micrometres) over one cell, materials coloured.
std::string previewSvg(const CadOptions& opt, const Stack& context, double y_um);

/// Tool description (JSON schema) of the check for an LLM tool loop (e.g. Nurbforge's Ollama
/// agent).
std::string llmToolSpec();

}  // namespace fg::cad
