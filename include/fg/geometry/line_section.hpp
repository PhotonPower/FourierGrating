// FourierGrating - line sections of 1D profiles (Plan section "M23: Linienschnitt", "Normalenfeld im
// Linienschnitt"): crossings of a horizontal line with rational Bezier profiles in the (x, y) plane,
// intervals with materials (parity rule, periodic folding, overlaps as interval operations), the
// normal field with linearly interpolated angle (Plan eq:cad1d_N). Independent of Nurbforge.
#pragma once

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/core/types.hpp"
#include "fg/geometry/contour.hpp"

namespace fg {

struct NormalFieldCoefficients;
class CoordinateMap;

/// Section of a 1D layer at one height: 0 = breaks[0] < ... < breaks[K] = dx, materials[k] on
/// [breaks[k], breaks[k+1]], normals[k] at breaks[k] (k < K; breaks[K] is breaks[0] + dx). breaks[k] is
/// a material boundary iff materials[k] != materials[k - 1] (cyclic); other normals are unused.
struct LineSection {
  std::vector<double> breaks;
  std::vector<std::string> materials;
  std::vector<Vec3D> normals;  ///< in-plane (n_x, n_y, 0), sign arbitrary (projector form)
};

/// Height-dependent 1D layer (Plan section "M23: Schichttyp Sectioned1DLayer"): CAD slicer and
/// analytic test providers.
class LineSectionProvider {
 public:
  virtual ~LineSectionProvider() = default;
  /// Section at y (micrometres, relative to the layer bottom, 0 < y < thickness).
  virtual LineSection section(double y_um) const = 0;
  virtual double thickness_um() const = 0;
};

namespace geometry {

/// Crossing of a profile with a horizontal line: x and the unit normal (Plan eq:cad1d_n).
struct Crossing {
  double x = 0;
  Vec3D n = Vec3D(1, 0, 0);
};

/// Crossings of a rational Bezier segment (points (x, y)) with the line y = c for t in [0, 1) (the end
/// t = 1 is the start of the next segment of a closed chain): roots of the Bernstein polynomial
/// sum_j B_j w_j (y_j - c), isolated by subdivision (sign changes of the coefficients), refined by
/// safeguarded Newton. Tangential crossings (double roots) are not separated; they lie on critical
/// heights, which are never sampled.
void lineCrossings(const BezierSegment& s, double c, std::vector<Crossing>& out);

/// Profile of one body in the (x, y) plane (x relative to the cell origin, y absolute): closed loops
/// in any orientation and nesting (parity rule).
struct ProfileBody {
  std::vector<Contour2D> loops;
  std::string material;
  std::string name;
};

struct LineSectionOptions {
  bool periodicUnion = false;  ///< overlapping periodic images of one body are united
  bool priority = false;       ///< overlapping bodies: earlier bodies win (B_i <- B_i \ U_{j<i} B_j)
  double tol = 1e-10;          ///< length tolerance relative to dx (touching, empty intervals)
};

/// Overlaps found by lineSection (the slicer's PERIODIC_OVERLAP / BODY_OVERLAP).
struct LineSectionFindings {
  std::vector<int> selfOverlap;               ///< bodies overlapping their own periodic image
  std::vector<std::array<int, 2>> bodyOverlap;  ///< (i, j), j < i, overlapping intervals
};

/// Line section y = const of the bodies in the cell [0, dx) (Plan section "M23: Linienschnitt"): parity
/// rule per body, folding into the cell; touching images of one body are merged without a boundary
/// (seam), overlaps are united (periodicUnion) or reported; overlapping bodies are resolved by priority
/// or reported; the background fills the gaps; equal neighbours merge. findings == nullptr: unresolved
/// overlaps throw; otherwise they are recorded and resolved as with periodicUnion / priority.
LineSection lineSection(const std::vector<ProfileBody>& bodies, double y, double dx, const std::string& background,
                        const LineSectionOptions& opt, LineSectionFindings* findings = nullptr);

/// Indices k of the material boundaries of s (materials[k] != materials[k - 1], cyclic).
std::vector<std::size_t> materialBoundaries(const LineSection& s);

/// Angle theta(x) of Plan eq:cad1d_N (x in [0, dx), physical coordinates): linear between consecutive
/// material boundaries, the difference reduced to (-pi/2, pi/2]; 0 (N = x) without boundaries.
double normalAngle(const LineSection& s, double dx, double x);

/// N = (cos theta, sin theta, 0) of Plan eq:cad1d_N at x.
inline Vec3D lineNormal(const LineSection& s, double dx, double x) {
  const double t = normalAngle(s, dx, x);
  return Vec3D(std::cos(t), std::sin(t), 0.0);
}

}  // namespace geometry

/// Coefficients of N_i N_j (and N_i, unsigned: hasSigned = false) of Plan eq:cad1d_N, |p| <= 2Nx;
/// piecewise-smooth Gauss-Legendre quadrature between the boundaries. xmap: transformed coordinates X
/// (Plan section "M17": physical normal at x = f(X)); null: identity.
NormalFieldCoefficients normalField(const LineSection& s, double dx_um, const IndexMap& map,
                                    const CoordinateMap* xmap = nullptr);

/// Provider of the CAD slicer (Plan section "M23: Bandklassifikation"): line sections of profile bodies
/// in the band [y0, y0 + h] (absolute y of the bodies); never on the band ends (margin 1e-8 of the cell).
class ProfileLineProvider final : public LineSectionProvider {
 public:
  ProfileLineProvider(std::shared_ptr<const std::vector<geometry::ProfileBody>> bodies, double y0, double h, double dx,
                      std::string background, geometry::LineSectionOptions opt);
  LineSection section(double y_um) const override;
  double thickness_um() const override { return h_; }

 private:
  std::shared_ptr<const std::vector<geometry::ProfileBody>> bodies_;
  double y0_, h_, dx_;
  std::string bg_;
  geometry::LineSectionOptions opt_;
  double margin_;
};

}  // namespace fg
