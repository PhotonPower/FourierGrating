// FourierGrating - open periodic NURBS profile of a 1D grating (Plan section "Erweiterung ab Version 2.37:
// NURBS-Profil in 1D-Gittern (M24)"): Bezier pieces (Piegl-Tiller A5.6), critical heights (Plan eq:nurbs1d_crit),
// checks (NURBS_*), line sections with the material rule of Plan eq:nurbs1d_rule. Independent of Nurbforge.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "fg/geometry/contour.hpp"
#include "fg/geometry/line_section.hpp"

namespace fg {
namespace geometry {

/// Clamped NURBS curve in the (x, y) plane (micrometres, y from the layer bottom); weights empty: all 1.
struct NurbsCurve {
  int degree = 3;
  std::vector<double> knots;
  std::vector<Point2> points;
  std::vector<double> weights;
};

/// Parameters t in (0, 1) where dy/dt = 0 inside a rational Bezier segment (roots of the numerator of Plan
/// eq:nurbs1d_crit, Bernstein subdivision and safeguarded Newton), ascending. A segment with constant y has none.
std::vector<double> yExtremaParameters(const BezierSegment& s);

/// A validated profile line (Plan section "M24: Konfiguration", "Prüfungen"): the Bezier pieces of one period,
/// running net in +x from (x0, y0) to (x0 + dx, y0), without self-intersection (also with its images at +-dx).
/// Errors start with the code ([NURBS_INPUT], [NURBS_NOT_PERIODIC], [NURBS_PERIOD], [NURBS_SELF_INTERSECTION]).
class NurbsProfile {
 public:
  NurbsProfile(const NurbsCurve& c, double dx_um, std::string below, std::string above);

  const std::vector<BezierSegment>& pieces() const { return pieces_; }
  double period() const { return dx_; }
  double yMin() const { return ymin_; }
  double yMax() const { return ymax_; }
  const std::string& below() const { return below_; }
  const std::string& above() const { return above_; }
  /// Critical heights (Plan section "M24: Kritische Höhen und Bänder"): interior y-extrema of every piece and
  /// the heights of all piece ends, ascending, merged within tol (absolute, micrometres).
  std::vector<double> criticalHeights(double tol) const;
  /// Line section at height y (never a critical height): crossings per piece on [0, 1), material right of a
  /// crossing from the sign of dy/dt (Plan eq:nurbs1d_rule), folded into [0, dx); the materials must alternate,
  /// else [NURBS_SELF_INTERSECTION].
  LineSection section(double y) const;

 private:
  std::vector<BezierSegment> pieces_;
  double dx_ = 0, ymin_ = 0, ymax_ = 0;
  std::string below_, above_;
};

}  // namespace geometry

/// Provider of profile1d "nurbs" (Plan section "M24: Linienschnitt"): sections of the layer [0, depth].
class NurbsProfileLineProvider final : public LineSectionProvider {
 public:
  NurbsProfileLineProvider(std::shared_ptr<const geometry::NurbsProfile> profile, double depth_um)
      : p_(std::move(profile)), h_(depth_um) {}
  LineSection section(double y_um) const override { return p_->section(y_um); }
  double thickness_um() const override { return h_; }
  /// Interior band boundaries 0 < eta_1 < ... < eta_B < depth (critical heights merged within 1e-12 depth).
  std::vector<double> bandBreaks() const;
  const geometry::NurbsProfile& profile() const { return *p_; }

 private:
  std::shared_ptr<const geometry::NurbsProfile> p_;
  double h_;
};

}  // namespace fg
