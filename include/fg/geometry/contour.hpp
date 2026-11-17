// FourierGrating - curved cross-section boundaries (Plan section "M13: Konturdarstellung").
// Rational Bezier segments in the (x, z) plane, closed contours, regions with holes.
#pragma once

#include <array>
#include <string>
#include <vector>

namespace fg::geometry {

using Point2 = std::array<double, 2>;  ///< (x, z) in micrometres

/// Rational Bezier curve of degree p <= 9: control points P_j, weights w_j > 0, t in [0, 1].
struct BezierSegment {
  std::vector<Point2> P;
  std::vector<double> w;  ///< empty: all weights 1 (polynomial)

  int degree() const { return static_cast<int>(P.size()) - 1; }
  double weight(std::size_t j) const { return w.empty() ? 1.0 : w[j]; }
  Point2 point(double t) const;
  /// Point and first derivative d(x, z)/dt at t.
  void eval(double t, Point2& r, Point2& dr) const;
  /// de Casteljau split at t in the homogeneous space (rational segments stay exact).
  void split(double t, BezierSegment& left, BezierSegment& right) const;
  /// Length of the control polygon (upper bound of the arc length).
  double controlLength() const;
  BezierSegment reversed() const;
};

/// Closed chain of segments (end of one = start of the next, tolerance 1e-12 of the extent).
/// Counter-clockwise in (x, z): material on the left (Plan eq:polygon).
struct Contour2D {
  std::vector<BezierSegment> segments;

  /// Signed enclosed area (> 0 counter-clockwise), exact up to quadrature (Plan eq:contour0).
  double signedArea() const;
  Contour2D reversed() const;
  /// Throws if the chain is not closed or a segment is invalid.
  void validate() const;
  /// Bounding box of the control points: {xmin, zmin, xmax, zmax}.
  std::array<double, 4> bounds() const;
  /// Points on the contour with chord deviation (sagitta) <= tol, first point not repeated.
  std::vector<Point2> flatten(double tol) const;

  static Contour2D polygon(const std::vector<Point2>& vertices);
  /// Exact circle: four rational quadratic segments (weights 1, sqrt(2)/2, 1).
  static Contour2D circle(Point2 center, double radius);
  /// Ellipse with semi-axes (a, b), axis a rotated by angle (radians) against x.
  static Contour2D ellipse(Point2 center, double a, double b, double angle);
  /// Clamped NURBS curve (closed: first = last point) decomposed into Bezier pieces by knot
  /// insertion (Piegl-Tiller A5.6). weights empty: all 1.
  static Contour2D nurbs(int degree, const std::vector<double>& knots, const std::vector<Point2>& points,
                         const std::vector<double>& weights);
};

/// Do the segments ab and cd cross properly, each end farther than tol from the other segment's line? Ends that
/// lie on the other segment within rounding (touching bodies, results of Boolean operations) do not count (v2.24).
bool properCrossing(const Point2& a, const Point2& b, const Point2& c, const Point2& d, double tol);

/// Do the closed polylines p and q share a boundary piece longer than tol (segments collinear within tol)?
/// A region sharing such a piece with its own periodic image spans the period: its periodic union wraps around the
/// cell and the shared piece is no material interface (v2.21, v2.24).
bool sharedBoundary(const std::vector<Point2>& p, const std::vector<Point2>& q, double tol);

/// Material region: outer contour (counter-clockwise) with holes (clockwise) after normalize().
struct Region {
  Contour2D outer;
  std::vector<Contour2D> holes;
  std::string material;

  /// Orients outer counter-clockwise and holes clockwise, validates, prepares contains().
  void normalize(double flattenTol);
  /// Point inside the region (not periodic; crossing test on the flattened contours).
  bool contains(double x, double z) const;
  /// Overlap probes (v2.24): a point eps inside the region next to the parameter midpoint of every boundary segment
  /// (exact normal; material on the left after normalize()). Together with the test for proper boundary
  /// crossings this finds overlaps whose boundaries run along each other (collinear edges), where crossings and
  /// boundary points alone fail; touching regions (common boundary, no common area) contain no probe of each other.
  std::vector<Point2> interiorProbes(double eps) const;
  std::array<double, 4> bounds() const { return outer.bounds(); }
  double area() const;

 private:
  std::vector<std::vector<Point2>> poly_;  ///< flattened outer + holes
};

}  // namespace fg::geometry
