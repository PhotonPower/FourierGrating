// FourierGrating - separable coordinate transformations x = f(X), z = g(Z) (Plan section "M17:
// Koordinatentransformationen", eq:ct_map) and the adaptive spatial resolution of Granet (eq:asr).
#pragma once

#include <memory>
#include <vector>

#include "fg/core/types.hpp"

namespace fg {

struct Stack;

/// Transformation of one axis (JSON "coordinate_transform": {"x": {...}, "z": {...}}).
struct AxisTransformSpec {
  enum class Type { None, ASR };
  Type type = Type::None;
  double eta = 0;                       ///< ASR strength, 0 <= eta < 1 (0: no transformation)
  std::vector<double> breakpoints_um;   ///< empty: "auto" (all material jumps along the axis)
};

struct CoordinateTransformSpec {
  AxisTransformSpec x, z;
  /// Tests only (C++): the transformed path (matrices A, G, numerical half-space modes) also for
  /// eta = 0, where it must reproduce the diagonal path. Off: eta = 0 is no transformation.
  bool forceMatrixPath = false;
};

/// x = f(X) on one axis with f(X + d) = f(X) + d (Plan eq:ct_map); ASR mapping of Granet (Plan eq:asr):
/// on [x_{j-1}, x_j], f(X) = X - eta Delta_j / (2 pi) sin(2 pi (X - x_{j-1}) / Delta_j).
class CoordinateMap {
 public:
  /// breakpoints: any order, any period images (reduced modulo the period); at least one.
  CoordinateMap(double period, std::vector<double> breakpoints, double eta);

  double period() const { return d_; }
  double eta() const { return eta_; }
  /// Breakpoints in [0, period), sorted (fixed points f(x_j) = x_j).
  const std::vector<double>& breakpoints() const { return br_; }

  double f(double X) const;
  double df(double X) const;   ///< f'(X)
  double inverse(double x) const;  ///< X = f^{-1}(x) (Newton with bisection, monotone for eta < 1)

  /// Quadrature over one period in X (nodes, weights including 1/period): composite Gauss-Legendre
  /// between the breakpoints and the extra X breaks, panels short enough for the frequency
  /// pmax * 2 pi / period and the analyticity strip of 1/f'. Exact up to rounding for piecewise smooth
  /// integrands whose kinks lie at those breaks.
  struct Nodes {
    std::vector<double> X, w;
  };
  Nodes nodes(int pmax, const std::vector<double>& extraBreaks = {}) const;

  /// Coefficients [[1/f']]_p, p = -pmax..pmax (vector index p + pmax), Plan eq:ct_deriv.
  std::vector<cplx> inverseDerivativeCoefficients(int pmax) const;
  /// C_{nm} = 1/d int e^{i alpha_m X} e^{-i alpha_n f(X)} f'(X) dX, |n|, |m| <= N (row/column n + N),
  /// alpha_n = alpha0 + 2 pi n / d in 1/um (Plan eq:ct_rayleigh): coefficients in X -> physical Fourier
  /// coefficients in x.
  MatC toPhysical(double alpha0, int N) const;
  /// D_{mn} = 1/d int e^{i alpha_n f(X)} e^{-i alpha_m X} dX (Plan eq:ct_inc): physical orders -> X.
  MatC fromPhysical(double alpha0, int N) const;

 private:
  /// Interval j of X (reduced into [br_0, br_0 + d)) and the period shift.
  void locate(double X, std::size_t& j, double& shift, double& u) const;
  double d_, eta_;
  std::vector<double> br_;
};

/// Transformation resolved for a stack: null = identity on that axis.
struct ResolvedTransform {
  std::shared_ptr<const CoordinateMap> x, z;
  bool force = false;  ///< CoordinateTransformSpec::forceMatrixPath
  bool active() const { return x || z || force; }
};

/// Resolves stack.transform (Plan section "ASR-Abbildung"): eta = 0 gives the identity (null),
/// "auto" breakpoints are the material jumps of the lamellar layers and the vertex coordinates of
/// polygons (polygon2d, polygonal contour2d regions) along the axis.
ResolvedTransform resolveTransform(const Stack& stack);

/// Maps a 1D piecewise constant function given on physical intervals 0 = x_0 < ... < x_K = d
/// (values[k] on [x_k, x_k+1]) to X (Plan section "Permittivitaet und Normalenfeld": X_j = f^{-1}(x_j)),
/// again as intervals 0 = X_0 < ... = d; the map's own breakpoints are added as extra breaks (same
/// value) so that smooth factors evaluated at f(X) are integrated piecewise. No-op for map == nullptr.
template <class T>
void mapIntervals(const CoordinateMap* map, std::vector<double>& breaks, std::vector<T>& values);

namespace detail {
/// Interval starts in [0, d) and the index of the physical interval each one belongs to.
void mappedStarts(const CoordinateMap& map, const std::vector<double>& breaks, std::vector<double>& starts,
                  std::vector<std::size_t>& index);
}  // namespace detail

template <class T>
void mapIntervals(const CoordinateMap* map, std::vector<double>& breaks, std::vector<T>& values) {
  if (!map) return;
  std::vector<double> starts;
  std::vector<std::size_t> index;
  detail::mappedStarts(*map, breaks, starts, index);
  std::vector<T> v;
  for (std::size_t k : index) v.push_back(values[k]);
  starts.push_back(map->period());
  breaks = std::move(starts);
  values = std::move(v);
}

}  // namespace fg
