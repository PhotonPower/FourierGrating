// FourierGrating - layer models (Plan section "Schichtmodelle und Scheibenzerlegung").
// All lengths here are in micrometres; conversion to normalized units happens in the solver setup.
#pragma once

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "fg/geometry/contour.hpp"
#include "fg/geometry/material.hpp"
#include "fg/geometry/section.hpp"

namespace fg {

/// FFF_iso: book (IV.32), Q = [[eps]] - Delta [[NN]]; FFF_sym: reciprocity-preserving
/// Q = [[eps]] - (Delta [[NN]] + [[NN]] Delta)/2 (Plan eq:Qsym, default for isotropic profiles).
enum class Factorization { Auto, Laurent, LiLamellar, FFF_iso, FFF_sym, FFF_rig, FFF_A, FFF_B };

struct HomogeneousLayer {
  double thickness_um = 0;
  std::string material;
};

struct Lamellar1DLayer {            ///< y-invariant, piecewise constant in x
  double thickness_um = 0;
  std::vector<double> breakpoints_um;  ///< 0 = x0 < x1 < ... < xK = dx
  std::vector<std::string> materials;  ///< K entries
};

struct Profile1DLayer {             ///< region y < g(x) is 'below', y > g(x) is 'above'
  enum class Shape { Sinus, PiecewiseLinear, Function };
  double depth_um = 0;                 ///< a = max g - min g; g in [0, a]
  Shape shape = Shape::Sinus;
  double x0_um = 0;                    ///< Sinus: g = a/2 (1 + cos(K (x - x0)))
  std::vector<double> xs_um, ys_um;    ///< PiecewiseLinear: periodic polygon, xs ascending in [0, dx)
  std::function<double(double)> g;     ///< Function: g(x) in [0, a], period dx
  std::function<double(double)> dg;    ///< Function: g'(x) (optional, else finite differences)
  std::string below, above;
  int slices = 0;                      ///< RCWA staircase slices (0 = solver default)
};

/// Settings of Algorithm NF (Plan alg:nfield) for y-invariant 2D cross-sections.
/// Tensor (default, Plan v1.5): the projector field n n^T of the nearest-boundary direction is
/// smoothed (sign-free, exact for z-invariant lamellae). Vector: the signed-distance gradient is
/// smoothed and normalized (original Plan v1.4 formulation, binary layouts only).
struct NormalFieldOptions {
  enum class Method { Tensor, Vector };
  Method method = Method::Tensor;
  int grid = 512;      ///< raster points per period (at least 8 (2N+1) are used)
  double sigma = 3.0;  ///< Gaussian width in raster points (along the shorter period)
  double tau = 1e-3;   ///< Vector method: fallback threshold relative to max |g|
};

struct Polygon2DLayer {             ///< y-invariant cross-section in the (x, z) cell
  double thickness_um = 0;
  std::string background;
  /// Simple polygon (x, z) in micrometres; may cross the cell boundary (periodic images must
  /// not overlap). Orientation is normalized to counter-clockwise internally.
  struct Inclusion { std::vector<std::array<double, 2>> vertices; std::string material; };
  std::vector<Inclusion> inclusions;   ///< non-overlapping
  NormalFieldOptions normalField;
};

struct Profile2DLayer {             ///< region y < g(x, z) is 'below', y > g(x, z) is 'above'
  /// SinusX: g = a/2 (1 + cos Kx (x - x0)) (z-invariant, V9);
  /// EggCrate: g = a/4 (1 + cos Kx (x - x0)) (1 + cos Kz (z - z0));
  /// Function: user g(x, z) in [0, a] (periodic), optional gradient.
  enum class Shape { SinusX, EggCrate, Function };
  double depth_um = 0;                 ///< a; g in [0, a]
  Shape shape = Shape::Function;
  double x0_um = 0, z0_um = 0;
  std::function<double(double, double)> g;
  std::function<std::array<double, 2>(double, double)> grad;  ///< (dg/dx, dg/dz), optional
  std::string below, above;
  int slices = 0;                      ///< RCWA staircase slices (0 = solver default)
  int zSamples = 0;                    ///< trapezoid nodes in z for eps(y) (0 = default)
};

/// y-invariant cross-section with curved boundaries (Plan section "M13: Neuer Schichttyp Contour2D"):
/// background and pairwise disjoint regions (circles, ellipses, polygons, NURBS, with holes);
/// exact coefficients from contour integrals. Factorizations as Polygon2D.
struct Contour2DLayer {
  double thickness_um = 0;
  std::string background;
  std::vector<geometry::Region> regions;
  NormalFieldOptions normalField;
};

/// Layer whose cross-section changes with y and cannot be written as a profile y = g(x, z) (Plan
/// section "M13: Hoehenabhaengige Schichtmodelle"; CAD slicer M14). The provider delivers regions and
/// boundary samples per height; coefficients, normal fields and Q(y) are computed on demand and cached
/// (key: y bit-exact). RK evaluates at the stage nodes, RCWA at the slice centres. Factorizations as
/// Contour2D. Not available in JSON (C++ API; M14: layer type "cad").
struct SectionedLayer {
  std::shared_ptr<const SectionProvider> provider;
  std::string background;
  int slices = 0;  ///< RCWA staircase slices (0 = solver default)
  NormalFieldOptions normalField;
  /// > 0: Q(y) from barycentric interpolation of the Q blocks at this many Chebyshev nodes (first kind,
  /// doubled until K and 2K agree to 1e-10 max|Q|, at most 256; then exact sections with a warning);
  /// 0: exact sections at every height (Plan M14 "Kosten und Interpolation in y").
  int chebyshevNodes = 0;
};

using LayerSpec = std::variant<HomogeneousLayer, Lamellar1DLayer, Profile1DLayer, Polygon2DLayer,
                               Profile2DLayer, Contour2DLayer, SectionedLayer>;

struct LayerEntry {
  LayerSpec spec;
  Factorization factorization = Factorization::Auto;
};

}  // namespace fg
