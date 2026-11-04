// FourierGrating - common solver interface and results.
// Both solvers (RK, RCWA) produce slice S-matrices in the reference basis; the stack assembly,
// global solution, efficiencies and near fields are shared (Plan sections 7, 10).
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "fg/core/index_map.hpp"
#include "fg/field/nearfield.hpp"
#include "fg/geometry/stack.hpp"
#include "fg/physics/incidence.hpp"
#include "fg/smatrix/smatrix.hpp"

namespace fg {

enum class SolverKind { RK, RCWA };

struct SolverOptions {
  SolverKind kind = SolverKind::RCWA;
  int Nx = 10, Nz = 0;
  std::string rkTableau = "gauss2";   ///< gauss2 | gauss3 | radau3
  double rkTol = 1e-8;                ///< step-doubling tolerance (Plan eq:stepcrit)
  int rcwaSlices = 40;                ///< staircase slices per profile layer
  double growthLimit = 1e8;           ///< slice criterion (Plan eq:slicecrit)
  bool storeInternal = false;         ///< keep S_L / S_U and slice records for near fields
};

struct OrderEfficiency {
  int n = 0, m = 0;
  double efficiency = 0;              ///< Plan eq:eff
  cplx amp_s = 0, amp_p = 0;          ///< s/p amplitudes of the order (isotropic media)
  double eff_s = 0, eff_p = 0;        ///< s/p parts of the efficiency (isotropic media)
};

struct Diagnostics {
  double energyBalance = 0;           ///< sum R + sum T (lossless media: 1)
  double maxConditionT11 = 0;
  long rkSteps = 0;
  std::vector<std::string> warnings;
};

struct Result {
  std::vector<OrderEfficiency> reflected, transmitted;
  double absorption = 0;
  Diagnostics diag;
};

struct Wavelength { double lambda_um; double phi, psi, delta; };  ///< angles in radians

/// Fourier coefficients of the field at height y (Plan section "Nahfeld"): F = ([Ex],[Ez],[Hx],[Hz])
/// plus [Ey], [Hy]; region -1 substrate, 0 layers, +1 superstrate.
struct FieldSample {
  double y_um = 0;
  int region = 0;
  VecC F, Ey, Hy;
};
struct FieldSampleSet {
  double k0 = 0;
  VecC alpha, gamma;   ///< normalized wave numbers of the orders
  std::vector<FieldSample> samples;
};

class Solver {
 public:
  virtual ~Solver() = default;
  /// Full diffraction computation for one wavelength / incidence.
  virtual Result solve(const Stack& stack, const Wavelength& wl) = 0;
  /// Field coefficients at the given heights (micrometres, y = 0 at the substrate interface).
  virtual FieldSampleSet fieldCoefficients(const Stack& stack, const Wavelength& wl,
                                           const std::vector<double>& y_um) = 0;
  /// Horizontal (y = const) and vertical (z = const) cuts in real space.
  virtual NearFieldResult nearField(const Stack& stack, const Wavelength& wl, const NearFieldRequest& req) = 0;
};

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt);

}  // namespace fg
