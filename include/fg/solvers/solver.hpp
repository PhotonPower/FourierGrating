// FourierGrating - common solver interface and results.
// Both solvers (RK, RCWA) produce slice S-matrices in the reference basis; the stack assembly,
// global solution, efficiencies and near fields are shared (Plan sections 7, 10).
#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fg/core/error.hpp"
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
  /// slice criterion (Plan eq:slicecrit); 0 = auto: eps^-1/2 of the working precision (FP64: 1e8 as
  /// before, FP32/mixed: 4e3, Plan section "M16: Angepasste Parameter")
  double growthLimit = 0;
  bool storeInternal = false;         ///< keep S_L / S_U and slice records for near fields
  std::string backend = "auto";       ///< "auto" | "cpu" | "cuda" (Plan section "CUDA-Backend", M10)
  int backendThreshold = 160;         ///< auto: matrix size from which the CUDA backend is used
  /// 1D in-plane isotropic problems: TE and TM as separate 2 n_f systems in solve() (Plan
  /// eq:TE, eq:TM, section "Spezialfaelle", v2.2). Near fields always use the full 4 n_f system.
  /// Off by default: results then bit-identical to v0.10 (AGENTS.md section 7).
  bool splitPolarizations = false;
  /// Sampled Fourier coefficients via fourier::fft for any sample count (Plan "Allgemeine
  /// Funktionen", v2.2); off: v0.10 paths (bit-identical).
  bool fftSampling = false;
  /// Thread model (Plan section "M15", step 2), 0 = auto: threads over scan points, over the
  /// layers/slices of one computation, and inside Eigen (blas). Their product never exceeds
  /// availableThreads() (FG_NUM_THREADS). Results do not depend on scan/slice threads.
  int threadsScan = 0, threadsSlices = 0, threadsBlas = 1;
  /// Star-product chain: false = sequential (v0.10 order, bit-identical), true = balanced tree of
  /// fixed shape (parallel; bit-identical for any thread count, <= 1e-12 from sequential).
  bool treeReduction = false;
  /// Scattered-field formulation (Plan section "Streufeldformulierung"): E = E_bg + E_s with the
  /// background stack (structured layers replaced by their background material) solved analytically and
  /// the source sigma = (M - M_bg) F_bg in every structured slice. Automatic for aperiodic structures
  /// under oblique incidence; this switch forces it (periodic: efficiencies of E_bg + E_s).
  bool scatteredField = false;
  /// Working precision of the O(n^3) operations (Plan section "M16"): "auto" | "fp64" | "mixed" |
  /// "fp32". auto: FP64 on the CPU; on the GPU by device policy (FP32 on a GeForce RTX 3090).
  /// FG_PRECISION overrides it.
  std::string precision = "auto";
  /// Repeat in FP64 if a lossless structure misses energy conservation by > 1e-3:
  /// -1 = auto (on for fp32), 0 = off, 1 = on (also for mixed).
  int precisionFallback = -1;
  /// Gradients (Plan section "M19: Schnittstellen", Speicher): keep only every k-th prefix S-matrix and no
  /// slice data in the forward pass; the backward pass recomputes them per segment (same result).
  int checkpointEvery = 1;
};

/// Resolved thread counts of a run with `scanPoints` points (Plan section "M15", step 2): auto
/// gives the scan level min(points, cores / blas) threads and the slice level the rest.
struct ThreadPlan {
  int scan = 1, slices = 1, blas = 1;
};
ThreadPlan threadPlan(const SolverOptions& opt, std::size_t scanPoints, int cores);

struct OrderEfficiency {
  int n = 0, m = 0;
  double efficiency = 0;              ///< Plan eq:eff
  cplx amp_s = 0, amp_p = 0;          ///< s/p amplitudes of the order (isotropic media)
  double eff_s = 0, eff_p = 0;        ///< s/p parts of the efficiency (isotropic media)
};

/// n and eps of a dispersive material at the scan wavelength (M22, Plan "Ausgaben"); component 0 for
/// isotropic materials, 1..3 for the principal values.
struct MaterialValue {
  std::string name;
  int component = 0;
  cplx n = 1.0, eps = 1.0;
};

struct Diagnostics {
  double energyBalance = 0;           ///< sum R + sum T (lossless media: 1)
  double maxConditionT11 = 0;
  long rkSteps = 0;
  std::string precision = "fp64";     ///< working precision actually used (M16)
  bool precisionFallback = false;     ///< repeated in FP64 (M16 safeguard)
  long refinementSteps = 0;           ///< iterative-refinement steps (mixed)
  long refinementFallbacks = 0;       ///< mixed solves that fell back to an FP64 LU
  std::vector<std::string> warnings;
  std::vector<MaterialValue> materials;  ///< dispersive materials at the wavelength (M22)
};

/// Aperiodic structures (Plan section "M18: Ausgaben"): power of the scattered field E_s = E - E_bg
/// through the planes y = top and y = 0 over the physical region, normalized to the incident intensity
/// (scattering width in um for 1D, cross section in um^2 for 2D); absorbed: net total-field flux into
/// the stack minus that of the background stack.
struct ScatteringResult {
  bool present = false;
  double up = 0, down = 0, absorbed = 0;
  /// Scattered power through the lateral sides of the box (window x [0, top]); up + down + lateral is the
  /// total scattered power, independent of the window (the split into up and down is not).
  double lateral = 0;
  double total() const { return up + down + lateral; }
};

/// Modal excitation (Plan "Einfall (b)", M18): coefficients of the guided modes of the structured
/// half spaces; power relative to the incident mode (flux over the physical region).
struct ModeCoefficient {
  int index = 0;            ///< guided-mode number (descending n_eff)
  cplx nEff = 0;            ///< normalized propagation constant along y
  cplx amplitude = 0;       ///< power-normalized amplitude (phase of the numerical eigenvector)
  double power = 0;         ///< |r|^2 resp. |t|^2
  double teFraction = 0;    ///< share of E_z in |E|^2 (1: TE, 0: TM in 1D)
};
struct ModalResult {
  bool present = false;
  ModeCoefficient incident;
  std::vector<ModeCoefficient> reflected, transmitted;
  double radiated = 0;      ///< 1 - sum of the guided powers (radiation into the half spaces and the zones)
};

/// Rotationally symmetric structures (M20, Plan section "M20: Ausgaben"): powers relative to the incident
/// beam, from the flux through y = top and y = 0 (Plan eq:cyl_flux); A = 1 - R - T also contains the power
/// that guided modes carry away laterally.
/// Angular spectrum of one harmonic on one side (Plan "M20: Ausgaben"): at the propagating nodes of the spectral
/// grid (theta = asin(k / n)), dP/dtheta per radian split into s (TE, (e_+ + e_-)/2) and p (TM, (e_+ - e_-)/2);
/// beam: relative to the incident power, plane wave: um^2 per radian. sum weight (s + p) is the power of the
/// harmonic on that side (C_Rm, C_Tm or C_Supm, C_Sdownm).
struct CylAngularSpectrum {
  int m = 0;
  int side = 1;  ///< +1 up into the superstrate, -1 down into the substrate
  std::vector<double> theta_deg, weight, s, p;  ///< weight: quadrature weight in theta (rad)
};

/// Far field of one side: dP/dOmega per steradian at (theta nodes) x (phi_k), s and p, row-major (theta outer);
/// beam: relative to the incident power, plane wave: differential cross section in um^2 / sr.
struct CylFarField {
  int side = 1;
  std::vector<double> theta_deg, phi_deg, s, p;
};

struct CylindricalResult {
  bool present = false;
  bool planeWave = false;  ///< plane wave (scattered field): cross sections instead of R, T, A
  double reflected = 0, transmitted = 0, absorbed = 0;  ///< beam: relative to the incident power
  /// Plane wave (Plan "M20: Ausgaben"): scattering cross sections in um^2, up into the superstrate,
  /// down into the substrate.
  double scatteredUp_um2 = 0, scatteredDown_um2 = 0, scattering_um2 = 0;
  std::vector<int> harmonics;
  /// Per harmonic: beam relative to the total incident power; plane wave the cross sections up/down (um^2).
  std::vector<double> reflectedPerHarmonic, transmittedPerHarmonic;
  int nodes = 0;  ///< N_k of the spectral grid
  /// M21 (Plan sec:cylP): order P of the rotational symmetry (0: axial, one harmonic at a time) and the truncation L
  /// of the classes m0 + l P; harmonics lists every m of the computed classes.
  int rotation = 0, classL = 0;
  /// Volume integrals over the rings (Plan "M20: Ausgaben", Gauss nodes in r and y, summed over the
  /// harmonics). Plane wave: absorption cross section sigma_abs = 1/(2 I_inc) int Im eps |E|^2 dV,
  /// extinction sigma_ext = sigma_sca + sigma_abs and, as a control, the optical theorem in volume form
  /// sigma_ext = 1/(2 I_inc) Im int (eps - eps_inf) E . conj(E_bg) dV (the forward-scattering amplitude);
  /// evaluated with lossy rings or cylindrical.optical_theorem, NaN otherwise and when the background
  /// (superstrate, substrate, homogeneous layers, outer ring material) is lossy.
  double absorption_um2 = 0, extinction_um2 = 0, extinctionTheorem_um2 = 0;
  /// Beam: power absorbed inside the rings relative to the incident power (A minus absorbedRings is the
  /// lateral outflow plus the absorption outside the rings).
  double absorbedRings = 0;
  /// Per harmonic and side (up, and down when the substrate is lossless); empty for harmonics without power.
  std::vector<CylAngularSpectrum> angularSpectrum;
  /// Up and down (lossless substrate) with cylindrical.far_field; empty otherwise.
  std::vector<CylFarField> farField;
};

struct Result {
  std::vector<OrderEfficiency> reflected, transmitted;
  double absorption = 0;
  CylindricalResult cylindrical;  ///< "cylindrical" stacks only (reflected/transmitted then empty)
  ScatteringResult scattering;  ///< aperiodic structures only (reflected/transmitted then empty)
  ModalResult modal;            ///< modal excitation only
  Diagnostics diag;
};

struct Wavelength { double lambda_um; double phi, psi, delta; };  ///< angles in radians

/// Fourier coefficients of the field at height y (Plan section "Nahfeld"): F = ([Ex],[Ez],[Hx],[Hz])
/// plus [Ey], [Hy]; region -1 substrate, 0 layers, +1 superstrate.
struct FieldSample {
  double y_um = 0;
  int region = 0;
  VecC F, Ey, Hy;
  /// Scattered-field formulation (M18): the background field (Ex, Ey, Ez, Hx, Hy, Hz) of the incidence
  /// order, ~ exp(i k0 (alpha0 x + gamma0 z)); empty otherwise.
  VecC bg;
};
struct FieldSampleSet {
  double k0 = 0;
  VecC alpha, gamma;   ///< normalized wave numbers of the orders
  double alpha0 = 0, gamma0 = 0;  ///< normalized wave numbers of the incident wave
  /// Scattered-field formulation: F contains the background (periodic, normal incidence); false for the
  /// Bloch-free scattered field of oblique incidence with PML (F = E_s, the background is in bg).
  bool backgroundIncluded = false;
  std::vector<FieldSample> samples;
};

/// S-matrix of the whole stack in the s/p basis of the plane-wave orders (Plan eq:Ssp), for one
/// wavelength and incidence direction. Column r is incidence in order r (another direction of the
/// Bloch family), but the truncation is centred on the incident order: only the central column
/// equals a solve() of its own direction exactly, the others to the truncation error. Index r = order
/// orders[r], s component at r, p component at n_f + r (s, p as in the efficiency output, Plan
/// section "Leistungsfluss und Effizienzen"); amplitudes referenced at the top of the stack
/// (superstrate) and at y = 0 (substrate). R_top: from above back up (S12), T_down: from above into
/// the substrate (S22), T_up: from below into the superstrate (S11), R_bottom: from below back
/// down (S21), Plan eq:Sdef.
struct OrderSMatrix {
  std::vector<std::pair<int, int>> orders;  ///< (n, m) per index r
  VecC alpha, gamma;                        ///< normalized tangential wave numbers (k0 = 1)
  VecC betaSup, betaSub;                    ///< normalized normal wave numbers (Im > 0: evanescent)
  MatC R_top, T_down, T_up, R_bottom;       ///< 2 n_f x 2 n_f each
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
  /// S-matrix in the s/p basis of the orders (Plan eq:Ssp): periodic stacks with homogeneous isotropic
  /// superstrate and substrate, without transformed coordinates.
  virtual OrderSMatrix orderSMatrix(const Stack& stack, const Wavelength& wl) {
    (void)stack;
    (void)wl;
    throw NotImplemented("orderSMatrix for this solver");
  }
};

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt);

/// Wavelength/angle scan, parallel over the scan points (OpenMP threads, one solver and - with the
/// CUDA backend - one stream per thread: batching over wavelengths, M10). Results in scan order,
/// identical to sequential solve() calls.
std::vector<Result> solveScan(const SolverOptions& opt, const Stack& stack, const std::vector<Wavelength>& scan);

}  // namespace fg
