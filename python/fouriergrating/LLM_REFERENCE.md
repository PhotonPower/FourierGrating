# fouriergrating: compact reference for agents

Get this text at runtime: `python -m fouriergrating` or `fouriergrating.llm_reference()`.
Every `json` block below is a complete configuration and every `python` block runs as is (both are tested).

Fourier modal method for periodic gratings (1D, 2D), isolated structures (PML) and axially symmetric structures.
Two solvers: `rcwa` (staircase, eigenmodes) and `rk` (differential method, implicit Runge-Kutta). Python wraps a
C++ core, and the core does all the physics. The input is one JSON-like `dict` (the same format as the CLI `fgsolve`).

## Conventions
- Units: lengths in µm, angles in degrees. Time dependence exp(-iωt): absorbing materials have Im n > 0 (Im ε > 0).
- Axes: x = grating period (`period_um.x`), z = second period (2D only), y = stacking direction pointing up.
  `layers` run from the bottom (substrate, y < 0) to the top; the superstrate is on top and must be lossless.
  Light comes from the superstrate.
- Incidence: `phi_deg` = polar angle from the y axis in the superstrate, `psi_deg` = azimuth of the plane of
  incidence measured from x (0 = classical mounting, otherwise conical), `delta_deg` = polarization
  (0 = s / TE, 90 = p / TM, anything else = linear mix).
- Orders (n, m): in-plane wave numbers α_n = α0 + n λ/dx and γ_m = γ0 + m λ/dz; m = 0 in 1D. Efficiency = power of
  the order divided by the incident power. `R` = reflected into the superstrate, `T` = transmitted into the
  substrate. For lossless stacks, Σ R + Σ T = 1, which `diagnostics["energy_balance"]` reports.
- Truncation: orders |n| ≤ Nx, |m| ≤ Nz. Convergence must be checked by increasing N (and RCWA slices).

## Minimal run
```json
{"dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 15},
 "period_um": {"x": 1.0},
 "incidence": {"lambda_um": 0.6328, "phi_deg": 20, "psi_deg": 0, "delta_deg": 90},
 "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "bar": {"n": 2.0}},
 "superstrate": "air", "substrate": "glass",
 "layers": [{"type": "lamellar1d", "thickness_um": 0.4,
             "breakpoints_um": [0.0, 0.45, 1.0], "materials": ["bar", "air"]}]}
```
```python
import fouriergrating as fg
cfg = {"dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 15}, "period_um": {"x": 1.0},
       "incidence": {"lambda_um": 0.6328, "phi_deg": 20, "psi_deg": 0, "delta_deg": 90},
       "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "bar": {"n": 2.0}},
       "superstrate": "air", "substrate": "glass",
       "layers": [{"type": "lamellar1d", "thickness_um": 0.4, "breakpoints_um": [0.0, 0.45, 1.0],
                   "materials": ["bar", "air"]}]}
res = fg.solve(cfg)                       # first scan point only
print(res.efficiency("R", -1), res.efficiency("T", 0), res.absorption)
print(res.reflected["n"], res.reflected["eff"])     # structured array, one row per propagating order
assert abs(res.diagnostics["energy_balance"] - 1) < 1e-8
```

## Configuration keys
Top level (`*` = required):
| key | value |
|---|---|
| `dimension`* | `"1D"` or `"2D"` (omitted only with `cylindrical`) |
| `solver` | `"rcwa"` or `"rk"` |
| `truncation`* | `{"Nx": int, "Nz": int}` (Nz only in 2D) |
| `period_um`* | `{"x": dx}` or `{"x": dx, "z": dz}` (2D: dz > 0) |
| `incidence`* | `lambda_um`, `phi_deg`, `psi_deg`, `delta_deg`: each a number, a list or `{"from": a, "to": b, "steps": k}`; the scan is the product λ × φ × ψ × δ |
| `materials`* | `{name: material}` (see below) |
| `superstrate`*, `substrate`* | material names (half spaces, homogeneous) |
| `layers`* | list from bottom to top (see below) |
| `sweep` | `{"parameter": "/json/pointer", "values": number, list or range}`: repeats the scan for each value |
| `rk` | `{"tol": 1e-8, "tableau": "gauss2"}` (`gauss2`, `gauss3`, `radau3`) |
| `rcwa` | `{"slices": 60}`: default staircase slices of y-variable layers |
| `nearfield` | `{"enabled": true, "horizontal_y_um": y, "vertical_z_um": z, "y_range_um": [y0, y1], "nx": 400, "ny": 300, "nz": 200, "periods": 2, "reconstruction": "fff" or "direct", "field": "total" or "scattered", "extent_um": 0}`; `extent_um` only with `cylindrical` (x, z in [-extent, extent], 0 = 1.5 R_s; the vertical cut is the plane z = `vertical_z_um`) |
| `output` | `{"dir": "out/x", "efficiencies": "efficiencies.csv", "materials": "materials.csv"}` (CLI files) |
| `coordinate_transform` | `{"x": {"type": "asr", "eta": 0.8, "breakpoints_um": "auto" or list}, "z": {...}}`: adaptive spatial resolution, for metals and TM polarization |
| `aperiodic` | `{"pml": {"x": {"width_um": w, "beta": 1.0}, "z": {...}}, "formulation": "auto" or "scattered", "excitation": "plane_wave" or {"mode": j} or {"n_eff": v}}`: an isolated structure in a cell with PML |
| `superstrate_structure`, `substrate_structure` | y-invariant section (`lamellar1d`, `polygon2d`, `contour2d` without thickness) for structured half spaces (PML) |
| `cylindrical` | axially symmetric structures about the y axis: `{"symmetry": "axial", "harmonics": "auto" or list, "spectrum": {"kappa_max": 6}, "factorization": "li_r" or "laurent", "optical_theorem": false}`; no `period_um`/`truncation`; layers `homogeneous` or `rings`; incidence `{"type": "gaussian_beam", "lambda_um", "waist_um", "focus_y_um", "polarization": "x"}` or `{"type": "plane_wave", ...angles}`. CSV sides: beam `C_R`, `C_T`, `C_A`; plane wave cross sections in µm² `C_Sup`, `C_Sdown`, `C_Ssca`, `C_Sabs`, `C_Sext` and the control `C_Sext_ot` (with `optical_theorem: true` or lossy rings; far from `C_Sext` = not converged) |
| `backend` | `"auto"`, `"cpu"` or `"cuda"`; `backend_threshold` (160) |
| `precision` | `"auto"`, `"fp64"`, `"mixed"` or `"fp32"`; `precision_fallback` |
| `threads` | `{"scan": "auto", "slices": "auto", "blas": 1}` |
| other | `split_polarizations`, `fft_sampling`, `reduction` (`"sequential"` or `"tree"`), `growth_limit`, `checkpoint_every`, `scattered_field` |

Materials (one key per entry). Each value is a number, `[re, im]` or a dispersive model:
| form | meaning |
|---|---|
| `{"n": 1.5}`, `{"n": [0.2, 3.4]}` | refractive index |
| `{"eps": [-10, 1]}` | permittivity |
| `{"eps_tensor": [[[re, im] x3] x3]}` | full 3×3 tensor (rows x, y, z) |
| `{"eps_principal": [e1, e2, e3], "euler_deg": [a, b, c]}` | principal values, rotated R diag Rᵀ (z-x-z) |
| `{"sellmeier": {"A": 1, "B": [...], "C_um2": [...], "range_um": [a, b]}}` | gives ε, as a value of `n`/`eps` |
| `{"cauchy": {"A": 1.45, "B_um2": 0.004, "C_um4": 0, "k": 0}}` | gives n |
| `{"drude_lorentz": {"eps_inf", "omega_p_eV", "f0", "gamma0_eV", "oscillators": [{"f", "omega_eV", "gamma_eV"}]}}` | gives ε |
| `{"table": {"data": [[λ, n, k], ...] or "file": "x.csv", "unit": "um", "nm" or "eV", "interpolation": "linear" or "cubic", "extrapolate": "error" or "constant"}}` | columns λ, n[, k] under `n`; λ, Re ε[, Im ε] under `eps` |

Layers (`type`, plus `factorization` where applicable):
| type | keys | dims |
|---|---|---|
| `homogeneous` | `thickness_um`, `material` | 1D, 2D |
| `lamellar1d` | `thickness_um`, `breakpoints_um` (must start at 0 and end at dx; this is not checked), `materials` (one per interval) | 1D, 2D (z-invariant) |
| `profile1d` | `depth_um`, `below`, `above`, `shape`: `sinus` (`x0_um`), `triangle` (`apex_x_um`) or `points` (`x_um`, `y_um`, periodic polygon); `slices` | 1D, 2D |
| `polygon2d` | `thickness_um`, `background`, `inclusions`: `[{"material", "vertices_um": [[x, z], ...]}` or `{"material", "rectangle_um": {"center": [x, z], "size": [wx, wz]}}]`, `normal_field` | 2D |
| `contour2d` | `thickness_um`, `background`, `regions`: `[{"material", "circle": {"center_um", "radius_um"}}`, `{"ellipse": {"center_um", "radii_um", "angle_deg"}}`, `{"polygon_um": [[x, z], ...]}` or `{"nurbs": {"degree", "knots", "points_um", "weights"}}`, with optional `"holes": [...]`]; regions must not overlap (images included) | 2D |
| `profile2d` | `depth_um`, `below`, `above`, `shape`: `egg_crate` or `sinus_x`, `x0_um`, `z0_um`, `slices`, `z_samples` | 2D |
| `cad` | STEP bodies, expanded into bands of the types above (needs `fg.cad.available()`); see "CAD" | 1D, 2D |
| `rings` | `thickness_um`, `radii_um`, `materials` (one more than radii; the last one is outside); with `cylindrical` | axial |

Factorizations: `auto` (recommended), `laurent`, `li_lamellar`, `FFF_iso`, `FFF_sym`, `FFF_rig` (profiles), `FFF_A`
and `FFF_B` (anisotropic). The default `auto` picks the right one per layer type: Li rule for lamellae,
`FFF_sym` / `FFF_rig` for profiles, `FFF_sym` / `FFF_B` for 2D sections.

## Python API
```text
fg.solve(cfg, profile=False) -> Result                 first scan point ("sweep" ignored)
fg.solve_scan(cfg, profile=False) -> ScanResult        all sweep values × scan points (parallel)
fg.near_field(cfg) -> {"horizontal": d, "vertical": d, "warnings": [...]} needs cfg["nearfield"]["enabled"]
fg.load_config(path) -> dict                           JSON file; relative output/table/STEP paths resolved
fg.material_eps(cfg, name, lambda_um) -> complex or array   (eps_principal: last axis of length 3)
fg.Simulation(cfg)                                     layer models built once, many solves
  .solve(lambda_um=None, phi_deg=None, psi_deg=None, delta_deg=None) -> Result   (None = first scan point)
  .set_parameter(pointer, value)                       e.g. "/layers/0/thickness_um"; rebuilds what changed
  .value_and_grad(objective, params, lambda_um=..., ...) -> (L, grad)   discrete adjoint; also sets .last_result
  .vjp(cotangents, params) -> ndarray
  .smatrix(lambda_um=None, phi_deg=None, psi_deg=None) -> OrderSMatrix
  .config (dict), .model_builds (int)
fg.set_backend("auto", threshold=160); fg.get_backend(); fg.cuda_available(); fg.core_variant(); fg.__version__
fg.FourierGratingError (RuntimeError): every configuration or computation error; NotImplementedError: not available yet
fg.ORDER_DTYPE                                         fields n, m, eff, eff_s, eff_p, amp_s, amp_p
fg.materials.from_refractiveindex(path_or_yaml_text) -> material dict   (needs PyYAML)
fg.cad.available(), fg.cad.check(step, cfg) -> dict, fg.cad.preview(step, cfg, y_um=None) -> svg, fg.cad.llm_tool_spec() -> dict
fouriergrating.torch.efficiencies(sim, params, values_tensor, orders) -> tensor (autograd)
fouriergrating.jax.efficiency_function(sim, params, orders) -> f(values) (custom_vjp)
fg.llm_reference() -> str                              this text
```
- `Result`: `reflected`, `transmitted` (ORDER_DTYPE arrays), `absorption`, `diagnostics`, `incidence` (scan
  points), `efficiency(side, n, m=0)`, plus `scattering` (PML: `up`, `down`, `lateral`, `absorbed`; widths in µm,
  in 2D cross sections in µm²), `modal` (mode excitation) and `cylindrical` (axial structures: dict with `reflected`,
  `transmitted`, `absorbed`, `absorbed_rings`, per-harmonic lists; plane wave `scattering_um2`, `absorption_um2`,
  `extinction_um2`, `extinction_theorem_um2`).
- `absorption` = 1 − Σ R − Σ T. `diagnostics`: `energy_balance` (Σ R + Σ T; 1 for lossless stacks, a convergence check), `warnings`, `rk_steps`,
  `max_condition_t11`, `precision`, `materials`, `timings`.
- `ScanResult`: `points` (list of Result, each with `.incidence`), `sweep_parameter`, `columns` (numpy arrays: `sweep_value`,
  `lambda_um`, `phi_deg`, `psi_deg`, `delta_deg`, `side`, `n`, `m`, `efficiency`, `eff_s`, `eff_p`, `amp_s`,
  `amp_p`), `to_pandas()`.
- `near_field` cut dict: `Ex`…`Hz` (complex, E normalized to |E_inc| = 1, H as Z0·H), `absE2`, `absH2`, `Sy`,
  `material`, `axis1_um`, `axis2_um`, `axis1_name`, `axis2_name`, `meta`. Vertical arrays have the shape (ny, nx).
- `OrderSMatrix`: `orders` ((n, m) per index), `alpha`, `gamma`, `beta_sup`, `beta_sub` (normalized; ×k0 = 2π/λ),
  `R_top`, `T_down`, `T_up`, `R_bottom` (s components at r, p components at n_f + r), `k0`, `incidence`,
  `index(n, m)`.
- Objectives and cotangents: `{("R", -1, 0): 1.0, ("T", 0, 0): -0.5}` or `[("R", -1, 0, 1.0)]`; PML:
  `{("scattering", "up"): w}`, `{("mode", "R", j): w}`.
- Gradient parameters are JSON pointers to numbers of the configuration: `/layers/k/thickness_um`,
  `/layers/k/breakpoints_um/j`, `/layers/k/depth_um`, `/layers/k/x0_um`, any number of `polygon2d`/`contour2d`
  geometry, `/materials/name/n` (`/n/0` = real part, `/n/1` = imaginary part), `/materials/name/eps`, model values such as
  `/materials/sio2/n/sellmeier/B/0`, and `/incidence/lambda_um`, `phi_deg`, `psi_deg`, `delta_deg`. The cost is
  about 1.3 to 2.5 forward solves for any number of parameters.

## Recipes
```python
import fouriergrating as fg
cfg = {"dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 10}, "period_um": {"x": 0.6},
       "incidence": {"lambda_um": {"from": 0.5, "to": 0.7, "steps": 5}, "phi_deg": 0, "psi_deg": 0,
                     "delta_deg": [0, 90]},
       "materials": {"air": {"n": 1.0}, "si": {"n": [3.88, 0.02]}},
       "superstrate": "air", "substrate": "si",
       "layers": [{"type": "profile1d", "shape": "sinus", "depth_um": 0.1, "below": "si", "above": "air"}]}
scan = fg.solve_scan(cfg)                           # 5 wavelengths × 2 polarizations
c = scan.columns
sel = (c["side"] == "R") & (c["n"] == 0)
print(c["lambda_um"][sel], c["delta_deg"][sel], c["efficiency"][sel])
cfg["sweep"] = {"parameter": "/layers/0/depth_um", "values": [0.05, 0.1, 0.15]}
print(fg.solve_scan(cfg).sweep_parameter)
sim = fg.Simulation(cfg)                             # optimization loop: reuse the layer models
L, g = sim.value_and_grad({("R", -1, 0): 1.0}, ["/layers/0/depth_um"], lambda_um=0.6, delta_deg=90)
sim.set_parameter("/layers/0/depth_um", 0.1 + 0.01 * g[0])
print(L, g, sim.solve(lambda_um=0.6, delta_deg=90).efficiency("R", -1))
```
2D pillars (circles of exact shape, with a hole and a square):
```json
{"dimension": "2D", "solver": "rcwa", "truncation": {"Nx": 4, "Nz": 4},
 "period_um": {"x": 0.5, "z": 0.5},
 "incidence": {"lambda_um": 0.633, "phi_deg": 10, "psi_deg": 30, "delta_deg": 45},
 "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "si": {"n": [3.88, 0.02]}},
 "superstrate": "air", "substrate": "glass",
 "layers": [{"type": "contour2d", "thickness_um": 0.15, "background": "air",
             "regions": [{"material": "si", "circle": {"center_um": [0.25, 0.25], "radius_um": 0.12},
                          "holes": [{"circle": {"center_um": [0.25, 0.25], "radius_um": 0.04}}]},
                         {"material": "si", "polygon_um": [[0.0, 0.0], [0.06, 0.0], [0.06, 0.06], [0.0, 0.06]]}]}]}
```

## CAD (STEP; `fg.cad.available()` must be True)
- Model in µm and write with `unit="um"`, one named solid per material. CAD +Z maps to FG y (stacking) and CAD +X to FG x; FG z = −CAD Y.
- Layer: `{"type": "cad", "file": "x.step" or "step": text, "unit": "auto", "mm" or "um", "axes": {"y": "+Z", "x": "+X"},
  "bodies": [{"name": "Line", "material": "si"}, {"index": 1, "material": ...}], "background": "superstrate",
  "cell_origin_um": [x0, z0] (1D: x0), "y_range_um": [y0, y1], "overlap": "error" or "priority",
  "periodic_union": false (1D only), "z_cut_um": z (1D), "slicing": {"rcwa_slices": 40,
  "y_interpolation": "none" or "chebyshev"}, "factorization": "auto", "normal_field": {"grid": 512}}`.
- 2D gratings: the bodies are cut in the (x, z) cell. 1D gratings (`"dimension": "1D"`): every body must be extruded along FG z
  (a line of a line grating); it is cut once at `z_cut_um`, and arcs stay exact.
- Always run `r = fg.cad.check(step_text, cfg)` first. Fix `r["errors"]` and look at `r["warnings"]` and `r["bands"]`.
  Codes: `SHELL_OPEN`, `SHELL_ORIENTATION`, `BODY_OVERLAP`, `PERIODIC_OVERLAP`, `OUT_OF_RANGE`, `UNIT_SUSPECT`,
  `FEATURE_SMALL`, `BODY_UNASSIGNED`, `STEP_UNSUPPORTED`, `NOT_Z_INVARIANT`.
- A body that overlaps its own periodic image is `PERIODIC_OVERLAP` in 2D, even with `periodic_union`; import line
  gratings as 1D. `priority` uses Booleans in 2D; if they fail (e.g. coplanar faces), the check reports
  `BODY_OVERLAP`. Then separate the bodies in CAD.

## Choices and pitfalls
- `rcwa` is fast and robust; y-variable layers become a staircase (`rcwa.slices`, `slices` per layer).
  `rk` integrates y-variable layers without a staircase (`rk.tol`) and is better for deep or slanted profiles.
- Converge N first: metals in TM converge slowly; use `coordinate_transform` ASR with eta 0.8 to 0.95.
  1D, ψ = 0, isotropic: `split_polarizations: true` is faster.
- A z-invariant structure in 2D (Nz = 0) equals 1D conical incidence; use `"dimension": "1D"` with `psi_deg`.
- `lamellar1d.breakpoints_um` must end exactly at `period_um.x`. Geometry may extend past the cell; periodic
  images must not overlap.
- The superstrate must be lossless at every wavelength. Dispersive models are evaluated per scan point.
- `fg.solve` ignores `sweep` and uses only the first scan point; use `fg.solve_scan` for everything.
- `Simulation` is not thread-safe; use one per thread. The GIL is released during computations.
- Efficiencies of evanescent orders are not listed; `efficiency(side, n, m)` returns 0 for them.
