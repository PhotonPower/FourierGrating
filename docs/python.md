# Python-Anbindung (M12)

Das Paket `fouriergrating` kann alles, was `fgsolve` kann, ohne Umweg über Dateien. Die Konfiguration ist
dasselbe JSON-Schema wie für `fgsolve` ([`config_format.md`](config_format.md)), als Python-`dict`; CLI und
Python teilen sich Parser und Semantik. Einheiten wie in JSON: µm und Grad.

## Installation

| Wo | Befehl |
|---|---|
| Windows 11, CPU | `pip install fouriergrating` (Wheel, MSVC-Build) |
| Windows 11, RTX 3090 | zusätzlich `pip install fouriergrating-cuda12` (zieht `nvidia-*-cu12`) |
| Linux (manylinux_2_28, z. B. AlmaLinux 8–10), CPU | `pip install fouriergrating` |
| Server AlmaLinux 10 mit CUDA | aus den Quellen, siehe unten |
| Repository (Entwicklung) | `pip install .` |

Aus den Quellen mit CUDA (Server, RTX 6000 Ada):

```bash
pip install . -C cmake.define.FG_ENABLE_CUDA=ON -C cmake.define.CMAKE_CUDA_ARCHITECTURES="89"
```

Unter Windows braucht der Quellbuild MSVC 2022 (Python ist mit MSVC gebaut), z. B.
`scripts\msvc.cmd python -m pip install .`; die CUDA-Variante: `scripts\msvc.cmd python -m pip install python/cuda12`.

`import fouriergrating` lädt den CUDA-Kern (`fouriergrating_cuda12._core`), wenn er installiert und eine GPU
nutzbar ist, sonst den CPU-Kern; `fg.core_variant()` meldet `"cpu"` oder `"cuda12"`. Die CUDA-DLLs werden nicht
mitgeliefert: Unter Windows werden die Verzeichnisse der Pakete `nvidia-*-cu12` und, falls gesetzt, `CUDA_PATH`
vor dem Import eingetragen. `FG_CORE=cpu` erzwingt den CPU-Kern.

## API

```python
import fouriergrating as fg

cfg = fg.load_config("examples/v3_lamellar_dielectric.json")   # dict; output.dir relativ zur Datei
res = fg.solve(cfg)                  # erster Scanpunkt ("sweep" wird ignoriert) -> fg.Result
res.reflected                        # numpy structured array: n, m, eff, eff_s, eff_p, amp_s, amp_p
res.efficiency("T", 0)               # Effizienz einer Ordnung (0, falls nicht propagierend)
res.transmitted, res.absorption
res.diagnostics                      # energy_balance, max_condition_t11, rk_steps, warnings, timings
fg.solve(cfg, profile=True).diagnostics["timings"]["scopes"]   # Zeit je Abschnitt (fg::prof, M15)
res.diagnostics["precision"]         # tatsächlich verwendete Genauigkeit (M16), "precision_fallback"

scan = fg.solve_scan(cfg)            # alle Sweep-Werte und Scanpunkte wie fgsolve (parallel, "threads")
scan.points                          # Liste von fg.Result (mit .incidence)
scan.columns                         # Spalten wie efficiencies.csv als numpy-Arrays
scan.to_pandas()                     # benötigt pandas

nf = fg.near_field(cfg)              # Block "nearfield": {"horizontal", "vertical"}
nf["vertical"]["Ex"], nf["vertical"]["absE2"], nf["vertical"]["axis1_um"], nf["vertical"]["meta"]

sim = fg.Simulation(cfg)             # Schichtmodelle (Fourierkoeffizienten, Normalenfelder) einmal aufbauen
sim.solve(lambda_um=0.6328, phi_deg=30.0)       # fehlende Werte: erster Scanpunkt der Konfiguration
sim.set_parameter("/layers/0/thickness_um", 0.25)   # JSON-Pointer wie bei Sweeps
sim.model_builds, sim.config

fg.set_backend("auto", threshold=1000)   # Standard für Konfigurationen ohne "backend"
fg.cuda_available(); fg.core_variant(); fg.__version__
```

* `set_parameter` baut nur die geänderten Schichten neu, wenn sich außerhalb von `/layers` (außer der
  Einfallsrichtung) nichts ändert; sonst alle Schichtmodelle. Das Ergebnis ist bitgleich zu einer neuen
  `Simulation` aus der geänderten Konfiguration (Test `test_set_parameter_equals_rebuild`).
* Während jeder Rechnung ist der GIL freigegeben. Mehrere Python-Threads dürfen parallel rechnen, jeder mit
  eigener `Simulation` (ein Objekt ist nicht threadsicher).
* Fehler: `fg.FourierGratingError` (Unterklasse von `RuntimeError`) für Konfigurations- und Rechenfehler,
  `NotImplementedError` für spezifizierte, aber noch nicht umgesetzte Funktionen.
* Typannotationen: `py.typed`, `_core.pyi`.

Beispiel: [`examples/python/scan_lamellar.py`](../examples/python/scan_lamellar.py) (Wellenlängen-Scan und
Tiefen-Sweep mit `Simulation`).

### Gradienten (M19)

```python
sim = fg.Simulation(cfg)
params = ["/layers/0/thickness_um", "/layers/0/breakpoints_um/1", "/materials/bar/n/1"]
value, grad = sim.value_and_grad({("R", -1, 0): 1.0, ("T", 0, 0): -0.5}, params)   # L = Σ w·e
v = sim.vjp([("R", -1, 0, 1.0)], params)          # Kotangenten der Effizienzen → Parameter
```

Diskrete Adjungierte (exakter Gradient der tatsächlich berechneten Abbildung): Aufwand etwa 1,3–2,5
Vorwärtsrechnungen, unabhängig von der Zahl der Parameter. Parameter sind Zahlen der Konfiguration (JSON-Pointer
wie bei `set_parameter`): Schichtdicken, `breakpoints_um/j` (Lamellen), `depth_um`/`x0_um` (Sinusprofil),
alle Zahlen einer `polygon2d`- oder `contour2d`-Geometrie (Ecken, NURBS-Kontrollpunkte und -Gewichte, Kreis),
`materials/<name>/n` bzw. `.../n/0|1` (Real-/Imaginärteil) und `.../eps`, auch für Superstrat (nur Realteil, es
bleibt verlustfrei) und Substrat, sowie `incidence/lambda_um`, `phi_deg`, `psi_deg`, `delta_deg` (am aktuellen
Scanpunkt, auch wenn dort eine Liste steht). Anisotrope Materialien: `.../eps_tensor/i/j/0|1`, `.../eps_principal/k`, `.../euler_deg/k` in
homogenen Schichten, Lamellen und Profilen (Laurent, `FFF_A`, `FFF_B`) und in 2D-Querschnitten mit Laurent. RCWA
und RK, periodische Gitter, auch mit Koordinatentransformation (ASR; Geometrieparameter dann nur mit expliziten
`breakpoints_um`, λ, Winkel und Superstrat/Substrat nicht); anisotrope 2D-Querschnitte auch mit `FFF_B`;
`profile2d` mit anisotropen Medien und aperiodische Strukturen: noch nicht (Fehlermeldung). Festgelegt (Plan M19 „Methode“): Scheibenzerlegung und RK-Schrittfolge; das Normalenfeld der Profile wird
mitdifferenziert, das Rasterfeld der 2D-Querschnitte (Algorithmus NF, FFF) bei Formparametern über eine Rasterzelle
(je Formparameter zwei NF-Rechnungen; unterhalb der Rasterzelle ist die Vorwärtsrechnung nicht glatt), das der
CAD-Schnitte bleibt eingefroren. CAD-Körper: C++-API
`fg::gradient(opt, stack, wl, objective, {StackParameter{...}})` mit dem neu geschnittenen Körper bei θ ± δ
(seitliche Bewegungen). `sim.last_result` enthält die Effizienzen der Rechnung.

Adapter (optional importierbar):

```python
from fouriergrating.torch import efficiencies           # torch.autograd.Function
x = torch.tensor([0.12, 0.22], dtype=torch.float64, requires_grad=True)
e = efficiencies(sim, params[:2], x, [("R", -1, 0)]); (-e[0]).backward()

from fouriergrating.jax import efficiency_function       # jax.custom_vjp + pure_callback
f = efficiency_function(sim, params[:2], [("R", -1, 0)]); jax.grad(lambda x: f(x)[0])(x0)
```

Beispiele: [`optimize_littrow.py`](../examples/python/optimize_littrow.py) (L-BFGS-B, Littrow −1. Ordnung, Vergleich
mit dem Rasterscan) und [`optimize_nurbs_pillar.py`](../examples/python/optimize_nurbs_pillar.py) (PyTorch-Adam auf
NURBS-Kontrollpunkten).

## Tests

```bash
pip install .[test]
pytest tests/python                         # Beispiele = Golden Files (1e-11), Simulation, Threads, GIL
FG_FGSOLVE=build/cpu-release/apps/fgsolve pytest tests/python   # zusätzlich Nahfeld = fgsolve-Dateien
                                            # (fgsolve mit demselben Compiler wie das Modul: 1e-12)
```

## Aufbau

* `python/bindings/module.cpp`: C++-Modul `_core` (pybind11); nur Umwandlung und Aufrufe, keine Physik
  (AGENTS.md, Abschnitt 7). Die Rechnung steckt in `fg::Simulation` (`include/fg/solvers/simulation.hpp`).
* `python/fouriergrating/`: Python-Paket (`__init__.py`, `_loader.py` für die Kernwahl).
* `python/cuda12/`: zweite Distribution `fouriergrating-cuda12` (nur `fouriergrating_cuda12._core`, mit
  `FG_ENABLE_CUDA=ON`, Architekturen `86;89`).
* `pyproject.toml` (scikit-build-core): setzt `FG_BUILD_PYTHON=ON`, ohne Tests und Programme; nur die
  CMake-Komponente `python` landet im Wheel.
* Wheels: `.github/workflows/wheels.yml` (cibuildwheel, `v*`-Tags).
