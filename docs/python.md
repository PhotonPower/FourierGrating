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
