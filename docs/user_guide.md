# Anwenderhandbuch

FourierGrating berechnet Beugungseffizienzen und Nahfelder periodischer Strukturen (1D- und 2D-Gitter,
isotrope und anisotrope Materialien) mit der Fourier-Modal-Methode. Dieses Handbuch beschreibt die
Benutzung; Formeln und Algorithmen stehen im Plan (`docs/implementation_plan/main.pdf`), das
Eingabeformat vollständig in `docs/config_format.md`.

## 1. Installation

```bash
sudo ./scripts/setup_almalinux10.sh            # Abhängigkeiten (--with-cuda für das CUDA-Toolkit)
cmake --preset cpu-release && cmake --build --preset cpu-release
ctest --preset cpu-release                      # ca. 15–20 min auf einem Kern
./build/cpu-release/apps/fgsolve --info
```

Mit GPU: `cmake --preset cuda-release` (CUDA ≥ 12.6 Update 2, siehe `docs/build_almalinux10.md`).
`fgsolve --info` zeigt, ob ein CUDA-Gerät erkannt wurde.

## 2. Erste Rechnung

```bash
./build/cpu-release/apps/fgsolve examples/v3_lamellar_dielectric.json
```

Ausgabe auf dem Terminal (Effizienz je Ordnung, Energiebilanz) und als CSV in `out/…/efficiencies.csv`.
Eine Konfiguration besteht aus:

| Abschnitt | Inhalt |
|---|---|
| `dimension`, `period_um` | `1D` (Periode x) oder `2D` (Perioden x und z) |
| `materials` | Brechzahl `n` (reell oder `[Re, Im]`), Tensor `eps_tensor` oder Hauptachsen + Euler-Winkel |
| `superstrate`, `substrate` | Einfallsmedium (oben) und Substrat (unten) |
| `layers` | Schichten von unten nach oben (siehe 3.) |
| `incidence` | λ in µm, φ (Polarwinkel zur Normalen y), ψ (Azimut gegen x), δ (0 = s/TE, 90 = p/TM) |
| `solver`, `truncation` | `rcwa` oder `rk`; Ordnungen −Nx…Nx (und −Nz…Nz in 2D) |

Koordinaten: y ist die Schichtnormale (y = 0 an der Substratgrenze, nach oben), x und z liegen in der
Gitterebene. Effizienzen beziehen sich auf den einfallenden Poynting-Fluss; die s/p-Anteile werden in
isotropen Medien zusätzlich ausgegeben.

## 3. Schichttypen

| Typ | Beschreibung | 1D | 2D |
|---|---|---|---|
| `homogeneous` | homogene Schicht | ✓ | ✓ |
| `lamellar1d` | Lamellen (Breakpoints in x, Material je Intervall), z-invariant | ✓ | ✓ |
| `profile1d` | Oberfläche y = g(x): `sinus`, `triangle` (Echelette), `points` (Spline) | ✓ | ✓ (z-invariant) |
| `polygon2d` | y-invarianter Querschnitt: Hintergrund + Polygone/Rechtecke | – | ✓ |
| `profile2d` | Oberfläche y = g(x, z): `egg_crate`, `sinus_x` | – | ✓ |

Profilschichten werden vom RK-Löser direkt integriert, von RCWA in `slices` Treppenstufen zerlegt.

## 4. Löser und Faktorisierung wählen

* **RCWA** (`"solver": "rcwa"`) ist für y-invariante Schichten (Lamellen, Polygone, homogen) erste Wahl.
  Für Profile ist die Treppenzerlegung ein zusätzlicher Fehler (∝ 1/`slices`).
* **RK** (`"solver": "rk"`) integriert Profile ohne Treppe; Genauigkeit über `rk.tol` (Standard 1e-8).
  In 2D teuer (Abschnitt 8), daher vor allem für Profile und Referenzrechnungen.
* **Faktorisierung** (`"factorization"` je Schicht; `auto` wählt den Standard): Lamellen `LiLamellar`,
  glatte Profile `FFF_rig`, Profile mit Ecken und Funktionen `FFF_sym`, 2D-Polygone `FFF_sym`, anisotrop
  `FFF_A` (Profile) bzw. `FFF_B` (Polygone). `laurent` dient Vergleichen; `FFF_iso` (Buchformel) versagt bei
  schwach absorbierenden Metallen (Warnung).
* **Konvergenz prüfen:** Ergebnis bei zwei bis drei Werten von N vergleichen. Richtwerte: 1D-Dielektrika
  N = 10–20, 1D-Metalle N = 20–40, 2D N = 5–10 je Richtung. Resonanzen (z. B. Plasmonen in Lochgittern)
  verschieben sich mit N – dort Spektren statt Einzelwerte vergleichen (`docs/validation.md`, V10).

## 5. Scans und Sweeps

Jeder Einfallsparameter darf eine Zahl, eine Liste oder ein Bereich `{"from", "to", "steps"}` sein; gerechnet
wird das kartesische Produkt λ × φ × ψ × δ, parallel auf `OMP_NUM_THREADS` Threads.

Beliebige Zahlenwerte der Konfiguration lassen sich per JSON-Pointer variieren:

```json
"sweep": { "parameter": "/layers/0/thickness_um", "values": { "from": 0.05, "to": 1.2, "steps": 24 } }
```

Beispiele für Pointer: `/layers/1/depth_um`, `/layers/0/breakpoints_um/1`, `/period_um/x`,
`/materials/si/n/0`, `/layers/0/inclusions/0/rectangle_um/size/0`. Die CSV erhält dann die Spalten
`sweep_parameter`, `sweep_value`. Darstellung:

```bash
fgsolve examples/sweep_lamellar_depth.json
scripts/plot_efficiencies.py out/sweep_lamellar_depth/efficiencies.csv --orders R:0 T:-1 T:0 --out sweep.png
```

## 6. Nahfeld

```json
"nearfield": { "enabled": true, "horizontal_y_um": 0.1, "vertical_z_um": 0.0,
               "y_range_um": [-0.3, 0.6], "nx": 400, "ny": 300, "nz": 200, "periods": 2,
               "reconstruction": "fff" }
```

* vertikaler Schnitt: Karte (x, y) bei z = z0; horizontaler Schnitt: Linie (1D) bzw. Karte (x, z) (2D) bei y = y0.
* `reconstruction: "fff"` (Standard) hält E-Tangential- und D-Normalkomponente an Materialgrenzen stetig;
  `"direct"` zeigt Gibbs-Schwingungen (Sprung von D_x an Wänden ≈ 300 %).
* Ausgabe in `out/…/nearfield_horizontal|vertical/`: `Ex…Hz.npy` (H als Z0·H), `absE2`, `absH2`, `Sy`,
  `material.npy`, `meta.json`; Darstellung mit `scripts/plot_nearfield.py`.

## 7. Ausgabe und Diagnose

* `efficiencies.csv`: eine Zeile je Ordnung und Scanpunkt; Spaltenliste in `docs/config_format.md`.
* Terminal: `sum R+T` und `absorption`. Bei verlustfreien Strukturen muss die Summe 1 sein (RCWA mit Laurent
  exakt, FFF-Varianten auf 1e-5 bis 5e-4). Negative Absorption oder R > 1 zeigen eine ungeeignete
  Faktorisierung oder zu kleines N an.
* Warnungen (z. B. `FFF_iso` mit Metall) erscheinen in der Ausgabe und in `Result::diag.warnings`.

## 8. Laufzeiten und Leistung

Richtwerte (eine Wellenlänge, ein CPU-Kern der Sandbox, v0.10; Zielsystem: `docs/benchmarks.md`):

| Fall | Größe | Zeit |
|---|---|---|
| 1D-Lamellen (V3), RCWA, N = 15 | 4 n_f = 124 | 0,04 s |
| 1D-Al-Sinus (V4), RK, N = 20, tol 1e-5 | 4 n_f = 164 | ≈ 37 s |
| 2D-Ag-Lochgitter, RCWA, N = 4 / 6 / 8 | 4 n_f = 324 / 676 / 1156 | 0,8 / 3,0 / 13,5 s |
| 2D-Eierkarton, RCWA 20 Scheiben, N = 5 | 4 n_f = 484 | ≈ 55 s |

Die Rechenzeit wächst etwa mit n³ (n = 4 n_f, n_f = (2Nx+1)(2Nz+1)); das Eigenwertproblem der RCWA dominiert.
Hinweise:

* Scans parallelisieren: `OMP_NUM_THREADS=<Kerne> fgsolve scan.json`.
* GPU: `"backend": "auto"` nutzt die GPU ab `backend_threshold` (Standard n = 1000), `"cuda"` immer,
  `"cpu"` nie. Lohnt sich vor allem in 2D ab N ≈ 8. Die Umgebungsvariablen `FG_BACKEND` und
  `FG_BACKEND_THRESHOLD` überschreiben die Konfiguration (z. B. `FG_BACKEND=cpu` zum Vergleich).
* 1D, ψ = 0 und nur isotrope Medien: `"split_polarizations": true` rechnet TE und TM getrennt mit 2 n_f
  (Faktor 2–5 bei RCWA, siehe `docs/benchmarks.md`); reines TE (δ = 0°) oder TM (δ = 90°) ist dann nochmals
  schneller. Nahfelder nutzen immer das volle System. Standard ist aus (bitgleich zu v0.10).
* y-invariante 2D-Schichten mit RCWA rechnen; RK in 2D nur für Profile.
* Für Profile mit RCWA genügen oft 20–40 Scheiben; Konvergenz in `slices` prüfen.

## 9. Programmierschnittstelle (C++)

```cpp
#include "fg/solvers/solver.hpp"
fg::Stack st;  st.dx_um = 0.5;
st.materials["air"] = fg::Material::fromIndex("air", 1.0);
st.materials["glass"] = fg::Material::fromIndex("glass", 1.5);
st.superstrate = "air";  st.substrate = "glass";
fg::LayerEntry le;  le.spec = fg::Lamellar1DLayer{0.3, {0.0, 0.25, 0.5}, {"glass", "air"}};
st.layers.push_back(le);
fg::SolverOptions o;  o.Nx = 15;
const fg::Result r = fg::makeSolver(o)->solve(st, {0.6328, 0.3, 0.0, 0.0});  // λ, φ, ψ, δ (rad)
```

Weitere Einstiegspunkte: `solveScan` (parallele Scans), `fieldCoefficients` / `nearField` (Nahfeld),
`fg::io::loadSweep` (Konfigurationen). Beispiele in `tests/validation/`.

## 10. Häufige Probleme

| Meldung / Symptom | Ursache, Abhilfe |
|---|---|
| Einschlüsse verschiedener Materialien berühren sich (`polygon2d`) | mit `FFF_B` unzuverlässig (Vorzeichen des Normalenfelds dort undefiniert, nicht geprüft) → Einschlüsse trennen oder `laurent` |
| Fehlermeldung zur Faktorisierung | Variante passt nicht zum Schichttyp (z. B. `FFF_rig` bei Polygonen) → `auto` |
| R + T > 1 bei Metallen | `FFF_iso` gewählt → `auto`/`FFF_sym`; oder N erhöhen |
| Ergebnis springt mit N | Resonanz; Spektrum rechnen, Peak-Lage/-Höhe vergleichen |
| RK sehr langsam | 2D oder große N; RCWA mit Treppe verwenden oder Toleranz lockern |
| `CUDA backend requested but not available` | ohne GPU gebaut/kein Gerät; `"backend": "auto"` |
