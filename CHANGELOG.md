# Changelog

## [Unreleased]

### Plan
- Plan v2.2: Optionen `split_polarizations` (Spezialfälle) und `fft_sampling` (Allgemeine Funktionen).
- Plan v2.0: Erweiterungsplan Version 2 (Teile 5–6, Meilensteine M12–M19): Python-Anbindung, Konturkoeffizienten und
  höhenabhängige Schichtmodelle, CAD-Slicer mit Nurbforge, Profiling/Parallelisierung, wählbare Genauigkeit
  (FP64/gemischt/FP32), Koordinatentransformationen und ASR, PML für isolierte Strukturen, diskrete Adjungierte.
  Roadmap, `AGENTS.md` und `CLAUDE.md` entsprechend ergänzt.
- Plan v2.1: Entscheidungen des Projektinhabers eingearbeitet (isolierte Strukturen mit Streufeldformulierung für
  schrägen Einfall, Windows-Wheels mit MSVC und CUDA-Variante, Geometrieprüfung für LLM-erzeugte Nurbforge-Geometrie,
  Einheit µm, FP32-Standard auf der RTX 3090).

### Hinzugefügt
- MSVC-Port (M12): Presets `win-msvc` und `win-msvc-cuda` (`86;89`), `scripts/msvc.cmd` (vcvars64), MSVC-Flags
  `/W4 /permissive- /utf-8 /bigobj`, externe Header `/external:W0`; MSYS2-Präfixe und Paketregistry im Preset
  ausgeschlossen. Alle Tests grün, `Backend.*` mit CUDA 13.4 auf der RTX 3090 grün. `docs/build_windows.md` Abschnitt 4.
- Golden Files `tests/regression/golden/` und Test `Regression.*` (Label `regression`): Effizienzen aller Beispiele mit
  Standardoptionen; 1e-13 bezogen auf die einfallende Leistung (Plan v2.3; rein relativ scheitern sehr kleine
  Ordnungen beim Compilerwechsel), `FG_REGRESSION_EXACT=1` bitgleich, `FG_REGRESSION_UPDATE=1` erzeugt neu.
- Option `"split_polarizations"` (`SolverOptions::splitPolarizations`, Standard aus): 1D in-plane isotrop
  (γ = 0, isotropes Substrat und isotrope Schichten) rechnet `solve()` TE und TM als getrennte 2n_f-Systeme
  (Plan eq:TE, eq:TM) – Grenzflächen, Schicht-S-Matrizen, RK-Transfer und Sternprodukte mit halber Größe;
  Einbettung blockdiagonal, nicht angeregte Polarisation (reines TE/TM) übersprungen. RCWA Faktor 1,9–5.
  Neu: `Polarization`, `polarizationBlocks`, `polarizationBasis`, `buildPropagationMatrix(Q, alpha, Polarization)`,
  `layerModesPolarization`.
- Option `"fft_sampling"` (`SolverOptions::fftSampling`, Standard aus): `fourier::fft` (beliebige Länge, Radix 2 /
  Bluestein) für `sampledCoefficients1D/2D` und `gridCoefficients2D`; Schalter `fourier::setFFTSampling` (thread-lokal).
- Tests `PolarizationSplit.*`, `Sampling.FFTEqualsDirectDFT`, `Sampling.FFTSwitchMatchesDFT`,
  `Sampling.SampledCoefficientsExactForTrigPolynomials`.
- Mit Standardoptionen sind alle `examples/` bitgleich zum Stand vor der Änderung (CSV-Vergleich).
- Umgebungsvariablen `FG_BACKEND` / `FG_BACKEND_THRESHOLD` überschreiben die Backend-Wahl (ganze Testsuite auf
  der GPU: `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest`); Test `Backend.EnvironmentOverride` (portabel, auch MinGW).
- V10 bei N = 11: Peak ≈ 1,3898 µm, ≈ 9,2 %.

### Geändert
- `layerModesAuto`: Klassifikation der TE/TM-Blöcke in eine gemeinsame Funktion ausgelagert (bitgleich).

## [0.10.0] – M11 (Scans, Sweeps, Anwenderdoku)

### Validierung (nachgetragen)
- V10 bei N = 9 und 10: Transmissionspeak bei 1,3876 / 1,3889 µm (Buch: 1,388 µm), Höhe 4,2 / 6,2 %
  (Buch 7,5 %); Tensor- und Vektorvariante von Algorithmus NF gleichwertig (N = 7, 8).

### Hinzugefügt
- Scans über λ, φ, ψ und δ (Zahl, Liste oder Bereich; kartesisches Produkt), parallel über `solveScan`.
- Parameter-Sweeps über beliebige Zahlenwerte der Konfiguration (`"sweep": {"parameter": <JSON-Pointer>,
  "values": ...}`), `fg::io::parseSweep`/`loadSweep`, CSV mit `sweep_parameter`, `sweep_value`.
- `scripts/plot_efficiencies.py` (Spektren und Sweeps aus `efficiencies.csv`).
- Beispiel `examples/sweep_lamellar_depth.json` (Smoke-Test), Tests `Scan.*`.
- Anwenderhandbuch `docs/user_guide.md` mit Laufzeittabelle; README-Status aktualisiert.

### Geändert
- `fgsolve` rechnet Sweeps (ein Block je Wert) und Scans parallel.
- Plan v1.9.

## [0.9.0] – M10 (CUDA-Backend, ohne GPU-Abnahme)

### Hinzugefügt
- `CudaBackend` (`src/cuda/cuda_backend.cpp`): cuBLAS `Zgemm`, cuSOLVER `Xgetrf`/`Xgetrs` (64 Bit) und
  `Xgeev` über die Host-API; RAII für Device- und Pinned-Speicher, Stream je Backend. Kein `nvcc` nötig.
- Dispatch (`linalg::setPreference`, `offload`, `device`): RCWA-Eigenprobleme, stabile Schicht-S-Matrix,
  Sternprodukt, T→S und RK-Schrittpropagator laufen ab der Schwelle auf der GPU; CPU-Pfad unverändert.
- `SolverOptions::backend` / `backendThreshold`, JSON `"backend"`, `"backend_threshold"`.
- `solveScan`: paralleler Wellenlängen-/Winkelscan (OpenMP, ein Löser und Stream je Thread); `fgsolve` nutzt ihn.
- `fgbench` (Benchmark-Tabellen), `docs/benchmarks.md` (Sandbox-Referenz, Zielsystem offen).
- `fgsolve --info` meldet „CUDA gebaut, kein Gerät“ statt abzubrechen.
- Tests `Backend.*` (GPU-Vergleiche werden ohne Gerät übersprungen).

### Geändert
- CMake: CUDA ohne `enable_language(CUDA)` (nur `find_package(CUDAToolkit 12.6)`), `cuda_backend.cu` → `.cpp`.
- Plan v1.8.

## [0.8.0] – M9 (Nahfeld 2D)

### Hinzugefügt
- `Solver::nearField` für 2D-Gitter: horizontale Karte u(x, y0, z) (`nz` Punkte je Periode, Achse `z`) und
  vertikaler Schnitt u(x, y, z0); Synthese mit separierten Phasen.
- FFF-konsistente Rekonstruktion in 2D: `LayerModel::projectorAt(x, z)`, `rank1Projector()`;
  isotrope Polygone mit dem Tensorfeld (E = e_t + d_n/ε), Variante B / anisotrop mit dem Rang-1-Feld,
  Profile2D mit eq:N2d; Punktwerte aus dem NF-Raster (bilinear, vorzeichenfrei).
- `epsAt(x, y, z)` und `normalAt(x, z)` für alle Schichtmodelle (Materialkarten der 2D-Schnitte).
- Validierung V12 (2D): z-invariant = 1D, Stetigkeit, Flusskonstanz, RK = RCWA, D_x-Sprung an Säulenwänden
  (isotrop und anisotrop), Schnitte und Dateien; Beispiel `examples/nearfield_pillars_2d.json` (Smoke-Test).

### Geändert
- Plan v1.7.

## [0.7.0] – M8 (2D anisotrop)

### Hinzugefügt
- Anisotrope Materialien in `Polygon2D` (Laurent exakt; `FFF_B`, Standard) und `Profile2D`
  (`FFF_A`, Standard; `FFF_B`; Laurent); anisotrope 1D-Profile in 2D-Stapeln.
- Variante B für Querschnitte: Rang-1-Normalenfeld aus Algorithmus NF (normierter geglätteter Gradient,
  Hauptachse des Tensors als Rückfall), Koeffizienten als glatter Anteil (FFT) plus Faltung der
  Materialdifferenz mit den exakten Polygonkoeffizienten (`fourier::convolveCoefficients`).
- `FFF_B` auch für isotrope `Polygon2D`, `FFF_A`/`FFF_B` auch für isotrope `Profile2D` (Vergleiche).
- Radix-2-FFT in `gridCoefficients2D` (komplexe Variante); NF-Raster auf Zweierpotenzen.
- Validierung V8 (2D), Spiegeltest mit anisotropem Tensor, Konvergenzstudie Co-Säulen und
  anisotroper Eierkarton.

### Geändert
- Anisotrope Faktorisierungen intern vereinheitlicht (`anisoEval`, `anisoQFromGrids`), 1D unverändert.
- Plan v1.6.

## [0.6.0] – M7 (2D isotrop)

### Hinzugefügt
- `Polygon2D`-Schichten: exakte Polygon-Koeffizienten von ε und 1/ε (Polygone dürfen den Zellrand kreuzen,
  Orientierung wird normalisiert); Faktorisierung `FFF_sym` (Standard), `FFF_iso`, Laurent.
- Algorithmus NF (Plan alg:nfield, v1.5): exakte Abstände zu den Materialgrenzen (periodische Bilder, innere
  Kanten entfernt, Mittelachsen gemittelt), periodische Gaußglättung, DFT; Tensorvariante (Standard) und
  Vektorvariante (v1.4) wählbar (`normal_field.method`).
- `Profile2D`-Schichten (`egg_crate`, `sinus_x`, Funktion): ε(y) zeilenweise exakt in x, in z Trapezregel bzw.
  Gauß–Legendre zwischen kritischen Linien; Normalenfeld nach eq:N2d (vorzeichenbehaftet, `FFF_rig` Standard).
  RK (M(y)) und RCWA-Treppe.
- Isotrope `profile1d`-Schichten in 2D-Stapeln (z-invariant).
- JSON: `polygon2d` (`vertices_um`, `rectangle_um`, `normal_field`), `profile2d`; Prüfung `period_um.z > 0` in 2D.
- `fourier::gridCoefficients2D` (separable DFT eines Rasters).
- Validierung V9 (z-invariante Strukturen = 1D konisch), V10 (Ag-Lochgitter, qualitativ), Spiegeltest x↔z,
  Eierkarton (Energie, RK = RCWA); Unit-Tests Algorithmus NF und 2D-Schichtmodelle.

### Geändert
- Plan v1.5: Algorithmus NF mit Tensorglättung (Befund: die Vektorglättung fällt an Mittelachsen auf ½ zurück
  und reproduziert Lamellen nicht exakt); Profile2D ohne Marching Squares.
- RK-Spektralradius für die Scheibenzerlegung: oberhalb n = 800 Potenziteration statt voller Eigenzerlegung.

### Befund
- V10 (Ag, ε = −79,9 + 1,79i): `FFF_iso` liefert R > 1 (N = 3, 6); der Transmissionspeak oberhalb des
  SPP-Pols wandert mit N nach Rot (1,383 → 1,386 → 1,391 → 1,404 µm für N = 5…8), Höhe 3–10 %;
  quantitativer Vergleich mit dem Buch (7,5 %) offen.

## [0.5.0] – M6 (Nahfeld 1D)

### Hinzugefügt
- Gemeinsame Schichtmodelle (`detail::LayerModel`) für Lösung und Feldrekonstruktion.
- Aufzeichnende Assemblierung: interne Amplituden (S_L/S_U), modale Scheiben, RK-Transfermatrizen
  an den angeforderten Höhen.
- `Solver::fieldCoefficients` und `Solver::nearField`: horizontaler und vertikaler Schnitt,
  direkte oder FFF-konsistente Rekonstruktion (Projektorform, iso und aniso).
- NPY/JSON-Ausgabe (`NearFieldWriter`), `fgsolve` schreibt Schnitte bei `nearfield.enabled`.
- Validierung V12 (Airy analytisch, Stetigkeit, Flusskonstanz, RK = RCWA, D_x-Stetigkeit).

## [0.4.0] – M5 (anisotrop, 1D)

### Hinzugefügt
- Anisotrope FFF: Variante A (Buch Kap. V/VIII, C_ε) und Variante B (Projektoren, Eigenherleitung);
  Fourierkoeffizienten stückweise glatter Funktionen per Gauß-Quadratur.
- Anisotrope Modenbasis homogener Medien (4×4 je Ordnung, Klassifikation); anisotrope Substrate,
  homogene Schichten, Lamellen (N = x konstant) und Profile; Transmission über den Modenfluss.
- Materialeingabe über Hauptachsen + Euler-Winkel; s/p-Anteile der Effizienzen (CSV-Spalten eff_s, eff_p).
- Validierung V7 (Kobalt), V8 (isotroper Tensor), einachsiges Substrat (analytisch), biaxiale Schicht.

### Behoben
- Eigen konjugiert bei `dot()`/`cross()` komplexer Vektoren: ξ und C_ε waren falsch, s/p-Amplituden konjugiert.

## [0.3.0] – M3, M4

### Hinzugefügt
- Profilschichten y = g(x): Sinus (exakte Schnittpunkte), stückweise linear (Dreieck, Echelette, Punkte;
  exakte Koeffizienten), Funktion (Nullstellensuche); Normalenfeld (NN und vorzeichenbehaftet).
- RK mit M(y) (exakte Koeffizienten an jeder Stufenhöhe), RCWA-Treppe mit glattem Normalenfeld.
- FFF-Varianten `FFF_sym` und `FFF_rig` (Plan eq:Qsym, eq:Qrig); Standard je nach Profilform;
  Warnung für die Buchvariante bei schwach absorbierenden Metallen.
- TE/TM-Entkopplung im in-plane-Fall (RK ≈ 6×, RCWA ≈ 2,5× schneller).
- Validierung V4 (Al-Sinus gegen C-Methode), V5, V6, konische Profile, V13 (Gold-Sinus gegen C-Methode),
  V14 (Gold-Dreieck, offen); Paarungs- und Profil-Unit-Tests.

### Behoben
- RK-Scheibenteilung: Spektralradius aus der abgeschnittenen Matrix statt Materialschätzung,
  plus Halbierung bei zu großem Wachstum (vorher Kontamination bei Gold).
- Modenklassifikation in metallischen Scheiben (Sortieren und Halbieren statt Schwellenwert).

### Befund
- Buch (IV.32) liefert für Gold unphysikalische Effizienzen (Σe > 1); FFF erhält die Energie nicht exakt.

## [0.2.0] – M1 abgeschlossen, M2

### Hinzugefügt
- Q-Matrizen: Laurent, Li-Lamellen, FFF isotrop (Plan eq:Qiso); Abtastung stetiger Funktionen.
- RCWA: reduziertes und allgemeines Eigenwertproblem mit Klassifikation, stabile Scheiben-S-Matrix.
- RK: Transfermatrix für konstantes M (Schrittverdopplung) und adaptiv für M(y); Scheiben-S aus T.
- Stapellöser für beide Verfahren: Substrat-/Referenzgrenzfläche, homogene und Lamellenschichten,
  Effizienzen aus dem Poynting-Fluss, s/p-Amplituden, Energiebilanz.
- JSON-Konfiguration (homogeneous, lamellar1d, Wellenlängen-/Winkel-Scans), `fgsolve` mit CSV-Ausgabe.
- Validierung V1 (Fresnel, konisch, verlustbehaftet), V2 (Airy, RK und RCWA), V3 (Energie,
  RK = RCWA, Reziprozität); 41 Tests.

### Geändert
- Plan v1.1: volles 4n_f-System auch in-plane (2n_f-Systeme als spätere Optimierung);
  Reziprozität im FMM nur bis auf Abschneidung.

## [0.1.0] – Grundgerüst

### Hinzugefügt
- Implementierungsplan (LaTeX/PDF) mit allen Formeln, Herleitungen, Tests und Meilensteinen.
- CMake-Projekt (Presets CPU/CUDA), Abhängigkeiten (Eigen3, nlohmann_json, GoogleTest), CI (CPU).
- Implementiert und getestet: `IndexMap` (korrigiert Buch IX.14), exakte Fourierkoeffizienten
  (Intervall, Rechteck, Polygon), Toeplitz-Matrizen (1D, 2D), S-Matrix (T→S, Sternprodukt,
  interne Amplituden), Butcher-Tableaus (Gauß 2/3, Radau IIA 3) und RK-Schrittpropagator,
  Propagationsmatrix `M` (μ = 1), isotrope Modenbasis, einfallende Welle (korrigiert Buch IX.1).
- Schnittstellen (noch `NotImplemented`): Q-Matrizen (FFF), anisotrope Moden, Löser RK/RCWA,
  Nahfeld, Konfiguration, CUDA-Backend.
