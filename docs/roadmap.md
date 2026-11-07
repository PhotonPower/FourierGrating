# Roadmap

Abnahme = alle genannten Tests grün (`ctest`), Details im Plan, Abschnitt „Meilensteine“.

- [x] **M0** Build-System, Abhängigkeiten, CI (CPU), Grundgerüst, Doku
  - Abnahme: `cmake --build` + `ctest` grün (CI grün; AlmaLinux 10 noch zu bestätigen)
- [x] **M1** Kernbausteine
  - [x] IndexMap, Fourierkoeffizienten, Toeplitz, Butcher + RK-Propagator, S-Matrix, isotrope Moden, `M`, Einfall
  - [x] Abtastung stetiger Funktionen (direkte DFT; FFT-Beschleunigung später, M11)
  - [x] Effizienzen aus Poynting-Fluss (Plan eq:eff), s/p-Zerlegung
  - [x] Q-Matrizen: Laurent, LiLamellar, FFF_iso
- [x] **M2** 1D isotrop: homogene Schichten, Lamellen (RCWA und RK), Gesamtlösung, JSON-Konfiguration, CSV-Ausgabe
  - Abnahme: V1, V2, V3 (TE/TM **und konisch**) – `tests/validation/test_v1_v3.cpp`
  - Hinweis: Es wird durchgehend das volle 4n_f-System verwendet (Plan, Abschnitt Spezialfälle)
- [x] **M3** 1D-Profile: Schnittpunkte, Normalenfeld, FFF-Varianten (iso/sym/rig); RK mit M(y); RCWA-Treppe mit FFF
  - Abnahme: V4, V5, V6, V13 grün; V14 (Gold-Dreieck) nur physikalische Grenzen – **offen**
  - Befund: Buchformel (IV.32) versagt bei schwach absorbierenden Metallen, siehe Plan „Varianten der isotropen FFF“
- [x] **M4** 1D konisch (4n_f-System)
  - [x] Lamellen und homogene Schichten konisch (mit M2 erledigt, V1–V3 konisch grün)
  - [x] Profilschichten konisch: RK = RCWA auf 1e-4, Energie auf 5e-4
  - [ ] s/p-Ausgabe gegen Buch Abb. IV.4 prüfen (optional)
- [x] **M5** 1D anisotrop: FFF_A, FFF_B, anisotrope Substrate/Schichten (4×4-Moden), Hauptachsen + Euler-Winkel
  - Abnahme: V7 (Kobalt), V8 (isotroper Tensor: Lamellen exakt, Profil 1e-3), einachsiges Substrat analytisch,
    biaxiale Schicht (Energie, Pfadgleichheit) – `tests/validation/test_v7_v8.cpp`
  - Offen: Literaturvergleich für Variante B
- [x] **M6** Nahfeld 1D (horizontale Linie, vertikale x–y-Karte), direkt/FFF-Rekonstruktion, Plot-Skript
  - Abnahme: V12 (1D) – `tests/validation/test_v12_nearfield.cpp`; Beispiel `examples/nearfield_lamellar_1d.json`
- [x] **M7** 2D isotrop: Polygon-Schichten, Normalenfeld (Algorithmus NF), RCWA; RK für Profile2D
  - Abnahme: V9 (exakt, 1e-12 / 1e-10), V10 qualitativ (Peak oberhalb des SPP-Pols, Polarisationsunabhängigkeit)
    – `tests/validation/test_v9_v10.cpp`, `tests/unit/test_normal_field2d.cpp`
  - Befund: Algorithmus NF als Tensorglättung (Plan v1.5); Profile2D zeilenweise exakt statt Marching Squares
  - V10 quantitativ (Buch 7,5 % bei 1,388 µm): Peak bis N = 8 unruhig, ab N = 9 bei 1,388 µm
    (N = 9: 4,2 %, N = 10: 6,2 %, N = 11: ≈ 9,2 %); Tensor- und Vektorvariante gleichwertig (N = 7, 8)
  - **Offen:** N = 11–15 auf dem Zielsystem (Konvergenz der Peakhöhe gegen ≈ 7,5 %?), Einfluss von σ
- [x] **M8** 2D anisotrop
  - Abnahme: V8 (2D) – `tests/validation/test_v8_2d.cpp`; Konvergenzstudie in `docs/validation.md`
  - `Polygon2D`: Laurent (exakt) und `FFF_B` (Rang-1-Feld aus Algorithmus NF, Faltung mit exakten
    Polygonkoeffizienten); `Profile2D`: `FFF_A` (Standard), `FFF_B`, Laurent
- [x] **M9** Nahfeld 2D (horizontal x–z, vertikal x–y)
  - Abnahme: V12 (2D) – `tests/validation/test_v12_2d.cpp`; Beispiel `examples/nearfield_pillars_2d.json`
  - FFF-konsistente Rekonstruktion in 2D: isotrope Polygone mit dem Tensorfeld (E = e_t + d_n/ε gilt für
    jedes P), Variante B/anisotrop mit dem Rang-1-Feld; Felder aus dem NF-Raster (bilinear)
  - Abnahme: V12 (2D)
- [~] **M10** CUDA-Backend (cuBLAS, cuSOLVER Xgetrf/Xgeev, cuFFT), Batching, Benchmarks
  - Abnahme: GPU = CPU auf 1e-10, Benchmark-Bericht in `docs/benchmarks.md`
  - Stand v0.9: `CudaBackend` (Host-API Zgemm, Xgetrf/Xgetrs, Xgeev), Dispatch der O(n³)-Stellen
    (`SolverOptions::backend`, Schwelle), paralleler Scan `solveScan` (OpenMP, ein Stream je Thread),
    `fgbench`; gegen CUDA 12.8 gebaut und gelinkt, ohne GPU laufen die CPU-Pfade (Tests
    `Backend.*`, GPU-Tests übersprungen)
  - Windows/MinGW + CUDA 13.4, RTX 3090: Build-Pfad `win-cuda` (experimentell), Backend-Tests laufen dort auf der GPU
  - Umgebungsvariablen `FG_BACKEND`, `FG_BACKEND_THRESHOLD` (ganze Testsuite auf der GPU)
  - **Offen:** Protokoll der GPU-Tests und `fgbench`-Tabelle in `docs/benchmarks.md`; Lauf der gesamten Suite mit
    `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0`; Schwelle festlegen; Abnahme auf dem Zielsystem (RTX 6000 Ada);
    cuFFT/eigene Kernel nur bei Bedarf
- [x] **M11** Parameter-Scans, Performance, Anwenderdoku
  - Scans: kartesisches Produkt λ × φ × ψ × δ (Zahl/Liste/Bereich), parallel (`solveScan`);
    Sweeps über beliebige Zahlenwerte per JSON-Pointer (`"sweep"`), CSV-Spalten `sweep_parameter`, `sweep_value`
    – `tests/unit/test_scan.cpp`, Beispiel `examples/sweep_lamellar_depth.json` (Smoke-Test)
  - `scripts/plot_efficiencies.py` (Spektren, Sweeps), Anwenderhandbuch `docs/user_guide.md` mit Laufzeittabelle
  - Performance: CPU-Referenzwerte in `docs/benchmarks.md`; GPU-Werte und Schwelle auf dem Zielsystem (M10)

### Optimierungen ohne Meilenstein (Plan v2.2, Optionen standardmäßig aus)

- [x] `"split_polarizations"`: 1D, γ = 0, nur isotrope Medien – `solve()` rechnet TE (Plan eq:TE) und TM (eq:TM)
  getrennt mit 2 n_f, Einbettung blockdiagonal, nicht angeregte Polarisation übersprungen; Nahfelder: volles System
  - Tests `PolarizationSplit.*` (Aufteilung = 4 n_f: RCWA 1e-10, RK 1e-8; konisch unverändert); Laufzeiten in
    `docs/benchmarks.md` (RCWA Faktor 1,9–5)
- [x] `"fft_sampling"`: FFT beliebiger Länge (Radix 2 / Bluestein) in `sampledCoefficients1D/2D`, `gridCoefficients2D`
  - Tests `Sampling.FFTEqualsDirectDFT`, `Sampling.FFTSwitchMatchesDFT`, `Sampling.SampledCoefficientsExactForTrigPolynomials`
- Standardoptionen: alle `examples/` bitgleich zu v0.10 geprüft (vor Anlage der Golden Files, Version 2)

## Version 2 (Plan v2.0, Abschnitte „Erweiterungen“ M12–M19)

Reihenfolge verbindlich (Abhängigkeiten im Plan, Abschnitt „Erweiterungen ab Version 2.0: Überblick“).
Vor dem ersten Eingriff: Golden Files aller `examples/` nach `tests/regression/golden/` (Stand v0.10), Test
`Regression.*` – mit Standardoptionen bleiben alle bisherigen Ergebnisse bitgleich.
Entscheidungen des Projektinhabers (5. 11. 2026, Plan „Querschnittsthemen“): iso = isolierte Strukturen (Steg, Graben,
Kontaktloch); Server AlmaLinux 10 + RTX 6000 Ada, Anwender Windows 11 + RTX 3090; CAD-Daten erzeugt ein Sprachmodell
mit Nurbforge; auf der RTX 3090 ist FP32 Standard.

- [x] Golden Files: `tests/regression/golden/` (alle Beispiele aus `tests/regression/examples.txt`, 17 Stellen),
  Test `Regression.*` (Label `regression`; Standard 1e-11 absolut auf den Effizienzen (Plan v2.4), `FG_REGRESSION_EXACT=1` bitgleich,
  `FG_REGRESSION_UPDATE=1` schreibt neu); erzeugt mit Standardoptionen auf dem Stand 5c037da (= v0.10, bitgleich geprüft)
- [x] **M12** Python-Anbindung (pybind11, scikit-build-core, Wheel manylinux_2_28, `fg.Simulation`)
  - [x] MSVC-Port: Presets `win-msvc`, `win-msvc-cuda` (`scripts\msvc.cmd`), Flags `/W4 /permissive- /utf-8`, ohne
    Warnungen; unit, validation, regression grün; `win-msvc-cuda` auf der RTX 3090: `Backend.*` (GPU = CPU) grün
  - [x] `fg::Simulation` (C++): Schichtmodelle über Rechnungen behalten, `setParameter` mit Teil-Neuaufbau
    (bitgleich zum Neuaufbau, `tests/unit/test_simulation.cpp`)
  - [x] Paket `fouriergrating` (`python/`, `pyproject.toml`, scikit-build-core): `solve`, `solve_scan`, `near_field`,
    `Simulation`, `set_backend`, `cuda_available`, `core_variant`; GIL frei; `FourierGratingError`; `py.typed`, `_core.pyi`
  - [x] Distribution `fouriergrating-cuda12` (`python/cuda12/`), Kernwahl beim Import; lokal (MSVC, CUDA 13.4,
    RTX 3090): `core_variant() == "cuda12"`, GPU = CPU auf 3e-15 (V10), pytest grün
  - [x] pytest `tests/python` (27 Tests, CPU- und CUDA-Kern grün): Beispiele = Golden Files (1e-11), Nahfeld =
    `fgsolve` (1e-12), `set_parameter` bitgleich, zwei Threads = sequentiell, GIL-Freigabe; `docs/python.md`,
    `examples/python/scan_lamellar.py`
  - [x] CI grün: `build-test`, `regression` (Linux-GCC 13 gegen Golden Files), `windows-msvc` (`/WX`), `python`;
    `wheels.yml`: manylinux_2_28 und Windows (CPython 3.10–3.13, Tests in frischer venv), `fouriergrating-cuda12`
    (CUDA 12.6.3)
  - [x] Anwenderrechner (Windows 11, RTX 3090): CI-Wheels in frischer venv mit `nvidia-*-cu12` (ohne `CUDA_PATH`):
    `core_variant() == "cuda12"`, V10 GPU = CPU auf 9e-15, pytest 27/27 – Protokoll in `docs/benchmarks.md`
  - Nicht geprüft: Installation des Linux-Wheels auf AlmaLinux 8–10 (nur im manylinux-Container getestet)
  - MSVC-Port (Wheels unter Windows nur mit MSVC; OpenMP `/openmp:llvm`), Windows-Wheels `fouriergrating` (CPU) und
    `fouriergrating-cuda12` (CUDA-Kern, nvidia-*-cu12 als Abhängigkeiten); Server: Build aus den Quellen
  - Abnahme: Python = `fgsolve` auf 1e-12 für alle Beispiele; GIL-Test; `set_parameter` = Neuaufbau; Wheels (Linux,
    Windows) aus CI installierbar; MSVC-CI grün; auf dem Windows-Rechner `fg.cuda_available()` und GPU = CPU;
    `docs/python.md`
- [x] **M13** Gekrümmte Querschnitte: Konturkoeffizienten (Plan eq:contour), `Contour2D`, Algorithmus NF mit
  Randstichproben und 3D-Normalen, `LayerModel` mit y-abhängigem Normalenfeld, `SectionedLayer` + `SectionProvider`
  - Abnahme: Unit-Tests (Polygon 1e-14, Kreis/Bessel 1e-13, NF bitgleich), V15 (Kreissäulen), V16a (Kegelstumpf
    gegen Profile2D 1e-6), Golden Files bitgleich
  - Erledigt: `geometry::BezierSegment/Contour2D/Region` (Kreis, Ellipse exakt, NURBS nach A5.6),
    `fourier::regionCoefficients` (adaptive Gauß-Quadratur, Matrixprodukte) – `tests/unit/test_contour.cpp`;
    Algorithmus NF aus Randstichproben (BVH statt Bucket-Gitter, 3D-Normalen, 3×3-Projektor; auf Polygonen bitgleich,
    512²/10⁴ in 0,47 s GCC / 0,68 s MSVC) – `test_normal_field_samples.cpp`; Schichttyp `contour2d` (JSON,
    `CrossSectionModel` gemeinsam mit `polygon2d`), Beispiel `examples/v15_circular_pillars.json` (Golden File);
    `LayerModel::normalAt/projectorAt/NN` mit y; `SectionedLayer` + `SectionProvider` (`fg/geometry/section.hpp`,
    Cache je Höhe)
  - V15: O(L⁻²) (Verhältnis 4,00), zweistufige Richardson ≤ 1,2e-10; RK = RCWA 1e-8 – `test_v15_contour.cpp`
  - V16a (Laurent, Plan v2.5): Kegelstumpf 5e-16, Kugelkappe 1e-15 gegen Profile2D, Energie 1e-10; FFF_sym mit
    3D-Normalen; Nahfeld folgt dem Querschnitt – `test_v16a_sectioned.cpp`. Golden Files bitgleich
- [x] **M14** CAD-Slicer mit Nurbforge
  - Vorleistung in Nurbforge v0.6.0 (eigener PR): `nf::section(Shell, Plane)`, `nf::critical_heights`,
    `nf::is_extrusion_band`, STEP-Namen lesen/schreiben, `StepLengthUnit::Micrometre`, Thread-Sicherheit
  - Geometrieprüfung mit Fehlercodes für die Prüfschleife des Sprachmodells (`fg.cad.check`, `llm_tool_spec`, `preview`)
  - FG: `FG_ENABLE_NURBFORGE`, Modul `fg::cad`, Achsen/Einheiten, Bandzerlegung, Materialzuordnung, Schichttyp `cad`,
    `--cad-report`, optionale Tschebyschow-Interpolation in y
  - Abnahme: V16b (Prisma/Chiralität), V17 (Zylinder), V18 (Kugelkappe, Kegelstumpf), zwei Materialien, Hohlraum;
    absichtlich fehlerhafte Dateien liefern den erwarteten Fehlercode; `docs/cad.md`, Beispiel `examples/cad_pillar_array.json`
  - Erledigt: Nurbforge v0.6.0 (PR #43, gemergt und getaggt; STEP-Leser war schon vorhanden): `section`, `critical_heights`,
    `is_extrusion_band`, `read_step_named`/`shell_names`, `Micrometre`, TSAN-Job. FG: `fg::cad` (`include/fg/cad`,
    `src/cad`), Schichttyp `cad` (expandiert beim Laden), `fgsolve --cad-report`, `fg.cad.check/preview/llm_tool_spec`,
    Tschebyschow-Interpolation der Q-Blöcke (`SectionedLayer::chebyshevNodes`), Presets `cpu-cad`/`win-cad`, CI-Job `cad`,
    Wheels mit Slicer, `scripts/make_cad_examples.py`, `tests/data/cad/`
  - Abnahme in `tests/cad/test_cad.cpp`: V16b (Koeffizienten 1e-13, RK/RCWA 1e-12, Spiegelung abgelehnt), V17 (Bessel
    1e-12), V18 (gegen Profile2D 1e-6, gegen M13-Provider 1e-9, Treppe konvergiert, Energie 1e-10, Tschebyschow 1e-9),
    zwei Materialien + Hohlraum (Bandfolge, Energie, RK ≈ RCWA), fehlerhafte Dateien → Codes; Golden File des Beispiels
  - Nicht umgesetzt: Reziprozitätstest (zwei Materialien), PML-Zone: mit M18 über `checkPmlZones` (CAD-Schnitte in den Zonen ⇒ Fehler)
- [x] **M15** Profiling und Parallelisierung
  - `fg::prof` und Profil vor M15 (Ryzen 9 5900X statt Zielsystem); Threadmodell (scan/slices/blas); parallele Scheiben-S;
    Sternprodukt als Baum fester Form; parallele Koeffizienten/NF/CAD-Schnitte; GPU-residente Ketten (`DeviceMatrix`)
  - Abnahme: Ergebnisse bitgleich für 1/2/4 Threads (alle Beispiele, sequentiell und Baum); Baum = sequentiell 1e-12;
    GPU-residente Kette = CPU 1e-10 (RTX 3090); TSAN-Job; Leistungsbericht (`docs/benchmarks.md`): Transferanteil 1,3 %
  - Offen: parallele Effizienz 36 % auf 12 Kernen (Ziel 60 % auf 16 Kernen; Eigenwertlöser bandbreitenbegrenzt) – auf dem
    Zielsystem nachmessen; RK-Stufen innerhalb einer Scheibe nicht parallel
- [x] **M16** Genauigkeit `fp64` | `mixed` | `fp32` am Backend (kein TF32), angepasstes Scheibenkriterium,
  Schutzmechanismus mit FP64-Rückfall, `CMAKE_CUDA_ARCHITECTURES` 86;89; Standard `auto`: FP32 auf der RTX 3090,
  FP64 auf dem Server und der CPU (`src/cuda/device_policy.cpp`, `FG_PRECISION`)
  - Abnahme: Backend-Unit-Tests (`test_precision`, CPU und RTX 3090), Genauigkeitsstudie (alle Beispiele + V9, V13;
    mixed/fp32 ≤ 2e-5 für V1–V9, Ziel festgeschrieben ≤ 1e-4), Laufzeiten RTX 3090 (N = 10/15/20), `AGENTS.md` angepasst
  - Offen: Laufzeiten RTX 6000 Ada (Zielsystem); V16–V18 in der Studie (CUDA-Build mit Slicer); Schichtmodell-Aufbau
    begrenzt die GPU-Laufzeit
- [x] **M17** Koordinatentransformationen (A = [[1/f']] α, G = [[1/g']] γ), ASR nach Granet für Metalle
  - Abnahme: η = 0 bitgleich; Matrixpfad mit f' = 1 = Diagonalpfad (1e-12); V1, V2 mit η = 0,5 auf 1e-8, V3, V9
    gemeinsamer Grenzwert und z-invariantes 2D = 1D (1e-10, Plan v2.9); V19 (Gold-Lamellen TM, Konvergenztabelle in
    `docs/validation.md`: 1e-3 ohne ASR ab N ≈ 60, mit η = 0,95 ab N ≤ 10)
- [x] **M18** PML (komplexe Streckung), isolierte Strukturen (Steg, Graben, Kontaktloch), strukturierte Ein-/Ausgangsmedien,
  Streufeldformulierung für schrägen Einfall (Plan eq:scat_source, eq:scat_ode)
  - Abnahme: V20a (Fresnel mit PML, 1e-10), V20b (Steg: Zellbreite, w, β < 1e-3), V20c (Wellenleiter: Durchgang 3e-6,
    abruptes Ende auf 1e-4 konvergiert), V20d (Streufeld = Totalfeld 1e-12, RK 1e-10), V20e (Normaleinfall 1e-10),
    V20f (Steg/Graben schräg: Gesamtstreuung 1e-4, Reziprozität 5e-2 bei L = 8 µm), V20g (Kontaktloch 2D: Symmetrie, CAD = Contour2D, 2D = 1D 1e-10);
    Plan v2.10, Messwerte in `docs/validation.md`
  - Offen: Literaturvergleiche (Silberstein 2001 für V20c, Hugonin/Lalanne 2005 für die PML),
    Konvergenzstudie V20g auf dem Zielsystem; Absorption verlustfreier Strukturen 2e-4 statt 1e-6
- [x] **M19** Differenzierbarkeit: diskrete Adjungierte (Plan „Konvention“, Rückwärtsregeln), Formableitungen
  (eq:shape2d, eq:shape3d), C++-API `gradient`/`vjp`, Python `value_and_grad`, PyTorch-/JAX-Adapter
  - Abnahme: Rückwärtsfunktionen gegen Differenzenquotienten 1e-6; Ende-zu-Ende V3, V9, V13, V15, V18 ≤ 7e-7
    (Laurent/Li ≤ 1e-7, Profil-Normalenfeld mitdifferenziert, Rasterfeld NF in v2.12 über eine Rasterzelle); Entartung 1e-9; RK = RCWA
    2e-9; Aufwand 1,3–2,3× Vorwärtsrechnung; Beispiele `optimize_littrow.py`, `optimize_nurbs_pillar.py`; Plan v2.11
  - Phase 2 umgesetzt: Wellenlänge, Winkel φ, ψ, δ, Brechzahl/Permittivität von Superstrat und Substrat (≤ 2e-7, RK 2e-6);
    `checkpoint_every` (bitgleich); anisotrope Materialien (Tensor, Hauptwerte, Euler-Winkel) mit Laurent/FFF_A/FFF_B in
    Lamellen und Profilen, Laurent in 2D-Querschnitten (≤ 3e-7); Koordinatentransformationen (ASR) mit Material-,
    Dicken- und Geometrieparametern (≤ 4e-7, RK 8e-6; Fixpunkte explizit für Geometrieparameter); anisotrope
    2D-Querschnitte mit FFF_B; Rasterfeld NF bei Formparametern mitbewegt (eingefroren war die Verschiebung eines
    Kreises um 50–150 % falsch, nur der Radius war geprüft); λ, Winkel, Superstrat/Substrat auch mit ASR (≤ 2,5e-7);
    `profile2d` isotrop und anisotrop (Materialien, Tiefe, Verschiebungen x0/z0; ≤ 3e-7)
    Aperiodisch (Plan v2.13, Stufe A): Streuflüsse oben/unten bei Normaleinfall, alle Parameterarten außer φ (≤ 3e-7);
    modaler Einfall, Leistungen der geführten Moden (≤ 1e-8); Streufeldformulierung (Stufe B, RCWA und RK): schräger
    Einfall, alle Parameterarten (≤ 1,4e-7)
  - Offen (Phase 2): GPU-Rückwärtsfunktionen; aperiodisch: seitlicher Fluss/Absorption (Stufe C, optional),
    Profilschichten mit PML
