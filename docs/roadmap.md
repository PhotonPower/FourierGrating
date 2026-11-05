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

## Version 2 (Plan v2.0, Abschnitte „Erweiterungen“ M12–M19)

Reihenfolge verbindlich (Abhängigkeiten im Plan, Abschnitt „Erweiterungen ab Version 2.0: Überblick“).
Vor dem ersten Eingriff: Golden Files aller `examples/` nach `tests/regression/golden/` (Stand v0.10), Test
`Regression.*` – mit Standardoptionen bleiben alle bisherigen Ergebnisse bitgleich.
Entscheidungen des Projektinhabers (5. 11. 2026, Plan „Querschnittsthemen“): iso = isolierte Strukturen (Steg, Graben,
Kontaktloch); Server AlmaLinux 10 + RTX 6000 Ada, Anwender Windows 11 + RTX 3090; CAD-Daten erzeugt ein Sprachmodell
mit Nurbforge; auf der RTX 3090 ist FP32 Standard.

- [ ] **M12** Python-Anbindung (pybind11, scikit-build-core, Wheel manylinux_2_28, `fg.Simulation`)
  - MSVC-Port (Wheels unter Windows nur mit MSVC; OpenMP `/openmp:llvm`), Windows-Wheels `fouriergrating` (CPU) und
    `fouriergrating-cuda12` (CUDA-Kern, nvidia-*-cu12 als Abhängigkeiten); Server: Build aus den Quellen
  - Abnahme: Python = `fgsolve` auf 1e-12 für alle Beispiele; GIL-Test; `set_parameter` = Neuaufbau; Wheels (Linux,
    Windows) aus CI installierbar; MSVC-CI grün; auf dem Windows-Rechner `fg.cuda_available()` und GPU = CPU;
    `docs/python.md`
- [ ] **M13** Gekrümmte Querschnitte: Konturkoeffizienten (Plan eq:contour), `Contour2D`, Algorithmus NF mit
  Randstichproben und 3D-Normalen, `LayerModel` mit y-abhängigem Normalenfeld, `SectionedLayer` + `SectionProvider`
  - Abnahme: Unit-Tests (Polygon 1e-14, Kreis/Bessel 1e-13, NF bitgleich), V15 (Kreissäulen), V16a (Kegelstumpf
    gegen Profile2D 1e-6), Golden Files bitgleich
- [ ] **M14** CAD-Slicer mit Nurbforge
  - Vorleistung in Nurbforge v0.6.0 (eigener PR): `nf::section(Shell, Plane)`, `nf::critical_heights`,
    `nf::is_extrusion_band`, STEP-Namen lesen/schreiben, `StepLengthUnit::Micrometre`, Thread-Sicherheit
  - Geometrieprüfung mit Fehlercodes für die Prüfschleife des Sprachmodells (`fg.cad.check`, `llm_tool_spec`, `preview`)
  - FG: `FG_ENABLE_NURBFORGE`, Modul `fg::cad`, Achsen/Einheiten, Bandzerlegung, Materialzuordnung, Schichttyp `cad`,
    `--cad-report`, optionale Tschebyschow-Interpolation in y
  - Abnahme: V16b (Prisma/Chiralität), V17 (Zylinder), V18 (Kugelkappe, Kegelstumpf), zwei Materialien, Hohlraum;
    absichtlich fehlerhafte Dateien liefern den erwarteten Fehlercode; `docs/cad.md`, Beispiel `examples/cad_pillar_array.json`
- [ ] **M15** Profiling und Parallelisierung
  - `fg::prof` und Profil vor M15 (Zielsystem); Threadmodell (scan/slices/blas); parallele Scheiben-S;
    Sternprodukt als Baum fester Form; parallele Koeffizienten/NF/CAD-Schnitte; GPU-residente Ketten (`DeviceMatrix`)
  - Abnahme: Ergebnisse bitgleich für 1/2/4 Threads; Baum = sequentiell 1e-12; TSAN-Job; Leistungsbericht
- [ ] **M16** Genauigkeit `fp64` | `mixed` | `fp32` am Backend (kein TF32), angepasstes Scheibenkriterium,
  Schutzmechanismus mit FP64-Rückfall, `CMAKE_CUDA_ARCHITECTURES` 86;89; Standard `auto`: FP32 auf der RTX 3090,
  FP64 auf dem Server und der CPU (`src/cuda/device_policy.cpp`, `FG_PRECISION`)
  - Abnahme: Backend-Unit-Tests, Genauigkeitsstudie V1–V18 (Tabelle in `docs/validation.md`), Laufzeiten RTX 3090 /
    RTX 6000 Ada; Regel in `AGENTS.md` angepasst
- [ ] **M17** Koordinatentransformationen (A = [[1/f']] α, G = [[1/g']] γ), ASR nach Granet für Metalle
  - Abnahme: η = 0 bitgleich; V1–V3, V9 mit η = 0,5 auf 1e-8; V19 (Metall-Lamellen TM, Konvergenztabelle)
- [ ] **M18** PML (komplexe Streckung), isolierte Strukturen (Steg, Graben, Kontaktloch), strukturierte Ein-/Ausgangsmedien,
  Streufeldformulierung für schrägen Einfall (Plan eq:scat_source, eq:scat_ode)
  - Abnahme: V20a (Fresnel mit PML), V20b (isolierter Steg), V20c (Wellenleiterende, Literatur), V20d (Streufeld =
    Totalfeld im periodischen Fall, ohne PML), V20e (Normaleinfall), V20f (Graben/Steg schräg), V20g (Kontaktloch 2D)
- [ ] **M19** Differenzierbarkeit: diskrete Adjungierte (Plan „Konvention“, Rückwärtsregeln), Formableitungen
  (eq:shape2d, eq:shape3d), C++-API `gradient`/`vjp`, Python `value_and_grad`, PyTorch-/JAX-Adapter
  - Abnahme: Rückwärtsfunktionen gegen Differenzenquotienten 1e-6; Ende-zu-Ende V3, V9, V13, V15, V18; Entartung;
    Aufwand ≤ 5× Vorwärtsrechnung; Beispiele `optimize_littrow.py`, `optimize_nurbs_pillar.py`
