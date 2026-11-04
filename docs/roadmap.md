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
  - **Offen:** V10 quantitativ (Buch 7,5 % bei 1,388 µm): Peak wandert bis N = 8 nach Rot, Höhe schwankt;
    Konvergenzstudie N = 10–15 auf dem Zielsystem (M10/M11), Vergleich Tensor/Vektor/Laurent
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
  - **Offen (Zielsystem):** `ctest -R Backend` mit GPU, Benchmark-Tabelle, Schwelle festlegen; cuFFT/eigene
    Kernel nur bei Bedarf (M11)
- [ ] **M11** Parameter-Scans, Performance, Anwenderdoku
