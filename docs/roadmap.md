# Roadmap

Abnahme = alle genannten Tests grün (`ctest`), Details im Plan, Abschnitt „Meilensteine“.

- [x] **M0** Build-System, Abhängigkeiten, CI (CPU), Grundgerüst, Doku
  - Abnahme: `cmake --build` + `ctest` grün (CI, AlmaLinux 10 noch zu bestätigen)
- [ ] **M1** Kernbausteine
  - [x] IndexMap, Fourierkoeffizienten, Toeplitz, Butcher + RK-Propagator, S-Matrix, isotrope Moden, `M`, Einfall
  - [ ] FFT-Abtastung stetiger Funktionen (PocketFFT/FFTW)
  - [ ] Effizienzen aus Poynting-Fluss (Plan eq:eff), s/p-Zerlegung
- [ ] **M2** 1D in-plane isotrop: homogene Schichten, Lamellen (RCWA; Laurent, LiLamellar), Gesamtlösung, JSON-Konfiguration
  - Abnahme: V1, V2, V3 (TE/TM)
- [ ] **M3** 1D-Profile: Nullstellen-Schnitte, Normalenfeld, FFF_iso; RK-Löser mit Schrittsteuerung; RCWA-Treppe mit FFF
  - Abnahme: V4, V5, V6; RK = RCWA auf Lamellen
- [ ] **M4** 1D konisch (4n_f-System)
  - Abnahme: V1–V3 konisch, Reziprozität
- [ ] **M5** 1D anisotrop: FFF_A, FFF_B, anisotrope Substrate (4×4-Moden)
  - Abnahme: V7, V8
- [ ] **M6** Nahfeld 1D (horizontale Linie, vertikale x–y-Karte), direkt/FFF-Rekonstruktion, Plot-Skript
  - Abnahme: V12 (1D)
- [ ] **M7** 2D isotrop: Polygon-Schichten, Normalenfeld (Algorithmus NF), RCWA; RK für Profile2D
  - Abnahme: V9, V10
- [ ] **M8** 2D anisotrop
  - Abnahme: V8 (2D), Konvergenzstudie
- [ ] **M9** Nahfeld 2D (horizontal x–z, vertikal x–y)
  - Abnahme: V12 (2D)
- [ ] **M10** CUDA-Backend (cuBLAS, cuSOLVER Xgetrf/Xgeev, cuFFT), Batching, Benchmarks
  - Abnahme: GPU = CPU auf 1e-10, Benchmark-Bericht in `docs/benchmarks.md`
- [ ] **M11** Parameter-Scans, Performance, Anwenderdoku
