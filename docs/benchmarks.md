# Benchmarks (M10)

Erzeugt mit `fgbench` (`apps/fgbench`):

```bash
fgbench --sizes 256,512,1024,1764 --orders 4,6,8,10,12 --repeat 2 > bench.md
```

Gemessen werden ZGEMM (n×n · n×n), LU-Lösen mit n rechten Seiten (`Xgetrf`/`Xgetrs` bzw. Eigen
`PartialPivLU`), das allgemeine Eigenwertproblem (`Xgeev` bzw. Eigen `ComplexEigenSolver`) und eine
vollständige RCWA-Rechnung des Ag-Lochgitters (V10, `Polygon2D`, `FFF_sym`) mit `backend = cpu` bzw. `cuda`.
GPU-Zeiten enthalten die Host↔Device-Transfers (das Backend arbeitet auf Eigen-Matrizen im Hostspeicher).

## Referenz: Sandbox (1 CPU-Kern, keine GPU), v0.9.0

| n | backend | ZGEMM n³ (s) | GFLOP/s | LU solve, n rhs (s) | eig (s) |
|---|---|---|---|---|---|
| 256 | cpu | 0,016 | 8,3 | 0,028 | 0,46 |
| 512 | cpu | 0,124 | 8,7 | 0,202 | 3,68 |
| 1024 | cpu | 1,07 | 8,0 | 1,56 | 33,8 |

RCWA, Ag-Lochgitter, eine Wellenlänge:

| N | 4 n_f | backend | Zeit (s) |
|---|---|---|---|
| 4 | 324 | cpu | 0,81 |
| 6 | 676 | cpu | 3,04 |
| 8 | 1156 | cpu | 13,5 |

Das Eigenwertproblem dominiert (≈ 30× ZGEMM gleicher Größe); bei N = 10 (reduziertes Problem 2 n_f = 882,
Schicht-S-Matrix 1764) sind auf einem Kern ≈ 45 s je Wellenlänge zu erwarten.

## Windows, RTX 3090 (sm_86), CUDA 13.4, MinGW/GCC – Messwerte **fehlen noch**

Mit dem Preset `win-cuda` gebaut (`docs/build_windows.md`); laut Commit laufen `Backend.CudaEqualsCpu` und
`Backend.CudaSolverEqualsCpu` dort gegen die GPU. Für die Abnahme hier einzutragen:

1. Ausgabe von `ctest --preset win-cuda -R Backend -V` (bestanden/übersprungen, Laufzeiten).
2. **Gesamte Validierung auf der GPU:** `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest --preset win-cuda`
   (PowerShell: `$env:FG_BACKEND="cuda"; $env:FG_BACKEND_THRESHOLD="0"`). Damit laufen alle O(n³)-Operationen
   aller Tests auf der GPU, nicht nur die beiden Vergleichstests (mit der Standardschwelle 1000 erreicht kaum ein
   Test die GPU). Erwartung: alle Tests grün.
3. `fgbench --sizes 256,512,1024,1764,2500 --orders 4,6,8,10,12 --repeat 2` (Tabelle unten einfügen).

Hinweis: Die RTX 3090 hat wie die RTX 6000 Ada nur 1/64 der FP32-Rate in FP64 (≈ 0,5 TFLOP/s); die Schwelle
`backend_threshold` daher erst nach dem Vergleich mit dem CPU-Wert (alle Kerne, `OMP_NUM_THREADS`) festlegen.

| n | backend | ZGEMM n³ (s) | GFLOP/s | LU solve (s) | eig (s) |
|---|---|---|---|---|---|
| | | | | | |

## Zielsystem: AlmaLinux 10, RTX 6000 Ada (sm_89) – **offen**

Auf dem Zielsystem auszuführen und hier einzutragen:

1. `cmake --preset cuda-release && cmake --build --preset cuda-release`
2. `ctest --preset cuda-release -R Backend` – `CudaEqualsCpu`, `CudaSolverEqualsCpu` müssen laufen (nicht
   übersprungen) und bestehen (GPU = CPU auf 1e-10); danach die gesamte Suite mit
   `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest --preset cuda-release`.
3. `build/cuda-release/apps/fgbench --sizes 256,512,1024,1764,2500 --orders 4,6,8,10,12 --repeat 2`
4. Schwelle `backend_threshold` (Standard 1000) anhand der Kreuzungspunkte CPU/GPU festlegen.
5. Parallelen Scan prüfen: `OMP_NUM_THREADS=1,2,4,8` mit `examples/v10_hole_array.json` (Scan über λ).

| n | backend | ZGEMM n³ (s) | GFLOP/s | LU solve (s) | eig (s) |
|---|---|---|---|---|---|
| | | | | | |

Hinweis zur FP64-Rate der RTX 6000 Ada (Plan: Größenordnung 1/64 der FP32-Rate): ZGEMM-GFLOP/s der GPU
mit dem CPU-Wert (alle Kerne) vergleichen, bevor die Schwelle gesenkt wird.

## Option `split_polarizations` (2n_f-Systeme, 1D in-plane isotrop), Plan v2.2

Windows 11, MSYS2 UCRT64 GCC (Preset `win-release`), ein Thread (`OMP_NUM_THREADS=1`), Gesamtlaufzeit von
`fgsolve` je Wellenlänge; `"split_polarizations": false` (volles 4n_f-System, Standard) gegen `true`.
Die Effizienzen stimmen überein (Tests `PolarizationSplit.*`: RCWA 1e-10, RK 1e-8).

| Fall | δ | 4 n_f (s) | 2 × 2 n_f (s) | Faktor |
|---|---|---|---|---|
| Al-Sinus (V4), RK, N = 20, tol 1e-5 | 90° (TM) | 20,6 | 13,9 | 1,5 |
| Al-Sinus (V4), RK, N = 20, tol 1e-5 | 45° | 20,8 | 21,2 | 1,0 |
| Al-Sinus (V4), RCWA 60 Scheiben, N = 20 | 90° (TM) | 1,56 | 0,59 | 2,7 |
| Al-Sinus (V4), RCWA 60 Scheiben, N = 20 | 45° | 1,57 | 0,84 | 1,9 |
| Lamellen (V3), RCWA, N = 60 | 90° (TM) | 0,48 | 0,09 | 5,1 |
| Lamellen (V3), RCWA, N = 60 | 45° | 0,48 | 0,15 | 3,1 |

Die Eigenprobleme waren schon vorher je Polarisation zerlegt (`layerModesAuto`, `rkTransferAuto`); der
Gewinn kommt aus den Schicht-S-Matrizen, Grenzflächen und Sternprodukten mit halber Größe und – bei reiner
TE- oder TM-Anregung – aus dem Überspringen der zweiten Polarisation. Bei RK mit gemischter Polarisation
dominiert die (schon zerlegte) Integration, daher kein Gewinn.

## M12: Abnahme der Wheels auf dem Anwenderrechner (Windows 11, RTX 3090), v0.11.0

CI-Artefakte von `wheels.yml` (Lauf 42026239722): `fouriergrating-0.11.0-cp311-win_amd64` und
`fouriergrating_cuda12-0.11.0-cp311-win_amd64` in einer frischen venv; Abhängigkeiten von PyPI
(`nvidia-cuda-runtime-cu12` 12.9.79, `nvidia-cublas-cu12` 12.9.2.10, `nvidia-cusolver-cu12` 11.7.5.82); `CUDA_PATH`
nicht gesetzt (die DLLs kommen nur aus den Pip-Paketen).

| Prüfung | Ergebnis |
|---|---|
| `fg.core_variant()`, `fg.cuda_available()` | `"cuda12"`, `True` |
| V10-Lochgitter (λ = 1,39 µm), `backend = cpu` gegen `cuda` (Schwelle 0) | max. Abweichung der Effizienzen 9,3·10⁻¹⁵ (T), 7,2·10⁻¹⁵ (R) |
| Laufzeit V10, eine Wellenlänge | CPU 16,9 s, GPU 5,3 s (Faktor 3,2) |
| `pytest tests/python` (mit `FG_FGSOLVE`) | 27 bestanden |

C++: `ctest --preset win-msvc-cuda -R Backend` (MSVC, CUDA 13.4) auf derselben GPU: `CudaEqualsCpu`,
`CudaSolverEqualsCpu` bestanden (GPU = CPU auf 1e-10).

## M15: Profil vor M15 und Leistungsbericht (v0.14.0)

Rechner: Windows 11, AMD Ryzen 9 5900X (12 Kerne, 24 Threads mit SMT, Dual-Channel DDR4, 128 GB), RTX 3090.
CPU: Preset `win-cad` (MinGW GCC 16, OpenMP); GPU: `win-msvc-cuda` (MSVC, CUDA 13.4). Konfigurationen in
`bench/` (je eine Wellenlänge), gemessen mit `fgsolve --profile` (Zeiten je Abschnitt summiert über alle Threads;
„wall“ = Laufzeit). Threadzahl über `FG_NUM_THREADS`.

### Profil vor M15 (ein Thread)

| Fall | wall | Hauptanteile |
|---|---|---|
| V3 Lamellen 1D, RCWA, N = 100 | 2,4 s | Scheiben-S 0,79 s, Sternprodukt 0,51 s, Eigenwerte 0,35 s |
| V9 `sinus_x` im 2D-Code, RK, N = 8 × 1 | 73 s | RK-Transfer (2 Stabilitätsscheiben) 73 s, davon Q(y) 9 s |
| V13 Gold-Sinus TM, RK, N = 16 | 32 s | RK-Transfer (10 Scheiben) 32 s, davon Q(y) 14,5 s, M 5,3 s |
| Eierkarton 2D, RCWA 32 Scheiben, N = 7 | 354 s | Eigenwerte 255 s, Scheiben-S 50 s, Sternprodukt 32 s |
| CAD-Säulen (Beispiel), RCWA 30 + 1 Scheiben, N = 8 | 881 s¹ | Eigenwerte 664 s, Scheiben-S 117 s, Sternprodukt 65 s, Schnitte + NF 30 s |
| V10 Lochgitter, N = 10 (eine Schicht) | 39 s | Scheiben-S 12 s, Sternprodukt 10 s, Eigenwerte 7 s, T→S 3 s |

¹ teilweise parallel zu anderen Läufen gemessen (zu hoch).

Priorisierung danach: Scheiben parallel (Eigenwerte dominieren 2D), Sternprodukte als Baum, BLAS-Ebene für
Einzelschichten, GPU-residente Kette. Koeffizienten, Normalenfeld und CAD-Schnitte liegen unter 5 % (trotzdem
parallel, wenn der Aufrufer sequentiell ist).

### Nach M15 (CPU)

| Fall | 1 Thread | 12 Threads | 16 Threads | 24 Threads |
|---|---|---|---|---|
| Eierkarton 2D, 32 Scheiben, `tree` | 354 s | 82,5 s | 74 s | 82 s |
| Eierkarton 2D, 32 Scheiben, `sequential` | – | – | 91 s | 108 s |
| CAD-Säulen, `tree` | (881 s) | 220 s | 188 s | 208 s |
| V13 Gold-Sinus RK (10 Scheiben) | 32 s | – | 5,0 s | 4,9 s |
| V9 RK (2 Scheiben) | 73 s | – | 52 s | 48 s |
| V10, `threads.blas = 24` (Eigen-GEMM parallel) | 39 s | – | 21 s | 19 s |

Parallele Effizienz (2D, 32 Scheiben, Baum): 4,3× auf 12 Kernen (36 %), 4,8× mit 16 Threads – unter dem Ziel
60 % auf 16 Kernen. Ursache: Eigens `ComplexEigenSolver` (Hessenberg/Schur, überwiegend Level-2) ist auf diesem
Rechner speicherbandbreitenbegrenzt; unter Last dauert ein Eigenproblem 2× (12 Threads) bis 4× (24 Threads) so
lange wie allein. Dazu kommt die Granularität (32 Scheiben auf 12 Threads = 3 Runden). Die sequentielle Kette
begrenzt zusätzlich (Sternprodukte ≈ 10 % der Arbeit); `"reduction": "tree"` gewinnt hier 15–25 %. Auf dem
Zielsystem (mehr Speicherkanäle) erneut zu messen. V9 (2 Stabilitätsscheiben) profitiert kaum; die parallele
Auswertung der RK-Stufen innerhalb einer Scheibe (Plan Schritt 3) ist nicht umgesetzt.

### GPU-residente Kette (RTX 3090, FP64)

| Fall | CPU 1 Thread | CPU bester Wert | GPU 1 Stream | GPU 4 Streams |
|---|---|---|---|---|
| Eierkarton 2D, 32 Scheiben, N = 7 | 354 s | 74 s | 128 s | 48 s |
| V10 Lochgitter, N = 10 | 39 s | 19 s | 11,0 s | 10,7 s |

Transferanteil (Nsight Systems 2026.3, V10): Host→Device 0,7 %, Device→Host ≈ 0,5 %, Device→Device 0,1 % der
GPU-Zeit (Ziel < 10 %). Die GPU-Zeit dominieren `Xgeev` (Bulge-Chasing) und ZGEMM. Auf der RTX 3090 hat FP64 nur 1/64
der FP32-Rate; der Gewinn durch FP32 folgt mit M16. GPU-Ergebnisse = CPU auf 1e-10 (`Resident.*`, `Backend.*`).
