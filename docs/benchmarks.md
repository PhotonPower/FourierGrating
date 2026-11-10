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
   aller Tests auf der GPU, nicht nur die beiden Vergleichstests (mit der früheren Standardschwelle 1000 erreichte kaum ein
   Test die GPU). Erwartung: alle Tests grün.
3. `fgbench --sizes 256,512,1024,1764,2500 --orders 4,6,8,10,12 --repeat 2` (Tabelle unten einfügen).

Hinweis: Die RTX 3090 hat wie die RTX 6000 Ada nur 1/64 der FP32-Rate in FP64 (≈ 0,5 TFLOP/s); die Schwelle
`backend_threshold` daher erst nach dem Vergleich mit dem CPU-Wert (alle Kerne, `OMP_NUM_THREADS`) festlegen.

| n | backend | ZGEMM n³ (s) | GFLOP/s | LU solve (s) | eig (s) |
|---|---|---|---|---|---|
| | | | | | |

## Zielsystem: AlmaLinux 10, RTX 6000 Ada (sm_89, 48 GB), CUDA 12.9 – Abnahme (v0.18.0)

Rechner: 32 logische CPU-Kerne, GCC 14, Preset `cuda-release`. CPU-Werte: Eigen, ein Thread (Eigen parallelisiert
Eigenwertproblem und LU nicht); GPU-Zeiten einschließlich Host↔Device-Transfers.

**Tests.** `ctest --preset cuda-release -R Backend`: 6/6 bestanden, `CudaEqualsCpu` und `CudaSolverEqualsCpu` laufen
(nicht übersprungen). Gesamte Suite mit `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0` (280 Tests, 95 min): im ersten Lauf
20 Fehler, davon 17 durch `cusolverDnXgeev` (Status 7, `CUSOLVER_STATUS_INTERNAL_ERROR`) bei gleichzeitigem Aufruf aus
mehreren Threads (nur mit `FG_NUM_THREADS` ≥ 2 und `FG_GPU_STREAMS` ≥ 2; jeder Thread hat eigenes Handle, eigenen
Stream und eigene Workspaces). Behoben durch einen Mutex um Xgeev (`CudaBackend::eigT`); Regressionstest
`Backend.CudaConcurrentEigenSolves` (schlägt ohne Fix fehl). Wiederholung der 20 Tests: 18 bestanden. Übrig:
* `Backend.DispatchSettings`: setzt die Einstellungen selbst und kollidiert mit den Umgebungsvariablen; ohne sie bestanden.
* `NormalFieldSamples.PolygonsBitIdenticalTensor/Vector`: kein GPU-Thema (reine Rasterrechnung). Ursache: GCC 14 auf
  AlmaLinux 10 zielt auf x86-64-v3 (FMA) und fusst mit `-ffp-contract=fast` a*b+c in den beiden Kopien von Algorithmus NF
  verschieden (156 FMA-Instruktionen in `normal_field.cpp.o`, Abweichungen ≤ 5e-15 in ca. 20 % der Punkte, Regionen gleich).
  Behoben durch `-ffp-contract=off` für `geometry/normal_field.cpp` (`src/CMakeLists.txt`).
* Übersprungen: `cad_pillar_array` (Nurbforge nicht gebaut).

**Mikrobenchmark** (`fgbench --sizes 256,512,1024,1764,2500 --orders 4,6,8,10,12 --repeat 2`):

| n | backend | ZGEMM n^3 (s) | GFLOP/s | LU solve, n rhs (s) | eig (s) |
|---|---|---|---|---|---|
| 256 | cpu | 0.010 | 13.6 | 0.055 | 1.217 |
| 256 | cuda | 0.002 | 88.0 | 0.008 | 0.476 |
| 512 | cpu | 0.056 | 19.2 | 0.513 | 21.192 |
| 512 | cuda | 0.006 | 184.3 | 0.028 | 1.161 |
| 1024 | cpu | 0.210 | 40.9 | 3.067 | 193.360 |
| 1024 | cuda | 0.037 | 229.1 | 0.087 | 3.661 |
| 1764 | cpu | 1.740 | 25.2 | 10.121 | 360.237 |
| 1764 | cuda | 0.187 | 234.3 | 0.334 | 9.208 |
| 2500 | cpu | 3.376 | 37.0 | 26.697 | 1026.417 |
| 2500 | cuda | 0.503 | 248.4 | 1.044 | 22.103 |

RCWA, Ag hole array (V10), lambda = 1.388 um, one solve:

| N | 4 n_f | backend | time (s) |
|---|---|---|---|
| 4 | 324 | cpu | 0.89 |
| 4 | 324 | cuda | 0.58 |
| 6 | 676 | cpu | 5.90 |
| 6 | 676 | cuda | 2.31 |
| 8 | 1156 | cpu | 22.52 |
| 8 | 1156 | cuda | 3.32 |
| 10 | 1764 | cpu | 77.29 |
| 10 | 1764 | cuda | 8.43 |
| 12 | 2500 | cpu | 216.73 |
| 12 | 2500 | cuda | 21.46 |

ZGEMM erreicht ≈ 250 GFLOP/s (FP64; Ada hat 1/64 der FP32-Rate). Die Eigenwerte dominieren: GPU 17–46× schneller
als Eigen (n = 1764: 9,2 s gegen 360 s). Gesamte V10-Rechnung je Wellenlänge bei N = 10: 8,4 s gegen 77 s (9×), N = 12: 21 s gegen 217 s.

**Kreuzungspunkte** (`fgbench --sizes 32,64,96,128,192 --repeat 5`): ZGEMM und LU ab n ≈ 64–96, Eigenwertproblem ab
n ≈ 128–192 (n = 128: 0,167 s CPU, 0,178 s GPU; n = 192: 0,54 s gegen 0,36 s). Daraus folgt die Standardschwelle
`backend_threshold` = 160 (vorher 1000; gilt für `"backend": "auto"`, Matrixgröße n der jeweiligen Operation). Rechnungen
mit Matrizen zwischen 160 und 1000 laufen damit auf der GPU und weichen von der CPU nur im Rundungsbereich ab.

**Paralleler Scan** (V10, N = 10, 8 Wellenlängen, `backend = cuda`, Schwelle 0): 1 Thread 88 s, 2 Threads 50 s,
4 Threads 35 s, 8 Threads 36 s (Spitzenspeicher des Prozesses 1,2 / 1,8 / 2,9 / 5,3 GB). Mehr als 4 Threads bringen
nichts: die Eigenwertprobleme sind serialisiert (Mutex) und die GPU ist ausgelastet. Zum Vergleich CPU mit 24 Threads: 82 s.
Effizienzen GPU gegen CPU: max. Abweichung 1·10⁻¹⁵ (24 Werte).

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

## M16: Laufzeiten nach Genauigkeit (RTX 3090, v0.15.0)

Ag-Lochgitter (V10-Geometrie) mit N = Nx = Nz, eine Wellenlänge, `backend = cuda`, Schwelle 0, ein Stream;
`win-msvc-cuda`. In Klammern: `solve.assemble` (GPU-Teil inkl. Aufbau von M) und `solve.layer_models` (FFF-Koeffizienten,
FP64 auf der CPU, ein Thread). T(0,0) = Transmission nullter Ordnung.

| N | 4 n_f | fp64 | mixed | fp32 | T(0,0) fp64 | Abw. fp32 |
|---|---|---|---|---|---|---|
| 10 | 1764 | 6,9 s (5,4 / 1,3) | 8,0 s | 3,5 s (2,1 / 1,2) | 0,0535779 | 2,0e-6 |
| 15 | 3844 | 48,6 s (34,0 / 13,6) | 52,7 s | 46,8 s (12,7 / 33,5¹) | 0,0237220 | 2,6e-5 |
| 20 | 6724 | 231 s (157 / 72) | 340 s | 157 s (57 / 98¹) | 0,0111617 | 3,4e-5 |

¹ Der Aufbau der Schichtmodelle schwankt von Lauf zu Lauf (13–33 s bei N = 15); er läuft auf der CPU in FP64.

FP32 beschleunigt den GPU-Teil um 2,6–2,7×; die Gesamtzeit begrenzt jetzt der Aufbau der Schichtmodelle
(Toeplitz-Inverse der FFF-Faktorisierung, n_f × n_f, FP64, ein Thread) – nächster Ansatzpunkt (BLAS-Threads oder GPU in
FP64). `mixed` ist auf der GPU langsamer als `fp64` (Host-API mit Transfers je Operation) und nicht genauer als
`fp32` (Studie in `docs/validation.md`). Vor M16 liefen die Produkte Ψ⁻¹Ψ als dichte FP64-GEMM auf der CPU: fp64 N = 20
518 s, fp32 379 s. RTX 6000 Ada: auf dem Zielsystem zu messen.
