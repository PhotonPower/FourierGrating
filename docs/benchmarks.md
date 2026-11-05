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
