# Build auf AlmaLinux 10

## 1. Pakete

```bash
sudo ./scripts/setup_almalinux10.sh              # Compiler, CMake, Eigen, JSON, GTest, Python, LaTeX
sudo ./scripts/setup_almalinux10.sh --with-cuda  # zusätzlich CUDA-Toolkit aus dem NVIDIA-Repo rhel10
```

Das Skript aktiviert EPEL und CRB und installiert mit `--skip-unavailable`. Fehlt ein Paket
(z. B. `eigen3-devel`, `json-devel`, `gtest-devel`), lädt CMake die Bibliothek per FetchContent
nach (Internetzugang beim ersten Konfigurieren nötig).

AlmaLinux 10 verwendet das NVIDIA-Repository `rhel10/x86_64` (Standard-GCC 14). Der
Grafiktreiber wird **nicht** vom Skript installiert – bitte nach dem NVIDIA Driver Installation
Guide (Abschnitt RHEL/Alma 10) vorgehen und mit `nvidia-smi` prüfen.

Mindestversionen: CMake 3.24, GCC 13+ (C++20), CUDA 12.6 Update 2 (wegen `cusolverDnXgeev`).

## 2. Konfigurieren, bauen, testen

```bash
cmake --preset cpu-release && cmake --build --preset cpu-release && ctest --preset cpu-release
cmake --preset cuda-release && cmake --build --preset cuda-release && ctest --preset cuda-release
./build/cuda-release/apps/fgsolve --info     # zeigt "cuda:NVIDIA RTX 6000 Ada Generation (sm_89, …)"
```

Optionen: `FG_ENABLE_CUDA`, `FG_ENABLE_OPENMP` (Standard an), `FG_ENABLE_HDF5`,
`FG_BUILD_TESTS`, `FG_WARNINGS_AS_ERRORS`. Die CUDA-Architektur ist auf `89` voreingestellt.

## 3. Hinweise zur GPU

Die RTX 6000 Ada hat eine im Verhältnis zu FP32 geringe FP64-Leistung. Der Kern rechnet in
doppelter Genauigkeit; die GPU lohnt sich daher vor allem bei 2D-Problemen (Matrizen ≥ 1000)
und bei gebatchten Wellenlängen-Scans. Vor M10 einen ZGEMM-/ZGETRF-Benchmark laufen lassen und
die Schwelle in `fg::linalg::makeBackend` anpassen.

## 4. Dokumentation

```bash
./scripts/build_docs.sh      # erzeugt docs/implementation_plan/main.pdf (pdflatex, 3 Läufe)
```
