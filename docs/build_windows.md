# Build auf Windows

Inoffizieller Pfad für Windows-Entwickler (MSYS2/MinGW-w64, UCRT64). Zielsystem laut
[`AGENTS.md`](../AGENTS.md), Abschnitt 7, bleibt AlmaLinux 10/GCC 14 (CI baut nur CPU dort); der
CPU-Referenzpfad (Preset `win-release`) muss aber auch hier immer funktionieren. CUDA
(Preset `win-cuda`, Abschnitt 3) ist experimentell und nicht Teil der CI.

## 1. Pakete

```powershell
.\scripts\setup_windows.ps1              # installiert MSYS2 (falls nötig) und das UCRT64-Toolchain
.\scripts\setup_windows.ps1 -SkipMsys2Install   # MSYS2 ist schon unter C:\msys64 installiert
```

Das Skript installiert über `pacman` in der `ucrt64`-Umgebung: `gcc`, `cmake`, `ninja`, `pkgconf`,
`eigen3`, `nlohmann-json`, `gtest`, sowie `python` mit `numpy`/`matplotlib` für die Plot-Skripte
(`scripts/plot_efficiencies.py`, `scripts/plot_nearfield.py`). Fehlt ein Paket, lädt CMake die
Bibliothek per `FetchContent` nach (Internetzugang beim ersten Konfigurieren nötig).

Mindestversionen: CMake 3.24, GCC 13+ (C++20). MSYS2 liefert zum Zeitpunkt dieser Anleitung
GCC 16, das reicht komfortabel.

## 2. Konfigurieren, bauen, testen

```powershell
cmake --preset win-release
cmake --build --preset win-release
ctest --preset win-release
```

Das Preset `win-release` setzt Compiler (`C:\msys64\ucrt64\bin\gcc.exe`/`g++.exe`) und `PATH`
(für die mingw-Laufzeit-DLLs wie `libstdc++-6.dll`) automatisch, unabhängig vom System-`PATH`.
Läuft eine gebaute Binary außerhalb von `cmake`/`ctest` (z. B. direkt aus einem neuen
PowerShell-Fenster), muss `C:\msys64\ucrt64\bin` im `PATH` stehen:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;" + $env:PATH
.\build\win-release\apps\fgsolve.exe --info
.\build\win-release\apps\fgsolve.exe examples\v3_lamellar_dielectric.json   # schreibt out/v3/efficiencies.csv
```

## 3. CUDA (experimentell)

```powershell
cmake --preset win-cuda
cmake --build --preset win-cuda
ctest --preset win-cuda
```

NVIDIAs Windows-Importbibliotheken (`cudart.lib`, `cusolver.lib`, ...) sind MSVC-gebaut und
referenzieren MSVC-CRT-Symbole für /GS-Stack-Cookies und SEH (`__security_check_cookie`,
`__GSHandlerCheck`, `__security_cookie`), die MinGW nicht mitbringt — Linken schlägt sonst mit
`undefined reference` fehl. `cmake/mingw_cuda_security_stub.cpp` stellt diese drei Symbole als
No-ops bereit (sicher: MinGW nutzt die zugehörigen Codepfade nie); `CMakeLists.txt` hängt sie bei
`FG_ENABLE_CUDA AND MINGW` per `link_libraries()` direkt an jede Executable (nicht an `fg_core`,
da ein über eine statische Bibliothek verlinktes Objekt zu einem Archiv-Member würde, das vor
Erreichen von `cudart.lib` in der Link-Zeile schon abgehakt ist und so nie gezogen wird).

Getestet mit CUDA Toolkit 13.4 und einer RTX 3090 (sm_86) — abweichend vom spezifizierten
Zielsystem (RTX 6000 Ada, sm_89, AlmaLinux 10). `Backend.CudaEqualsCpu`/`CudaSolverEqualsCpu`
laufen dabei gegen die echte GPU statt übersprungen zu werden.

## 4. MSVC 2022 (M12: Python-Wheels, offizieller Windows-Pfad)

Python unter Windows ist mit MSVC gebaut; die Wheels (M12) und die CUDA-Host-API brauchen daher MSVC.
Voraussetzung: Visual Studio 2022 oder die Build Tools mit „Desktopentwicklung mit C++“ (enthält CMake und
Ninja). `scripts\msvc.cmd` führt einen Befehl in der `vcvars64`-Umgebung aus:

```powershell
scripts\msvc.cmd cmake --preset win-msvc
scripts\msvc.cmd cmake --build --preset win-msvc
scripts\msvc.cmd ctest --preset win-msvc -L unit
scripts\msvc.cmd ctest --preset win-msvc -L validation
scripts\msvc.cmd ctest --preset win-msvc -L regression   # cmd: kein "a|b" in -L (Pipe)
```

Unter MSVC ignoriert `cmake/FGDependencies.cmake` die MSYS2-Präfixe und die CMake-Paketregistry, sonst würden
MinGW-Header bzw. Eigen-Kopien anderer Projekte gefunden; Eigen, nlohmann_json und GoogleTest kommen per
`FetchContent`. Flags: `/W4 /permissive- /utf-8 /bigobj`, externe Header (`<...>`) mit `/external:W0`,
`_USE_MATH_DEFINES`, `NOMINMAX`. OpenMP: MSVC-Standard (`-openmp`, OpenMP 2.0); für Tasks (M15) `/openmp:llvm`.
`win-msvc-cuda` baut zusätzlich das CUDA-Backend (`86;89`) und ersetzt den experimentellen MinGW-CUDA-Pfad.

Unterschiede zu GCC liegen bei ≤ 1,2·10⁻¹⁴ in den Effizienzen (Test `Regression.*`).

## 5. Einschränkungen

* **Kein LaTeX-Setup.** `./scripts/build_docs.sh` (PDF des Plans) ist für Linux/`pdflatex`
  geschrieben; unter Windows z. B. über eine MiKTeX-Installation und Anpassung des Skripts möglich,
  aber nicht Teil dieser Anleitung.
* Dieser Pfad ist nicht Teil der CI; vor einem Commit bleibt `ctest --preset cpu-release` auf dem
  Linux-Zielsystem (oder äquivalent hier `ctest --preset win-release`) maßgeblich (siehe
  [`AGENTS.md`](../AGENTS.md), Abschnitt 3).
