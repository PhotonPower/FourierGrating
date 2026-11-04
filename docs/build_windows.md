# Build auf Windows

Inoffizieller CPU-Pfad für Windows-Entwickler (MSYS2/MinGW-w64, UCRT64). Zielsystem laut
[`AGENTS.md`](../AGENTS.md), Abschnitt 7, bleibt AlmaLinux 10/GCC 14 (CI baut nur CPU dort); dieser
Pfad verzichtet auf CUDA, der CPU-Referenzpfad muss aber auch hier immer funktionieren.

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

## 3. Einschränkungen

* **Kein CUDA-Backend.** `FG_ENABLE_CUDA` ist im Preset `win-release` fest auf `OFF`; für GPU-Läufe
  (M10) siehe [`docs/build_almalinux10.md`](build_almalinux10.md).
* **Kein LaTeX-Setup.** `./scripts/build_docs.sh` (PDF des Plans) ist für Linux/`pdflatex`
  geschrieben; unter Windows z. B. über eine MiKTeX-Installation und Anpassung des Skripts möglich,
  aber nicht Teil dieser Anleitung.
* Dieser Pfad ist nicht Teil der CI; vor einem Commit bleibt `ctest --preset cpu-release` auf dem
  Linux-Zielsystem (oder äquivalent hier `ctest --preset win-release`) maßgeblich (siehe
  [`AGENTS.md`](../AGENTS.md), Abschnitt 3).
