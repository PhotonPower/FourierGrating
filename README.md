# FourierGrating

Löser für die Beugung elektromagnetischer Wellen an periodischen Strukturen nach der
**Fourier-Modal-Methode (FMM)** in zwei Ausprägungen:

* **Differentialmethode** – numerische Integration der Propagationsgleichung `dF/dy = M(y) F`
  mit einem **impliziten Runge–Kutta-Verfahren** (Gauß–Legendre, optional Radau IIA),
  stabilisiert durch die S-Matrix,
* **RCWA** – Scheibenzerlegung mit Eigenmoden je Scheibe.

Beide nutzen die schnelle Fourierfaktorisierung (FFF, Normalenvektor-Methode) nach
Nevière/Popov, *Light Propagation in Periodic Media* (2003), mit den in `docs/book_errata.md`
dokumentierten Korrekturen.

| | isotrop | anisotrop (voller ε-Tensor) |
|---|---|---|
| **1D** (eine Fouriertransformation) | M2–M4 | M5 |
| **2D** (zwei Fouriertransformationen) | M7 | M8 |
| **Nahfeld** (horizontaler + vertikaler Schnitt) | M6 (1D), M9 (2D) | M6, M9 |
| **CUDA** (RTX 6000 Ada, sm_89) | M10 | M10 |

> **Status (v0.2):** M0–M2 erledigt. 1D-Gitter aus homogenen Schichten und Lamellen
> (isotrop, in-plane und konisch) lassen sich mit **beiden Lösern** (RK und RCWA) rechnen;
> Validierung V1 (Fresnel), V2 (Airy), V3 (Energie, RK = RCWA, Reziprozität) ist grün.
> Als Nächstes: M3 (Profile y = g(x) mit FFF). Nicht implementierte Teile werfen
> `fg::NotImplemented` mit Verweis auf den Meilenstein. Fortschritt: `docs/roadmap.md`.

## Dokumentation

| Dokument | Inhalt |
|---|---|
| [`docs/implementation_plan/main.pdf`](docs/implementation_plan/main.pdf) | **Verbindliche Spezifikation**: alle Formeln mit Herleitungen, Algorithmen, Tests, Meilensteine (LaTeX-Quelle daneben) |
| [`AGENTS.md`](AGENTS.md) | Arbeitsregeln für KI-Agenten und Entwickler |
| [`docs/architecture.md`](docs/architecture.md) | Module, Datenfluss, Schnittstellen |
| [`docs/roadmap.md`](docs/roadmap.md) | Meilensteine M0–M11 mit Abnahmekriterien und Status |
| [`docs/validation.md`](docs/validation.md) | Referenzfälle mit Sollwerten |
| [`docs/conventions.md`](docs/conventions.md) | Zeitkonvention, normierte Einheiten, Koordinaten, Indizes |
| [`docs/config_format.md`](docs/config_format.md) | JSON-Eingabeformat |
| [`docs/build_almalinux10.md`](docs/build_almalinux10.md) | Installation auf AlmaLinux 10 mit CUDA |
| [`docs/book_errata.md`](docs/book_errata.md) | Fehler im Buch, die die Implementierung betreffen |

## Schnellstart

```bash
# AlmaLinux 10: Abhängigkeiten (einmalig, siehe docs/build_almalinux10.md)
sudo ./scripts/setup_almalinux10.sh            # --with-cuda für das CUDA-Toolkit

cmake --preset cpu-release                      # oder cuda-release
cmake --build --preset cpu-release
ctest --preset cpu-release

./build/cpu-release/apps/fgsolve --info
./build/cpu-release/apps/fgsolve examples/v3_lamellar_dielectric.json   # schreibt out/v3/efficiencies.csv
```

Plan als PDF neu erzeugen: `./scripts/build_docs.sh`.

## Verzeichnisse

```
include/fg/   öffentliche Header (core, fourier, geometry, physics, smatrix, solvers, field, linalg, io)
src/          Implementierungen, src/cuda/ CUDA-Kernel und -Backend
apps/         Kommandozeilenprogramm fgsolve
tests/        GoogleTest-Unit-Tests, Validierungsfälle
examples/     JSON-Konfigurationen der Referenzfälle
scripts/      Setup, Doku-Build, Plot-Skripte (Python)
docs/         Spezifikation und Begleitdokumente
```

## Lizenz

MIT, siehe [`LICENSE`](LICENSE).
