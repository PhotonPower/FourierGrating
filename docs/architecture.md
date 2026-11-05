# Architektur

## Datenfluss einer Rechnung

```
JSON-Konfiguration (µm, Grad)
   └─ io::loadConfig ──► Stack (Materialien, Schichten), SolverOptions, Wavelength-Scan
        └─ Normierung (k0 = 2π/λ, alle Längen · k0)
             ├─ IndexMap (n,m) ↔ r, α_n, γ_m                         core/, physics/incidence
             ├─ je Schicht / Höhe y:
             │    Fourierkoeffizienten ε, 1/ε, ε_ij, N_iN_j           fourier/, geometry/
             │    Q_ε (Laurent | LiLamellar | FFF_iso | FFF_A | FFF_B) physics/qmatrix
             │    M(y) aus Q_ε, α, γ                                  physics/propagation_matrix
             ├─ Scheiben-S-Matrix
             │    RK:   ∏ P_h (implizites RK) → T_F → Ψ⁻¹T_FΨ → S       solvers/ (Löser A)
             │    RCWA: Eigenmoden W, Λ → stabile S = H G⁻¹             solvers/ (Löser B)
             ├─ Grenzflächen Substrat/Superstrat ↔ Referenzmedium      physics/modes
             ├─ Sternprodukt über den Stapel                           smatrix/
             ├─ Lösung: reflektierte/transmittierte Amplituden
             ├─ Effizienzen (Poynting), s/p-Zerlegung, Diagnosen      physics/efficiency
             └─ optional Nahfeld: interne Amplituden, Rekonstruktion,
                 Synthese auf Gitter, Schnitte → .npy + meta.json      field/
```

## Module

| Verzeichnis | Inhalt | Status |
|---|---|---|
| `core` | Typen, Fehler, `IndexMap` | fertig |
| `fourier` | Koeffizienten (Intervall, Rechteck, Polygon), Toeplitz 1D/2D | fertig (FFT-Pfad fehlt, M3/M7) |
| `geometry` | Material, Schichttypen, Stapel, Normalenfelder (1D-Profile, Algorithmus NF inkl. Rang-1-Feld, 2D-Profile) | fertig |
| `physics` | Einfall, `M`, isotrope Moden, Q-Matrizen, anisotrope Moden, Effizienzen | fertig bis auf FFF_A/B, anisotrope Moden (M5) |
| `smatrix` | T→S, Sternprodukt, interne Amplituden | fertig |
| `solvers` | Butcher, RK-Transfer, RCWA-Moden + stabile S, Schichtmodelle (`detail::LayerModel`, inkl. `Polygon2D`, `Profile2D`), Stapellöser mit Aufzeichnung, Feldkoeffizienten, Schnitte | 1D und 2D fertig (inkl. Nahfeld, M9) |
| `field` | Schnitte (`FieldCut`), NPY/JSON-Ausgabe | 1D und 2D fertig |
| `linalg`, `cuda` | Backend-Schnittstelle, CPU (Eigen), CUDA (cuBLAS/cuSOLVER Host-API), Dispatch | fertig; GPU-Abnahme auf dem Zielsystem offen |
| `io` | JSON-Konfiguration, CSV-Ausgabe | alle Schichttypen |
| `geometry/contour`, `fourier/contour` | Rationale Bézier-Konturen, Regionen (Kreis, Ellipse, NURBS, Löcher), Koeffizienten aus Konturintegralen (eq:contour) | M13 |
| `geometry/section`, `geometry/boundary` | `SectionProvider`/`CrossSection` für höhenabhängige Querschnitte (M14: CAD), Randstichproben für Algorithmus NF | M13 |
| `cad` | CAD-Slicer (`fg::cad`): STEP über Nurbforge (nur hier eingebunden), Prüfung, Bandzerlegung, Schichten (`docs/cad.md`) | M14 |
| `solvers/simulation` | `fg::Simulation`: Konfiguration und Schichtmodelle über viele Rechnungen, `setParameter` mit Teil-Neuaufbau | M12 |
| `python/bindings` | pybind11-Modul `_core` (nur Umwandlung, GIL freigegeben), Pakete `fouriergrating`, `fouriergrating_cuda12` (`docs/python.md`) | M12 |

## Größenordnungen

| Fall | n_f | Zustandsvektor 4n_f | S-Block n_h | Speicher je S-Matrix |
|---|---|---|---|---|
| 1D in-plane, N = 30 | 61 | (2·61) | 61 | 0,2 MB |
| 1D konisch, N = 30 | 61 | 244 | 122 | 1 MB |
| 2D, N = 10 | 441 | 1764 | 882 | 50 MB |
| 2D, N = 15 | 961 | 3844 | 1922 | 236 MB |
| 2D, N = 20 | 1681 | 6724 | 3362 | 723 MB |

Nahfeld benötigt je innere Grenzfläche zwei S-Matrizen (S_L, S_U), bei 2D und vielen Scheiben
also mehrere GB – auf der RTX 6000 Ada (48 GB) machbar, sonst Zwei-Durchlauf-Variante (Plan).

## Parallelität und Profil (M15)

* `fg::prof` (`fg/core/prof.hpp`): `FG_PROF_SCOPE("rcwa.eig")` misst Abschnitte, aggregiert über alle Threads
  (Summe, Anzahl, Maximum); abgeschaltet eine atomare Abfrage. Ausgabe: `fgsolve --profile`, Python
  `solve(..., profile=True)`.
* `fg::parallelFor` (`fg/core/parallel.hpp`): OpenMP, ohne OpenMP `std::thread` (ThreadSanitizer-Job);
  Ergebnisse nach Index, erste Ausnahme (kleinster Index) wird weitergereicht. `innerThreads()` verhindert
  Überbelegung (innere Schleifen nur, wenn der Aufrufer sequentiell ist).
* `StackSolver::assemble` erzeugt je Schicht/Scheibe eine Aufgabe, rechnet sie blockweise parallel und verkettet
  in Stapelreihenfolge über `BasicStarChain` (sequentiell oder Baum fester Form). Lazy aufgebaute Schichtdaten
  (`SectionedModel`: Schnitt-Cache, Tschebyschow-Interpolant) sind mit Mutex bzw. `call_once` geschützt.
* GPU-residente Kette (`fg/solvers/resident.hpp`): `linalg::DeviceMatrix` (RAII, Speicherpool je Backend),
  Moden, Scheiben-S und Sternprodukte auf dem Gerät; ein Stream je Thread, Synchronisation beim Übergeben.

## Backend-Auswahl

`fg::linalg::makeBackend("auto", n)` wählt CUDA ab einer Matrixgröße (Standard 1000), sonst CPU.
Der CPU-Pfad ist Referenz; GPU-Ergebnisse müssen auf 1e-10 relativ übereinstimmen.
