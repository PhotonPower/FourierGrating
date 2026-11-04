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
| `linalg` | Backend-Schnittstelle, CPU (Eigen), CUDA | CPU fertig, CUDA M10 |
| `io` | JSON-Konfiguration, CSV-Ausgabe | alle Schichttypen |

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

## Backend-Auswahl

`fg::linalg::makeBackend("auto", n)` wählt CUDA ab einer Matrixgröße (Standard 1000), sonst CPU.
Der CPU-Pfad ist Referenz; GPU-Ergebnisse müssen auf 1e-10 relativ übereinstimmen.
