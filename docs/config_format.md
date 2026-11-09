# Konfigurationsformat (JSON)

Einheiten: Längen in µm, Winkel in Grad. Brechzahlen als Zahl oder `[Re, Im]`.
Eingelesen von `fg::io::loadConfig` bzw. `fg::io::loadSweep`; dieses Dokument ist die Spezifikation.
Eine Einführung mit Beispielen steht im Anwenderhandbuch `docs/user_guide.md`.

```json
{
  "dimension": "1D",                     // "1D" | "2D"
  "solver": "rcwa",                      // "rcwa" | "rk"
  "truncation": { "Nx": 15, "Nz": 0 },
  "period_um": { "x": 0.5, "z": 0.0 },
  "incidence": {
    "lambda_um": 0.6328,                 // Zahl, Liste [..] oder {"from":..,"to":..,"steps":..}
    "phi_deg": 30, "psi_deg": 0,         // Polarwinkel zur y-Achse, Azimut gegen x (ebenso Zahl/Liste/Bereich)
    "delta_deg": 90                      // 0 = s (TE), 90 = p (TM); Scan = kartesisches Produkt λ × φ × ψ × δ
  },
  "sweep": {                             // optional (M11): Parameter-Sweep über einen beliebigen Zahlenwert
    "parameter": "/layers/0/thickness_um",   // JSON-Pointer (RFC 6901) in diese Datei
    "values": { "from": 0.1, "to": 0.5, "steps": 9 }  // Zahl, Liste oder Bereich
  },
  "materials": {
    "air":   { "n": 1.0 },
    "al":    { "n": [1.38, 7.62] },
    "co":    { "eps_tensor": [[[-8.19,16.38],[0,0],[-0.495,-0.106]],
                              [[0,0],[-8.19,16.38],[0,0]],
                              [[0.495,0.106],[0,0],[-8.19,16.38]]] },
    "lnb":   { "eps_principal": [4.98, 4.98, 4.6], "euler_deg": [0, 90, 0] }  // R diag R^T, z-x-z
  },
  "superstrate": "air",
  "substrate": "al",
  "layers": [                            // von unten (Substrat) nach oben
    { "type": "profile1d", "shape": "sinus",            // sinus | triangle | points
      "depth_um": 0.2, "below": "al", "above": "air",
      "slices": 0, "factorization": "auto" },     // isotrop: auto | laurent | li_lamellar | FFF_iso | FFF_sym | FFF_rig
                                                  // anisotrop: auto (= FFF_A) | FFF_A | FFF_B | laurent
    { "type": "homogeneous", "thickness_um": 0.05, "material": "air" },
    { "type": "lamellar1d", "thickness_um": 0.3,
      "breakpoints_um": [0.0, 0.2, 0.5], "materials": ["al", "air"] },
    { "type": "polygon2d", "thickness_um": 0.2, "background": "ag",   // nur "dimension": "2D"
      "inclusions": [ { "material": "air",                           // Polygon (x, z), beliebige Orientierung,
                        "vertices_um": [[0.325,0.325],[0.575,0.325],[0.575,0.575],[0.325,0.575]] },
                      { "material": "air", "rectangle_um": { "center": [0.0, 0.0], "size": [0.2, 0.1] } } ],
      "normal_field": { "method": "tensor", "grid": 512, "sigma": 3 },  // Algorithmus NF; method: tensor | vector
      "factorization": "auto" },     // auto: FFF_sym (isotrop) / FFF_B (anisotrop) | FFF_sym | FFF_iso | FFF_B | laurent
    { "type": "contour2d", "thickness_um": 0.2, "background": "air", // gekrümmte Ränder, exakte Koeffizienten (M13)
      "regions": [
        { "material": "si", "circle": { "center_um": [0.25, 0.25], "radius_um": 0.1 } },
        { "material": "si", "ellipse": { "center_um": [0.6, 0.6], "radii_um": [0.12, 0.06], "angle_deg": 30 } },
        { "material": "sio2", "polygon_um": [[0, 0], [0.1, 0], [0.1, 0.1]] },
        { "material": "si", "nurbs": { "degree": 2, "knots": [...], "points_um": [[x, z], ...], "weights": [...] },
          "holes": [ { "circle": { "center_um": [0.6, 0.2], "radius_um": 0.02 } } ] } ],
      "normal_field": { "grid": 512, "sigma": 3 },
      "factorization": "auto" },     // wie polygon2d: auto: FFF_sym / FFF_B | FFF_iso | laurent (kein FFF_rig)
    { "type": "profile2d", "shape": "egg_crate", "depth_um": 0.1,     // egg_crate | sinus_x
      "below": "si", "above": "air", "x0_um": 0, "z0_um": 0,
      "slices": 0, "z_samples": 0, "factorization": "auto" }        // auto: FFF_rig (isotrop) / FFF_A (anisotrop)
  ],
  "rk":   { "tableau": "gauss2", "tol": 1e-8 },
  "rcwa": { "slices": 60 },
  "backend": "auto",                     // "auto" | "cpu" | "cuda" (auto: CUDA ab backend_threshold, falls vorhanden)
  "backend_threshold": 1000,             // Matrixgröße, ab der auto die GPU nutzt
  "split_polarizations": false,          // 1D, psi = 0, isotrop: TE/TM getrennt mit 2 n_f (schneller); false = 4 n_f
  "fft_sampling": false,                 // abgetastete Koeffizienten per FFT jeder Länge; false = v0.10 (DFT)
  "threads": { "scan": "auto", "slices": "auto", "blas": 1 },  // Threadmodell (M15), siehe unten
  "reduction": "sequential",             // Sternprodukt-Kette: "sequential" (v0.10) | "tree" (parallel, M15)
  "precision": "auto",                   // auto | fp64 | mixed | fp32 (M16), siehe unten
  "precision_fallback": true,            // fp32: verlustfrei und |Σe − 1| > 1e-3 -> in fp64 wiederholen
  "growth_limit": 0,                     // Scheibenkriterium; 0 = auto (fp64 1e8, fp32/mixed 4e3)
  "checkpoint_every": 1,                 // Gradienten (M19): nur jedes k-te Präfixprodukt speichern, Rest neu rechnen
  "nearfield": {
    "enabled": true,
    "horizontal_y_um": 0.1,              // horizontaler Schnitt (1D: Linie in x; 2D: Karte x-z, nz Punkte je Periode)
    "vertical_z_um": 0.0,                // vertikaler Schnitt (x-y-Ebene)
    "y_range_um": [-0.2, 0.5],
    "nx": 400, "ny": 300, "nz": 200, "periods": 2,
    "quantities": ["E", "H", "absE2", "poynting"],
    "reconstruction": "fff"              // "direct" | "fff"
  },
  "output": { "dir": "out/example", "efficiencies": "efficiencies.csv" }
}
```

Formen für `profile1d.shape`: `sinus` (g = a/2 (1 + cos Kx)), `triangle` (Parameter
`blaze_deg`, `apex_deg`), `points` (Stützstellen `x_um`, `y_um`, periodischer kubischer Spline).
Formen für `profile2d.shape`: `egg_crate` (g = a/4 (1 + cos Kx(x−x0)) (1 + cos Kz(z−z0))), `sinus_x`
(g = a/2 (1 + cos Kx(x−x0)), z-invariant). Polygone dürfen über den Zellrand hinausragen (periodische Bilder
dürfen sich nicht überlappen); Einschlüsse dürfen sich nicht überlappen. Bei `"dimension": "2D"` ist
`period_um.z > 0` Pflicht; `lamellar1d`- und isotrope `profile1d`-Schichten sind in 2D-Stapeln erlaubt (z-invariant).
`cad` (M14, Build mit `FG_ENABLE_NURBFORGE`): STEP-Körper, beim Laden in Bänder expandiert (homogen,
`contour2d`, y-variabel); Schlüssel und Prüfcodes in [`cad.md`](cad.md).

`contour2d` (M13): Regionen mit äußerer Kontur und Löchern, paarweise disjunkt (auch mit den periodischen
Bildern; geprüft, sonst Fehler). Formen: `circle`, `ellipse` (exakt als rationale Quadratiken), `polygon_um`,
`nurbs` (geklemmter Knotenvektor, geschlossen: erster = letzter Punkt; Zerlegung in Bézierstücke). Die
Orientierung ist beliebig. Koeffizienten aus Konturintegralen (Plan eq:contour, ≈ 1e-14), Normalenfeld aus
Randstichproben der Kurven (Algorithmus NF, M13). Ein `contour2d`-Polygon ist identisch mit `polygon2d`.
Kommentare (`//`) dienen nur der Erläuterung; echte Dateien sind reines JSON.

**Threads (M15).** Drei Ebenen: Scanpunkte (`threads.scan`), Schichten/Scheiben einer Rechnung
(`threads.slices`: RCWA-Scheiben, RK-Stabilitätsscheiben, homogene Schichten; dazu Normalenfeld-Raster,
Konturkoeffizienten und Tschebyschow-Knoten, wenn die Ebene darüber sequentiell ist) und Eigen innerhalb einer
Matrixoperation (`threads.blas`, Standard 1). Das Produkt bleibt ≤ der Threadzahl `FG_NUM_THREADS` (sonst
`OMP_NUM_THREADS`, sonst alle Kerne). `auto` gibt der Scanebene min(Scanpunkte, Kerne) Threads und den
Scheiben den Rest; eine einzelne Wellenlänge nutzt also alle Kerne für ihre Scheiben. Die Ergebnisse hängen
nicht von der Threadzahl ab (bitgleich). `"reduction": "tree"` verkettet die Scheiben als balancierten
Binärbaum fester Form (parallel; bitgleich für jede Threadzahl, ≤ 1e-12 von `sequential`); Standard ist die
sequentielle Kette (bitgleich zu v0.10). Mit der GPU (`backend` cuda/auto über der Schwelle) bleiben
Scheiben-S-Matrizen und Kette auf dem Gerät; höchstens `FG_GPU_STREAMS` (Standard 4) Scheiben gleichzeitig.
Benutzerfunktionen (`profile2d` mit Funktion in C++, eigene `SectionProvider`) werden parallel aufgerufen und
müssen threadsicher sein. Sweep-Werte laufen nacheinander.

**Genauigkeit (M16).** `precision` wählt die Genauigkeit der O(n³)-Operationen (GEMM, LU, Eigenwerte); alles andere
(Koeffizienten, Q, M, Exponentiale, stabile S-Formeln, Effizienzen) bleibt FP64. `fp64`: wie bisher (bitgleich).
`mixed`: GEMM und Eigenwerte in einfacher Genauigkeit, LU in einfacher Genauigkeit mit FP64-Nachiteration.
`fp32`: alles in einfacher Genauigkeit; `rk.tol` mindestens 1e-5. `auto` (Standard): auf der CPU immer `fp64`; auf
der GPU `fp32` auf einer GeForce RTX 3090, sonst `fp64` (Tabelle `src/cuda/device_policy.cpp`). `FG_PRECISION`
überschreibt die Konfiguration. Bei `fp32` (Standard) bzw. `precision_fallback: true` wird eine verlustfreie
Struktur mit |Σe − 1| > 1e-3 in fp64 wiederholt (Warnung). Das Ergebnis meldet `diagnostics.precision`,
`precision_fallback`, `refinement_steps`. TF32 wird nie benutzt; zusätzlich kann `NVIDIA_TF32_OVERRIDE=0` gesetzt
werden. Genauigkeit und Laufzeiten: `docs/validation.md`, `docs/benchmarks.md`.

**Koordinatentransformation (M17).** Adaptive räumliche Auflösung nach Granet (ASR) für Metalle und scharfe Kanten:

```jsonc
"coordinate_transform": {
  "x": {"type": "asr", "eta": 0.8, "breakpoints_um": "auto"},   // oder eine Liste [x0, x1, ...]
  "z": {"type": "asr", "eta": 0.5, "breakpoints_um": [0.1, 0.4]} // nur 2D
}
```

`eta` (0 ≤ η < 1) verdichtet das Fourier-Raster an den Bruchstellen um den Faktor 1 − η; `eta = 0` ist keine
Transformation (bitgleich). `"auto"`: alle Materialsprünge der Lamellenschichten und die Eckkoordinaten der Polygone
(`polygon2d`, polygonale `contour2d`-Gebiete) entlang der Achse; ohne solche (z. B. nur Profile) muss eine Liste
angegeben werden. Mit einer einzigen Bruchstelle ist f analytisch (schnellste Konvergenz der Halbräume); bei mehreren
ungleich langen Intervallen springt f‴ (algebraische Konvergenz ~N⁻⁵). Gilt für alle Schichttypen; `split_polarizations`
wird dabei nicht genutzt. Empfehlung für Metall-Lamellen in TM: η = 0,8 … 0,95 (V19 in `docs/validation.md`).
Nahfeld: die Ausgabe bleibt auf dem physikalischen Raster in x, z; die Feldkoeffizienten (`fieldCoefficients`)
beziehen sich dann auf die transformierten Koordinaten X, Z.
**Isolierte Strukturen, PML (M18).** Ein einzelnes Objekt (Steg, Graben, Kontaktloch, Wellenleiterende) wird in
einer Rechenzelle mit absorbierenden Randzonen (PML, komplexe Koordinatenstreckung) gerechnet:

```jsonc
"aperiodic": {
  "pml": {"x": {"width_um": 0.3, "beta": 1.0}, "z": {"width_um": 0.3}},  // z nur 2D
  "formulation": "auto",          // oder "scattered" (Streufeldformulierung auch bei Normaleinfall)
  "excitation": "plane_wave"      // oder {"mode": 0} bzw. {"n_eff": 1.9}: geführte Mode des Superstrats
},
"superstrate_structure": {"type": "lamellar1d", "breakpoints_um": [0, 1.3, 1.7, 3.0], "materials": ["air", "core", "air"]},
"substrate_structure": {...},     // y-invariante Querschnitte (lamellar1d, polygon2d, contour2d), Dicke entfällt
"scattered_field": false          // Streufeldformulierung für periodische Gitter (Test, V20d)
```

Physikalischer Bereich: [w, d − w] (bzw. das Rechteck in 2D); Strukturen dürfen nicht in die Zonen reichen
(geprüft, auch für CAD-Schnitte; Filme dürfen durchlaufen). Normaleinfall: Totalfeld minus Hintergrundstapel (jede
strukturierte Schicht durch ihr Material am Zellrand ersetzt). Schräger Einfall: automatisch Streufeldformulierung
(Quelle σ = (M − M_bg) F_bg in den strukturierten Schichten, Streufeld ohne Blochphase). Ausgabe statt der
Effizienzen (`efficiencies.csv`, Seite `S_up`, `S_down`, `S_abs`; Python `result.scattering`): Leistung des
Streufeldes durch die Ebenen y = Oberkante und y = 0 über den physikalischen Bereich, bezogen auf die Intensität der
einfallenden Welle (1D: Streubreite in µm, 2D: Streuquerschnitt in µm²), `lateral`: durch die Seiten der Box
(physikalischer Bereich × Schichten), `absorbed`: Fluss des Totalfeldes in die Box minus Hintergrund. Nur die Summe
up + down + lateral ist von der Zellbreite unabhängig; die Aufteilung auf oben/unten konvergiert bei schrägem Einfall
langsam (~1/L, streifende Anteile). Modenanregung: `efficiencies.csv` mit `MR`/`MT` (n = Modennummer,
efficiency = Leistung) und `Mrad`; Python `result.modal`. Nahfeld: `"nearfield": {"field": "scattered"}` gibt das
Streufeld statt des Totalfeldes aus. PML und `coordinate_transform` (ASR) schließen sich aus.

**Rotationssymmetrische Strukturen (M20, in Arbeit).** Strukturen um die Achse y (x = z = 0) auf einem Schichtstapel,
angeregt durch einen Strahl auf der Achse; `period_um`, `truncation` und `dimension` entfallen:

```jsonc
"cylindrical": {
  "symmetry": "axial",
  "harmonics": "auto",                                    // oder Liste, z. B. [-1, 1]
  "spectrum": {"kappa_max": 6, "dk_factor": 1.0, "nodes_per_panel": 8,
               "deform": "auto"},                          // "auto" | false | {"delta": 0.05} (relativ zu k_b - k_a)
  "factorization": "li_r"                                  // oder "laurent"
},
"incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.8,
              "polarization": "x"},  // x | z | circular_plus | circular_minus | radial | azimuthal
"layers": [{"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["sin", "air"]}]
```

`rings`: n Radien und n + 1 Materialien von innen nach außen (das letzte gilt außerhalb). Erlaubt sind `homogeneous`
und `rings`, isotrope Materialien, homogene Halbräume. Ausgabe in `efficiencies.csv`: Seiten `C_R`, `C_T`, `C_A`
(Leistungen bezogen auf den einfallenden Strahl; `C_A` = 1 − R − T enthält auch den seitlichen Abfluss geführter
Moden) und je Harmonischer `C_Rm`, `C_Tm` (n = m). Ebene Wellen, Nahfeld, `rings_profile` und CAD folgen.

Bei einem Sweep beginnen die CSV-Zeilen mit `sweep_parameter, sweep_value`; das Nahfeld wird für den ersten
Scanpunkt des ersten Sweep-Werts berechnet. Darstellung: `scripts/plot_efficiencies.py`.

Ausgabe: `efficiencies.csv` (Spalten: lambda_um, phi, psi, delta, side, n, m, efficiency,
eff_s, eff_p, amp_s_re, amp_s_im, amp_p_re, amp_p_im; s/p nur in isotropen Medien), `diagnostics.json`, Nahfeld je Schnitt in
`nearfield_<schnitt>/` (`*.npy`, `meta.json`).
