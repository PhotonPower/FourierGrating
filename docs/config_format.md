# Konfigurationsformat (JSON)

Einheiten: Längen in µm, Winkel in Grad. Brechzahlen als Zahl, `[Re, Im]` oder dispersives Modell (M22, unten).
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
    "lnb":   { "eps_principal": [4.98, 4.98, 4.6], "euler_deg": [0, 90, 0] },  // R diag R^T, z-x-z
    "sio2":  { "n": { "sellmeier": { "A": 1.0, "B": [0.6961663, 0.4079426, 0.8974794],      // M22: Modelle,
                                     "C_um2": [0.00467914826, 0.0135120631, 97.9340025] } } }  // siehe unten
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
  "backend_threshold": 160,              // Matrixgröße, ab der auto die GPU nutzt
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
    "quantities": ["E", "H", "absE2"],     // derzeit ohne Wirkung: es werden stets alle Größen geschrieben
    "reconstruction": "fff",             // "direct" | "fff"
    "extent_um": 0                       // nur "cylindrical" (M20): x, z in [-extent, extent]; 0 = 1,5 R_s
  },
  "output": { "dir": "out/example", "efficiencies": "efficiencies.csv",
              "materials": "materials.csv" }    // optional (M22): n, ε je Scan-Wellenlänge
}
```

Formen für `profile1d.shape`: `sinus` (g = a/2 (1 + cos K(x − x0)), optional `x0_um`), `triangle` (Parameter
`apex_x_um`: Flanke von (0, 0) zur Spitze (`apex_x_um`, a) und zurück nach (d, 0)), `points` (Stützstellen `x_um`,
`y_um`, x aufsteigend in [0, d), periodisches Polygon, stückweise linear).
`lamellar1d`: `breakpoints_um` beginnt bei 0 und endet bei der Periode d; das wird derzeit **nicht geprüft** (eine
falsche letzte Stelle, z. B. nach einer Änderung von `period_um.x`, liefert ohne Meldung ein falsches Ergebnis).
Formen für `profile2d.shape`: `egg_crate` (g = a/4 (1 + cos Kx(x−x0)) (1 + cos Kz(z−z0))), `sinus_x`
(g = a/2 (1 + cos Kx(x−x0)), z-invariant). Polygone dürfen über den Zellrand hinausragen (periodische Bilder
dürfen sich nicht überlappen); Einschlüsse dürfen sich nicht überlappen. Bei `"dimension": "2D"` ist
`period_um.z > 0` Pflicht; `lamellar1d`- und isotrope `profile1d`-Schichten sind in 2D-Stapeln erlaubt (z-invariant).
`cad` (M14, Build mit `FG_ENABLE_NURBFORGE`): STEP-Körper, beim Laden in Bänder expandiert (homogen,
`contour2d`, y-variabel); in 1D-Stapeln (M23) z-invariante Körper, ein Schnitt bei `z_cut_um`, Bänder homogen,
`lamellar1d` oder y-variabel. Schlüssel und Prüfcodes in [`cad.md`](cad.md).

`contour2d` (M13): Regionen mit äußerer Kontur und Löchern, paarweise disjunkt (auch mit den periodischen
Bildern; geprüft, sonst Fehler). Formen: `circle`, `ellipse` (exakt als rationale Quadratiken), `polygon_um`,
`nurbs` (geklemmter Knotenvektor, geschlossen: erster = letzter Punkt; Zerlegung in Bézierstücke). Die
Orientierung ist beliebig. Koeffizienten aus Konturintegralen (Plan eq:contour, ≈ 1e-14), Normalenfeld aus
Randstichproben der Kurven (Algorithmus NF, M13). Ein `contour2d`-Polygon ist identisch mit `polygon2d`.
Kommentare (`//`) dienen nur der Erläuterung; echte Dateien sind reines JSON.

**Dispersive Materialien (M22, Plan Abschnitt „Dispersive Materialien“).** Ein Materialwert unter `"n"`, `"eps"`
oder in einem Eintrag von `"eps_principal"` darf statt einer Zahl ein Modell sein; es wird an jedem Scanpunkt
bei dessen Wellenlänge ausgewertet. Die Geometrie hängt nicht von λ ab.

```json
"sio2": { "n": { "sellmeier": { "A": 1.0, "B": [0.6961663, 0.4079426, 0.8974794],
                                "C_um2": [0.00467914826, 0.0135120631, 97.9340025],
                                "range_um": [0.21, 6.7] } } },          // range_um optional
"tio2": { "n": { "cauchy": { "A": 2.40, "B_um2": 0.03, "C_um4": 0.0, "k": 0.0 } } },
"au":   { "n": { "table": { "file": "materials/au.csv", "unit": "um",  // oder "data": [[λ, n, k], ...]
                            "interpolation": "cubic", "extrapolate": "error" } } },
"ag":   { "eps": { "table": { "data": [[0.40, -4.4, 0.21], [0.50, -9.8, 0.31]] } } },
"al":   { "eps": { "drude_lorentz": { "eps_inf": 1.0, "omega_p_eV": 14.98, "f0": 0.523, "gamma0_eV": 0.047,
            "oscillators": [ { "f": 0.227, "omega_eV": 0.162, "gamma_eV": 0.333 }, ... ] } } },  // gekürzt: Rakić hat 4
"lnb":  { "eps_principal": [ { "sellmeier": {...} }, { "sellmeier": {...} }, 4.6 ], "euler_deg": [0, 90, 0] }
```

| Modell | liefert | Formel / Werte |
| --- | --- | --- |
| `sellmeier` | ε | ε = A + Σ B_k λ²/(λ² − C_k), C_k in µm², A Standard 1; reell. `range_um` (optional): außerhalb Fehler |
| `cauchy` | n | n = A + B/λ² + C/λ⁴ + i k (`B_um2`, `C_um4`, `k` Standard 0) |
| `drude_lorentz` | ε | ε = ε∞ − f₀ωp²/(ω(ω + iΓ₀)) + Σ f_j ωp²/(ω_j² − ω² − iωΓ_j), Energien in eV, ħω = 1,239841984/λ |
| `table` | Schlüssel | unter `"n"`: Spalten λ, n[, k]; unter `"eps"`/`"eps_principal"`: λ, Re ε[, Im ε] |

* Liefert ein Modell n, ist ε = n²; liefert es ε, gilt das auch unter `"n"`. Verluste: Im ε > 0, k ≥ 0.
* Tabellen: `"file"` (CSV; Trennzeichen Komma, Semikolon oder Leerraum; `#`-Zeilen und eine Kopfzeile ohne
  Zahlen werden übersprungen; relativ zum Verzeichnis der Konfiguration) **oder** `"data"`. `"unit"` der ersten
  Spalte: `um` (Standard), `nm`, `eV`. `"interpolation"`: `linear` (Standard) oder `cubic` (Not-a-knot-Spline, ab
  4 Zeilen); jede Spalte für sich in λ. Außerhalb der Tabelle: Fehler, mit `"extrapolate": "constant"` Randwert
  und Warnung in der Diagnose.
* Dispersive Einträge in `eps_tensor` sind nicht vorgesehen (`eps_principal` mit `euler_deg` verwenden).
* Frühe Prüfung beim Laden (auch jeder Sweep-Wert): Bereich, Sellmeier-Pol, verlustfreies Superstrat an jedem
  Scanpunkt. Die Meldung nennt den JSON-Pointer, z. B. `/materials/au/n/table: lambda 0.25 um outside [0.3, 1.9] um`.
* Modellwerte sind Zahlen der Konfiguration: Sweeps, `set_parameter` und Gradienten sprechen sie per Pointer an,
  z. B. `/materials/sio2/n/sellmeier/B/0`, `/materials/ag/eps/table/data/3/2` (Tabellenwerte nur mit `"data"`;
  die Wellenlängen-Spalte ist kein Gradientenparameter).
* Ausgabe: `Result::diag.materials` (n, ε der dispersiven Materialien je Scanpunkt), Terminal `material …`,
  optional `materials.csv` mit den Spalten `lambda_um, material, component, n_re, n_im, eps_re, eps_im`
  (`component`: `iso` oder Hauptachse `1`–`3`; alle Materialien außer `eps_tensor`).

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
"scattered_field": false          // Streufeldformulierung für periodische Gitter (Test, V20d) und Strahlen (M20, V21e)
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
angeregt durch einen Strahl auf der Achse oder eine ebene Welle; `period_um`, `truncation` und `dimension` entfallen:

```jsonc
"cylindrical": {
  "symmetry": "axial",
  "harmonics": "auto",                                    // oder Liste, z. B. [-1, 1]
  "spectrum": {"kappa_max": 6, "dk_factor": 1.0, "nodes_per_panel": 8,
               "deform": "auto"},                          // "auto" | false | {"delta": 0.05} (relativ zu k_b - k_a)
  "factorization": "li_r",                                 // oder "laurent"
  "optical_theorem": false                                 // ebene Welle: optisches Theorem als Kontrolle (teurer)
},
"incidence": {"type": "gaussian_beam", "lambda_um": 1.0, "waist_um": 0.5, "focus_y_um": 0.8,
              "polarization": "x"},  // x | z | circular_plus | circular_minus | radial | azimuthal
"layers": [{"type": "rings", "thickness_um": 0.3, "radii_um": [0.4], "materials": ["sin", "air"]}]
```

`rings`: n Radien und n + 1 Materialien von innen nach außen (das letzte gilt außerhalb). Erlaubt sind `homogeneous`
und `rings`, isotrope Materialien, homogene Halbräume. Ausgabe in `efficiencies.csv`: Seiten `C_R`, `C_T`, `C_A`
(Leistungen bezogen auf den einfallenden Strahl; `C_A` = 1 − R − T enthält auch den seitlichen Abfluss geführter
Moden) und je Harmonischer `C_Rm`, `C_Tm` (n = m). Sind Ringmaterialien verlustbehaftet, steht in Python
`result.cylindrical["absorbed_rings"]` die in den Ringen absorbierte Leistung (Volumenintegral); `C_A` minus
dieser Wert ist der seitliche Abfluss plus die Absorption außerhalb der Ringe. Mit `"scattered_field": true`
rechnet auch ein Strahl in der Streufeldformulierung (Hintergrund = Stapel mit dem äußeren Material jeder
Ringschicht, Quelle σ = (M − M_bg) F_bg); die Ergebnisse sind dieselben wie beim Totalfeld (Test V21e).

Ebene Welle (Streufeldformulierung, Hintergrund = Stapel mit dem äußersten Ringmaterial): Winkel wie beim
kartesischen Einfall, Polarisation über `delta_deg`; Strahlparameter entfallen.

```jsonc
"incidence": {"type": "plane_wave", "lambda_um": 1.0, "phi_deg": 30, "psi_deg": 0, "delta_deg": 0}
```

Harmonische `"auto"`: bei Normaleinfall m = ±1, sonst |m| ≤ m_cut = ⌈x + 4x^{1/3} + 2⌉ mit x = k₀ n_sup sin φ R_s
(R_s: größter Ringradius; darüber verschwindet die Quelle). Ausgabe: Streuquerschnitte in µm² `C_Sup` (ins
Superstrat), `C_Sdown` (ins Substrat), `C_Ssca` (Summe), je Harmonischer `C_Supm`, `C_Sdownm` (n = m); Python
`result.cylindrical["scattering_um2"]` usw. Dazu `C_Sabs` (Absorptionsquerschnitt: ½∫Im ε |E|² dV über die Ringe,
bezogen auf die einfallende Intensität), `C_Sext` = `C_Ssca` + `C_Sabs` und als Kontrolle `C_Sext_ot` (optisches
Theorem in Volumenform, ½ Im ∫(ε − ε_∞) E·Ē_bg dV). Die Volumenintegrale laufen nur bei verlustbehafteten Ringen
(sonst σ_abs = 0, σ_ext = σ_sca) oder mit `"optical_theorem": true` (etwa doppelte Rechenzeit); `C_Sext_ot` fehlt
sonst, ebenso bei verlustbehaftetem Hintergrund (Python dann `nan`). Python
`absorption_um2`, `extinction_um2`, `extinction_theorem_um2`. `rings_profile` und CAD folgen.

Nahfeld (`"nearfield"`, M20): horizontaler Schnitt bei `horizontal_y_um` als Karte x–z, vertikaler Schnitt x–y in der
Ebene z = `vertical_z_um` (z = 0: Meridianebene); x und z laufen symmetrisch über [−`extent_um`, `extent_um`] mit
`nx` bzw. `nz` Punkten einschließlich der Ränder (Standard 1,5 R_s; `periods` entfällt). Strahlen liefern das
Totalfeld, ebene Wellen Total- oder Streufeld (`"field"`). `"reconstruction": "fff"` setzt bei Li in r das Feld
aus den stetigen Hilfsgrößen G± zusammen (D_r stetig an den Wänden). Auf einem deformierten Weg sind Punkte mit
δ r > ln 10³ `nan` (Warnung).

Bei einem Sweep beginnen die CSV-Zeilen mit `sweep_parameter, sweep_value`; das Nahfeld wird für den ersten
Scanpunkt des ersten Sweep-Werts berechnet. Darstellung: `scripts/plot_efficiencies.py`.

Ausgabe: `efficiencies.csv` (Spalten: lambda_um, phi_deg, psi_deg, delta_deg, side, n, m, efficiency,
eff_s, eff_p, amp_s_re, amp_s_im, amp_p_re, amp_p_im; s/p nur in isotropen Medien), Nahfeld je Schnitt in
`nearfield_<schnitt>/` (`Ex…Hz.npy`, `absE2.npy`, `absH2.npy`, `Sy.npy`, `material.npy`, `meta.json`), mit
`fgsolve --profile` zusätzlich `profile.json`, mit `output.materials` die Materialwerte (M22, siehe oben). Eine Datei `diagnostics.json` wird nicht geschrieben; die Diagnose steht
auf dem Terminal bzw. in Python in `Result.diagnostics`. Relative `output.dir` gelten bei `fgsolve` gegen das
Arbeitsverzeichnis, in Python (`load_config`) gegen die Konfigurationsdatei.
