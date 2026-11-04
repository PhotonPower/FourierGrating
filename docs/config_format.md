# Konfigurationsformat (JSON)

Einheiten: Längen in µm, Winkel in Grad. Brechzahlen als Zahl oder `[Re, Im]`.
Die Implementierung erfolgt in M2 (`fg::io::loadConfig`); dieses Dokument ist die Spezifikation.

```json
{
  "dimension": "1D",                     // "1D" | "2D"
  "solver": "rcwa",                      // "rcwa" | "rk"
  "truncation": { "Nx": 15, "Nz": 0 },
  "period_um": { "x": 0.5, "z": 0.0 },
  "incidence": {
    "lambda_um": 0.6328,                 // Zahl oder {"from":..,"to":..,"steps":..}
    "phi_deg": 30, "psi_deg": 0,         // Polarwinkel zur y-Achse, Azimut gegen x
    "delta_deg": 90                      // 0 = s (TE), 90 = p (TM)
  },
  "materials": {
    "air":   { "n": 1.0 },
    "al":    { "n": [1.38, 7.62] },
    "co":    { "eps_tensor": [[[-8.19,16.38],[0,0],[-0.495,-0.106]],
                              [[0,0],[-8.19,16.38],[0,0]],
                              [[0.495,0.106],[0,0],[-8.19,16.38]]] }
  },
  "superstrate": "air",
  "substrate": "al",
  "layers": [                            // von unten (Substrat) nach oben
    { "type": "profile1d", "shape": "sinus",            // sinus | triangle | points
      "depth_um": 0.2, "below": "al", "above": "air",
      "slices": 0, "factorization": "auto" },     // auto | laurent | li_lamellar | FFF_iso | FFF_sym | FFF_rig
    { "type": "homogeneous", "thickness_um": 0.05, "material": "air" },
    { "type": "lamellar1d", "thickness_um": 0.3,
      "breakpoints_um": [0.0, 0.2, 0.5], "materials": ["al", "air"] },
    { "type": "polygon2d", "thickness_um": 0.2, "background": "ag",
      "inclusions": [ { "material": "air",
                        "vertices_um": [[0.325,0.325],[0.575,0.325],[0.575,0.575],[0.325,0.575]] } ] },
    { "type": "profile2d", "shape": "egg_crate", "depth_um": 0.1, "below": "si", "above": "air" }
  ],
  "rk":   { "tableau": "gauss2", "tol": 1e-8 },
  "rcwa": { "slices": 60 },
  "backend": "auto",                     // "auto" | "cpu" | "cuda"
  "nearfield": {
    "enabled": true,
    "horizontal_y_um": 0.1,              // horizontaler Schnitt
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
Kommentare (`//`) dienen nur der Erläuterung; echte Dateien sind reines JSON.

Ausgabe: `efficiencies.csv` (Spalten: lambda_um, phi, psi, delta, side, n, m, efficiency,
amp_s_re, amp_s_im, amp_p_re, amp_p_im), `diagnostics.json`, Nahfeld je Schnitt in
`nearfield_<schnitt>/` (`*.npy`, `meta.json`).
