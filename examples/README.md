# Beispiele / Referenzfälle

JSON-Konfigurationen nach `docs/config_format.md`. Das optionale Feld `reference` enthält die
Sollwerte aus `docs/validation.md`; Validierungstests (`tests/validation/`) lesen es aus.

| Datei | Fall |
|---|---|
| `v1_fresnel_interface.json` | V1, ebene Grenzfläche, konischer Einfall |
| `v4_al_sinus_te.json`, `v4_al_sinus_tm.json` | V4, Al-Sinusgitter |
| `v5_glass_deep_te.json` | V5, tiefes Glasgitter in Transmission |
| `echelette_au_scan.json` | Gold-Echelette (Dreieck), RK, Winkel-Scan |
| `v7_cobalt_anisotropic.json` | V7, anisotropes Co-Gitter |
| `v10_hole_array.json` | V10, Ag-Lochgitter (2D) |
| `nearfield_lamellar_1d.json` | Nahfeld-Beispiel 1D (horizontaler + vertikaler Schnitt) |
| `sweep_lamellar_depth.json` | Parameter-Sweep: Stegtiefe einer Glas-Lamelle 0,05–1,2 µm, TE und TM (Littrow-nah) |
| `nearfield_pillars_2d.json` | Nahfeld-Beispiel 2D: Si-Säulen, horizontale Karte (x, z) und vertikaler Schnitt (x, y) |
