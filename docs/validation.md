# Validierung

Sollwerte für die Referenzfälle (Plan, Abschnitt „Validierung und Tests“). „Intern“ bedeutet:
in der Prüfung der Buchkapitel mit unabhängigen Python-Implementierungen berechnet, nicht aus
der Literatur. Konfigurationen liegen in `examples/`.

| ID | Fall | Sollwert | Toleranz |
|---|---|---|---|
| V1 | ebene Grenzfläche, alle Winkel, TE/TM, beliebiger Azimut | Fresnel (analytisch) | 1e-12 |
| V2 | homogene Schicht auf Substrat (Airy), auch absorbierend | analytisch | 1e-12 |
| V3 | dielektrisches Lamellengitter, TE/TM/konisch | Energie = 1; RK = RCWA; Reziprozität (Abweichung fällt mit N, da die Ordnungsmengen beider Richtungen verschoben sind) | 1e-10 / 1e-8 / < 1e-4 bei N = 40 |
| V4 | Al-Sinus: d = 0,5 µm, g = (a/2)(1+cos Kx), a = 0,2 µm, λ = 0,6328 µm, Φ = 30°, ν = 1,38+7,62i, Ordnung −1 | TE 0,3235; TM 0,8696 (intern, C-Methode) | 2e-3 |
| V5 | Glas-Sinus: d = 0,5 µm, a = 1,0 µm, ν = 1,5, Φ = 30°, Transmission −1, TE | 0,820 (intern, 120 Scheiben; Buch Abb. III.4) | 5e-3 |
| V6 | dielektrischer Sinus ν = 2,5, a/d = 1, d = 1 µm, 30°, TM (Buch Abb. IV.3) | FFF konvergiert ab N ≈ 8 auf 1e-4, klassisch nicht | qualitativ |
| V7 (M5: 0,05673 / 0,14355 / 1,275e-7 / 2,372e-6) | Co-Sinus anisotrop: d = 0,6, a = 0,5 µm, λ = 0,6328 µm, 30°, TM-Einfall, ε_ii = −8,19+16,38i, ε_xz = −ε_zx = −0,495−0,106i | TM₀ 0,0567; TM₋₁ 0,1435; TE₀ 1,3e-7; TE₋₁ 2,37e-6 (intern, N = 20, 100 Scheiben) | 3 % |
| V8 | anisotroper Code, isotroper Tensor | Lamellen: = isotroper Code; Profil (Al-Sinus, N = 20): FFF_A = FFF_rig. 2D (M8): Laurent mit Tensor = isotrop (Polygon, Eierkarton); `FFF_B`-Streifen = Li; z-invariante anisotrope Strukturen = 1D | 1e-10 / 1e-3 |
| V7b | einachsiges Substrat (Achse ∥ y, ε_o = 2,4, ε_e = 3,1), TE/TM, alle Winkel | r_s = (β₀ − k_o)/(β₀ + k_o), r_p = (β₀ − k_e/ε_o)/(β₀ + k_e/ε_o) | 1e-12 |
| V7c | gedrehte biaxiale verlustfreie Schicht, konisch | Energie; homogener Pfad = Lamellenpfad (RCWA, RK) | 1e-12 / 1e-10 / 1e-8 |
| V9 | 2D-Code, z-invariante Struktur (Lamellen, Polygon-Streifen, `sinus_x`, 1D-Profil im 2D-Stapel), RK und RCWA | = 1D konisch (M7: N_z = 0 1e-12; RCWA N_z > 0 1e-10; RK N_z > 0 auf RK-Toleranz) | 1e-10 |
| V10 | Ag-Lochgitter: Film 0,2 µm auf Glas, Quadratlöcher 0,25 µm, Periode 0,9 µm, ν_Ag = 0,1+8,94i, Normaleinfall | T(1,388 µm) ≈ 7,5 % (Buch); SPP-Pol eben: λ = 1,369406 µm. M7: Test qualitativ (Peak oberhalb des Pols, ≥ 10× Untergrund, polarisationsunabhängig) | 10 % (offen) |
| V11 | Al-Sinus 40°, ν = 1,3+7,6i: Li-Lamellentreppe vs. FFF (Buch Abb. VI.10) | Treppe langsam, FFF schnell | qualitativ |
| V13 | Gold-Sinus TM: d = 1 µm, a = 0,25 µm, λ = 0,6328 µm, ν = 0,2+3,4i, Φ = 20°, Ordnung −1 und 0 | C-Methode (intern): e₋₁ = 0,26614, e₀ = 0,02362 | FFF_rig 1e-3, FFF_sym 3e-3; Buchvariante versagt |
| V14 | Gold-Dreieck TM (Spitze bei x = 0,7 µm), sonst wie V13 | **offen** – keine belastbare Referenz; geprüft: Σe < 1 | – |
| V15 | Säulen mit Kreisquerschnitt (n = 3, h = 0,2 µm, R = 0,14 µm, Periode 0,5 µm, N = 4): `contour2d` (exakter Kreis) gegen `polygon2d` mit L = 16…512 Ecken, Laurent | Fehler ∝ L⁻² (Verhältnis 4,00); zweistufige Richardson (L = 128, 256, 512) = Contour2D; RK = RCWA | 1e-8 (gemessen ≤ 1,2e-10); 1e-8 |
| V16a | Geneigte Wand ohne CAD: `SectionedLayer` (Kegelstumpf R 0,17 → 0,07 µm, h 0,2 µm; Kugelkappe ρ 0,25 µm, h 0,12 µm; n = 2 auf Glas) gegen `Profile2D` mit g(x, z), RK tol 1e-9, Laurent, N = 2 | gleiche Effizienzen, Energie verlustfrei | 1e-6 (gemessen ≤ 1e-15); 1e-10 |
| V19 | Gold-Lamellen TM (M17): ν = 0,2+3,4i, Periode 1 µm, Breite 0,5 µm, Höhe 0,3 µm auf Gold, λ = 0,6328 µm, Φ = 20° | ohne ASR bis N = 800 (Richardson in 1/N): e₋₁ = 0,013591, e₀ = 0,79435 (±5e-5); mit ASR gleicher Grenzwert | 1e-4 (η = 0,95, N = 40: 8e-5); Tabelle unten |
| V12 | Nahfeld: Airy-Schicht analytisch (RK, RCWA); Stetigkeit an allen Element-/Scheibengrenzen; Flusskonstanz in verlustfreien Lamellenstapeln; RK-Feld = RCWA-Feld; D_x-Stetigkeit an einer Lamellenwand | analytisch / Stetigkeit / Fluss = Σ T | 1e-12; 1e-6; 1e-12 (RCWA), 1e-8 (RK); 1e-7; FFF < 0,1 (direkt > 1) |

Zusätzliche Kontrollen aus der Kapitelprüfung (Al-Sinus V4, TM, FFF-Treppe):
60 / 150 / 300 Scheiben → 0,853 / 0,865 / 0,869. Treppe mit Li-Lamellenformel bei N = 30:
M = 1 → 0,565, M = 20 → 0,804, M = 200 → 0,774 (keine monotone Konvergenz).

**Energie bei FFF:** Keine FFF-Variante erhält bei verlustfreien Medien im TM-Fall die Energie exakt
(Abweichung 1e-5 bis 5e-4, Laurent exakt). Tests prüfen die Energie bei Profilen deshalb mit Toleranz 5e-4.

**Ergebnisse M3 (RK, Gauß 2, tol 1e-5):** V4 TM N = 20: FFF_iso 0,86951, FFF_rig 0,87075, FFF_sym 0,86565;
N = 30 (FFF_iso): 0,87000. TE N = 20: 0,32335. RCWA 300 Scheiben N = 20: FFF_iso 0,86866, FFF_rig 0,86996.
V5 TE: 0,81996 (RCWA), 0,81997 (RK).

**Nahfeld in FFF-Profilschichten:** Der Leistungsfluss ist innerhalb einer FFF-Profilschicht nur auf etwa
3e-4 konstant (dieselbe Ursache wie die nicht exakte Energiebilanz der FFF); in homogenen Schichten und
Lamellen exakt.

**V10 (M7, RCWA, `FFF_sym`, NF-Tensor, σ = 3, G = 512):** Transmissionspeak (golden-section-Suche)

| N | λ_Peak (µm) | T_max | T(1,388 µm) | Rechenzeit je λ (1 Kern) |
|---|---|---|---|---|
| 5 | 1,3828 | 4,1 % | 1,2 % | 2 s |
| 6 | 1,3860 | 3,4 % | 2,4 % | 4,5 s |
| 7 | 1,3906 | 10,3 % | 3,9 % | 10 s |
| 8 | 1,404 | 5,3 % | 0,3 % | 21 s |
| 9 | 1,3876 | 4,2 % | 4,2 % | 55 s |
| 10 | 1,3889 | 6,2 % | 5,4 % | 80 s |
| 11 | ≈ 1,3898 | ≈ 9,2 % | 4,9 % | 140 s |

Der Peak liegt stets oberhalb des SPP-Pols der Ag/Glas-Grenzfläche (1,3694 µm). Bis N = 8 wandert er nach Rot
und seine Höhe schwankt; ab N = 9 liegt er bei der Buch-Wellenlänge 1,388 µm (1,3876 / 1,3889 µm), die Höhe
steigt auf 6,2 % (N = 10) bzw. ≈ 9,2 % (N = 11, Suche nicht ganz abgeschlossen): die Höhe schwankt um die
Buchangabe 7,5 %, konvergiert aber noch nicht. Der Peak ist nur ≈ 2–3 nm breit, daher reagiert
T bei fester Wellenlänge stark auf N. Vektorvariante von Algorithmus NF (N = 7: 1,3904 µm / 10,26 %;
N = 8: 1,4024 µm / 5,53 %) weicht kaum von der Tensorvariante ab (1,3906 / 10,32 %; 1,4038 / 5,35 %) – die
N-Abhängigkeit dominiert. Laurent schwankt
stärker, `FFF_iso` liefert R > 1. Offene Konvergenzstudie (N = 10–15, Tensor/Vektor, σ) auf dem Zielsystem.

**Konvergenzstudie M8 (RCWA, λ = 0,6328 µm, Φ = 30°, Ψ = 20°, δ = 60°, Nx = Nz = N):**

Co-Säulen (Tensor V7, Quadrat 0,3 µm, Periode 0,6 × 0,6 µm, Höhe 0,1 µm) auf Glas, Σ T:

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| `FFF_B` | 0,5485 | 0,5431 | 0,5487 | 0,5495 | 0,5493 | 0,5508 |
| Laurent | 0,4165 | 0,4600 | 0,4744 | 0,4765 | 0,4963 | 0,4911 |

Anisotroper Eierkarton (Tiefe 0,25 µm, Periode 0,6 µm, einachsig 2,75/2,21 gedreht, verlustfrei), r(0,0), 20 Scheiben:

| N | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| `FFF_A` | 0,006098 | 0,006033 | 0,006022 | 0,006010 | 0,006003 |
| `FFF_B` | 0,006089 | 0,006090 | 0,006021 | 0,006006 | – |
| Laurent | 0,005605 | 0,005688 | 0,005742 | 0,005780 | – |

Energiefehler der FFF-Varianten 1e-3 (N = 2) bis 3e-5 (N = 6). Laufzeit (1 Kern) N = 5: Säulen 3 s, Eierkarton 55 s.

**V12 (2D, M9):** z-invariante Strukturen (Polygon-Streifen, `sinus_x`) = 1D-Nahfeld (1e-9; horizontale Karten =
1D-Linie · e^{iγ0(z−z0)}), Stetigkeit von F an allen Element- und Scheibengrenzen (RK, RCWA) 1e-6, Flusskonstanz
(verlustfrei) 1e-12, RK = RCWA im Säulenfeld 1e-7. D_x-Sprung an der Wand einer Säule (ε = 4, N = 6, Raster 1 nm):
direkt 293 %, FFF-konsistent 5,3 % (isotrop, Tensorfeld) bzw. 3,0 % (anisotrop, `FFF_B`).

## Genauigkeitsstudie M16 (fp64 / mixed / fp32)

`tests/validation/test_precision_study.cpp`: alle Beispiele mit Golden File (V1, V3, V4, V5, V7, V10, V15, Nahfeld,
Echelette, Sweep) und die Profilfälle V9, V13 (`bench/`), je Modus gegen `fp64`. Angegeben ist die größte
Abweichung einer Effizienz (Anteil der einfallenden Leistung) über alle Ordnungen und Scanpunkte, die
Energiebilanz |Σe − 1| der verlustfreien Punkte in `fp32` und ob der Rückfall auf fp64 ausgelöst wurde.
CPU: Eigen in einfacher Genauigkeit (`FG_BACKEND=cpu`); GPU: RTX 3090, alles auf dem Gerät
(`FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0`). V16–V18 (CAD) sind im MSVC-Build ohne Slicer nicht enthalten.

| Fall | mixed CPU | fp32 CPU | mixed GPU | fp32 GPU | |Σe−1| fp32 GPU | Rückfall |
|---|---|---|---|---|---|---|
| V1 Fresnel | 7,9e-10 | 4,5e-8 | 7,9e-10 | 5,0e-9 | < 1e-6 | nein |
| V3 Lamellen RCWA | 4,4e-6 | 4,5e-6 | 1,4e-6 | 3,4e-6 | 3e-6 | nein |
| V3 Lamellen RK (Scan) | 1,6e-6 | 3,3e-6 | 1,6e-6 | 2,5e-6 | 2e-6 | nein |
| V4 Al-Sinus TE | 4,2e-7 | 3,5e-7 | 5,6e-7 | 2,7e-7 | – | nein |
| V4 Al-Sinus TM | 1,4e-7 | 3,3e-7 | 1,6e-6 | 1,2e-6 | – | nein |
| V5 Glas tief TE | 1,2e-6 | 4,1e-6 | 1,9e-6 | 1,8e-5 | 1,8e-5 | nein |
| V7 Kobalt anisotrop | 1,3e-6 | 1,3e-6 | 3,3e-7 | 6,1e-8 | – | nein |
| V9 `sinus_x` 2D RK | 6,0e-7 | 1,7e-6 | 1,3e-6 | 9,4e-7 | – | nein |
| V10 Ag-Lochgitter N = 10 | 1,1e-4 | 1,1e-4 | 4,6e-5 | 1,1e-5 | – | nein |
| V13 Gold-Sinus RK | 1,2e-7 | 3,6e-7 | 1,2e-7 | 2,7e-7 | – | nein |
| V15 Kreissäulen | 3,8e-5 | 3,8e-5 | 2,6e-5 | 2,6e-5 | – | nein |
| Echelette Au (Scan) | 2,7e-6 | 7,1e-6 | 3,1e-6 | 3,8e-6 | – | nein |
| Nahfeld 1D / 2D (Effizienzen) | 6,2e-6 / 1,7e-5 | 6,2e-6 / 1,7e-5 | 6,2e-6 / 1,4e-5 | 1,9e-5 / 1,0e-5 | – | nein |
| Sweep Lamellentiefe | 5,8e-6 | 5,9e-6 | 1,1e-5 | 1,1e-5 | 1e-5 | nein |

Befund: `mixed` ist kaum genauer als `fp32`, weil die Eigenprobleme (in beiden Modi einfach genau) den Fehler
bestimmen, nicht die LU-Zerlegungen. Festgeschriebene Zielwerte (Plan v2.8): `mixed` und `fp32` ≤ 1e-4 für V1–V9
(gemessen ≤ 2e-5); resonante 2D-Fälle (V10, V15) bis ≈ 1e-4. Für Ergebnisse auf 1e-6 und besser: `fp64`.

## Koordinatentransformationen und ASR (M17)

`tests/validation/test_v19_asr.cpp`. η = 0 ist bitgleich zum Pfad ohne Transformation. Der transformierte Pfad
(dichte A, G, numerische Halbraummoden, Umrechnung eq:ct_inc / eq:ct_rayleigh) mit f′ = g′ = 1 stimmt mit dem
Diagonalpfad überein (V2, V3 dielektrisch/metallisch, RK, V9-Polygonstreifen 2D, Sinusprofil: ≤ 1e-12, RK 1e-10).

| Fall (η = 0,5) | Ergebnis |
|---|---|
| V1/V2 eben, eine Bruchstelle | = analytisch ≤ 1e-12 ab N = 6 (f analytisch) |
| V1/V2 eben, zwei Bruchstellen | 4,8e-5 (N = 4), 2,0e-7 (8), 5,3e-9 (16), 2,0e-10 (32): f‴ springt |
| V3 dielektrisch, Abweichung von N = 120 ohne ASR | N = 10/20/40: 4,1e-5 / 6,5e-6 / 1,1e-6 (ohne ASR 3,4e-4 / 2,9e-5 / 4,2e-6) |
| V3 dielektrisch, Energiebilanz | −9e-6 (N = 10), −1,4e-7 (20), −1,2e-8 (40): Umrechnung der Ordnungen |
| V3 metallisch (0,5+3i) | N = 40: 1,9e-4 (ohne ASR 4,8e-4) |
| V9 z-invariant (Polygonstreifen, sinus_x, 1D-Profil im 2D-Stapel) | = 1D konisch mit gleicher Transformation ≤ 1e-10 (auch mit z-Transformation) |
| Nahfeld Airy-Schicht | = ohne Transformation ≤ 1e-10 |

Die Abweichung vom untransformierten Ergebnis bei gleichem N ist der Abschneidefehler beider Rechnungen; der
Planwert 1e-8 gilt nur für die ebenen Fälle (V1, V2). Für V3/V9 prüft der Test den gemeinsamen Grenzwert und die
Übereinstimmung der z-invarianten 2D-Rechnung mit 1D.

**V19 (Gold-Lamellen, TM).** e₀ und Abweichung vom Grenzwert 0,79435:

| N | ohne ASR | ASR η = 0,8 | ASR η = 0,95 |
|---|---|---|---|
| 10 | – | 0,79305 (1,3e-3) | 0,79422 (1,3e-4) |
| 20 | – | 0,79378 (5,7e-4) | 0,79421 (1,4e-4) |
| 25 | 0,78713 (7,2e-3) | – | – |
| 40 | – | 0,79408 (2,7e-4) | 0,79427 (8e-5) |
| 50 | 0,79322 (1,1e-3) | – | – |
| 80 | – | 0,79420 (1,6e-4) | 0,79429 (5,6e-5) |
| 100 | 0,79365 (7,0e-4) | – | – |
| 160 | – | 0,79422 (1,3e-4) | 0,79430 (5e-5) |
| 200 / 400 / 800 | 0,79382 / 0,79406 / 0,79421 | – | – |

Benötigtes N für 1e-3: ohne ASR ≈ 60, η = 0,8 ≈ 15, η = 0,95 ≤ 10. Der Rest von ≈ 5e-5 zwischen ASR (η = 0,95) und
der Richardson-Extrapolation liegt in der Unsicherheit der Extrapolation (ohne ASR konvergiert e₀ nur wie 1/N).
Literaturvergleich (Granet 1999) nicht durchgeführt: die Daten sind bibliographisch nicht geprüft.
