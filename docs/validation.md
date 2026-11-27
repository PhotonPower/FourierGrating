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

## PML, isolierte Strukturen und Streufeld (M18)

`tests/validation/test_v20_pml.cpp` (`test_validation_v20`), CAD-Vergleich in `tests/cad/test_cad.cpp` (`CAD_V20g`).
λ = 0,6328 µm, Luft über Glas (n = 1,5); Steg n = 1,5 (0,5 × 0,3 µm), Graben 0,5 µm in einem Film n = 2,1 (0,2 µm),
Kontaktloch r = 0,25 µm in einem Film n = 2,1 (0,3 µm). N ≈ 12 L/λ.

| Test | Ergebnis |
|---|---|
| V20a Film mit PML, Normaleinfall | Feldkoeffizienten = Airy (N = 0) ≤ 1e-10, Streubreiten 0 |
| V20b Steg, Zellbreite 4/6/8 × Stegbreite (w = 0,15 L) | TE oben 0,00909 / 0,01030 / 0,01067, unten 0,9072 / 0,9101 / 0,9104; TM oben 0,00515 / 0,00516 / 0,00517, unten 0,7788 / 0,7790 / 0,7791 |
| V20b w, β bei festem Bereich (2,4 µm) | Streubreiten ändern sich < 2e-5 |
| V20b Absorption (verlustfrei) | ≤ 2e-4 (N = 50…150; N = 200: 3e-4, Kondition der Eigenzerlegung mit komplexem A) |
| V20c Wellenleiter → gleicher Wellenleiter | R = 0, T = 1 ± 3e-6 (TE0, TM0) |
| V20c abruptes Ende (Kern n = 2,0, 0,4 µm, TE0, n_eff = 1,9047) | R(TE0) = 0,15999, R(TE2, n_eff 1,067) = 0,0191, abgestrahlt 0,8210; N = 30/50/80, Zellbreite, w, β: ± 1e-4 |
| V20d Streufeld = Totalfeld (periodisch, ohne PML) | Lamellen TE/TM/konisch: RCWA 1e-15, RK 3e-12; Metallprofil, Contour2D-Säulen ≤ 1e-10; Nahfeld ≤ 1e-10 |
| V20e Normaleinfall mit PML | Streufeldformulierung = Totalfeld − Hintergrund: alle Ausgaben gleich (1e-10) |
| V20f Gesamtstreuung, Fenster 2,4 / 3,2 / 4,0 / 4,8 µm | Steg TE 60°: 1,55441 / 1,55433 / 1,55433 / 1,55435; Steg TM 60°: 0,97681 / 0,97722 / 0,97700 / 0,97702; Graben TM 30°: 1,31077 / 1,31098 / 1,31103 / 1,31107; Graben TE 0°: 1,44502 / 1,44512 / 1,44508 / 1,44514 |
| V20f oben/unten | konvergieren nur ~1/L (Steg TE 60° unten 1,1801 / 1,1910 / 1,1979 / 1,2006): streifende Streuung kreuzt die Ebenen außerhalb des Fensters |
| V20f w, β | alle Ausgaben < 5e-4; Absorption ≤ 7,6e-4 (Graben TM) |
| V20g z-invariantes 2D (PML in x) | = dz × 1D ≤ 1e-10 (normal und 30°) |
| V20g Kontaktloch, x- und z-Polarisation | gleich (Symmetrie) ≤ 1e-10; CAD-Körper = Contour2D ≤ 1e-10 |
| V20g Konvergenz in N (Zelle 2 µm) | Gesamtstreuung 0,501 / 0,568 / 0,570 bei N = 4/6/8, Absorption 3e-2 / 1e-2 / 1e-3 |
| V20f Reziprozität f(α_s; α_i) = f(−α_i; −α_s), Film und Steg nebeneinander (asymmetrisch), 20° → −35° | komplexe Abweichung TE 5,1e-2, TM 2,5e-2 bei L = 8 µm (L = 4 / 6 µm: 0,15–0,41 / 0,03–0,14); Fernfeld aus dem Streufeld auf einer Linie über der Struktur (Winkelspektrum über das Fenster, cos²-Fenster), Genauigkeit durch die endliche Linie begrenzt; Test-Toleranz 0,1 |

Offen: Literaturvergleich V20c (Silberstein et al. 2001, Wert nicht bibliographisch geprüft), Vergleich der PML mit
Hugonin/Lalanne (2005, Formulierung nicht geprüft), Konvergenzstudie V20g in Zellgröße und N auf dem Zielsystem.

## Gradienten (M19)

`tests/unit/test_adjoint.cpp` (jede Rückwärtsregel gegen zentrale Differenzenquotienten, Richtungsableitung mit
Zufallsmatrizen, h = 1e-7, 1e-6), `tests/validation/test_gradient.cpp` (`test_validation_gradient`), CAD in
`tests/cad/test_cad.cpp` (`CAD_M19`), Python `tests/python/test_gradient.py`. Angegeben ist die relative
Abweichung Adjungierte gegen Differenzenquotient der Vorwärtsrechnung (h = 1e-6 relativ).

| Fall | Parameter | Abweichung |
|---|---|---|
| V3-Lamellen, konisch, Li und Laurent | Dicken, Lamellenkanten, n (Re, Im), ε eines Films | ≤ 6,6e-8 (Kanten), sonst ≤ 1e-8 |
| Symmetrisches Lamellengitter, Normaleinfall (exakte Entartung) | Dicke, n, Kante | ≤ 1e-9 |
| Zentrierte quadratische Säule, Normaleinfall (entartete Paare, von der Störung gekoppelt; Plan v2.53) | zwei Eckpunkte, Dicke, n | ≤ 3,7e-9 (vorher Eckpunkte 15 % bzw. 9 %) |
| V9/V15: Polygon- und Kreissäulen (2D, Laurent) | Ecken, Radius, Mittelpunkt, Dicke, n | ≤ 2,5e-9 |
| NURBS-Säule (rationaler Kreis, Laurent) | Kontrollpunkte / Gewicht | ≤ 9e-9 / 4,5e-8 |
| Kreissäulen, FFF_sym (v2.11: Normalenfeld NF eingefroren) | Radius / n | 1,8e-7 / 2e-9 (N = 3); Mittelpunkt: −0,44 gegen Trend ≈ 0 (46 % bei N = 3, 150 % bei N = 5, 29 % bei N = 8), Polygonecke 0,7 % |
| Säulen FFF_sym, v2.12: Rasterfeld NF über eine Rasterzelle mitbewegt | Radius, n / Mittelpunkt, Polygonecken | 2e-7, 2e-9 / Mittelpunkt −5e-6 (Trend −7e-4 bei h = 2e-3), Ecken 4e-4–6e-4 gegen den Quotienten bei h = 1e-3 |
| Anisotrope 2D-Säulen, FFF_B | Tensoreinträge, Dicke, Radius / Mittelpunkt, Ecken | ≤ 4e-7 / Mittelpunkt 2e-6 (Trend −6e-5), Ecken ≤ 2,5e-3 gegen h = 1e-3 |
| V13-Typ Gold-Sinus TM, Laurent | Tiefe, Verschiebung | 3e-10, Verschiebung ≈ 0 (Translation) |
| V13-Typ, FFF_rig (Profil-Normalenfeld mitdifferenziert) | Tiefe | 7e-7 (eingefroren wäre der Gradient −7,0 statt −83,0) |
| RK (Lamellen, Sinusprofil) | Dicke, Kante, n | ≤ 9e-7 (der Quotient sieht andere adaptive Schritte); RK = RCWA ≤ 2e-9 |
| V18-Typ Kegelstumpf, Schnittschicht (analytisch) / CAD-Körper | oberer Radius | 1e-9 / 8e-9 |
| Phase 2: konische Lamellen, Gold-Sinus (verlustbehaftetes Substrat), 2D-Säulen | λ, φ, ψ, δ, n Superstrat, n Substrat (Re, Im) | ≤ 3e-8, λ der 2D-Säulen 1,8e-7 |
| Phase 2: anisotrope Lamellen (Laurent, FFF_A, FFF_B) | ε-Tensoreinträge (Re, Im), Kanten, isotropes Material derselben Schicht, Dicke, λ; homogener Film: Hauptwerte, Euler-Winkel | ≤ 2e-8, Euler-Winkel bei Laurent 8e-7 (Gradient 1,5e-6) |
| Phase 2: anisotropes Sinusprofil (FFF_A, FFF_B) | Tiefe, Verschiebung, Euler-Winkel, Hauptwert | ≤ 3e-7 (Verschiebung ≈ 0: Translation) |
| Phase 2: anisotrope 2D-Säule (Laurent) | Radius, Tensoreintrag | 4e-11 / 4e-8 |
| Phase 2: RK (Lamellen) | λ, φ, n Superstrat/Substrat | ≤ 1,7e-6 (der Quotient sieht andere adaptive Schritte) |
| Phase 2: ASR (η = 0,5) konische Lamellen (Li, Laurent) | Dicken, Kanten (auf und neben den Fixpunkten), n (Re, Im), ε | ≤ 4e-7; Wirkungsgrade = Löser 1e-12 |
| Phase 2: ASR Sinusprofil (Laurent, FFF_rig), anisotrop (FFF_A, FFF_B) | Tiefe, Verschiebung, n, Euler-Winkel | ≤ 1,2e-7, Euler-Winkel 6e-7 (Gradient 3e-5) |
| Phase 2: ASR anisotrope Lamellen (Laurent, FFF_A, FFF_B) | Tensoreintrag, Kanten, Dicke | ≤ 7e-9 |
| Phase 2: ASR in x und z, 2D-Säulen (Laurent) | Radius, Mittelpunkt, Ecken, n, Dicke, Tensoreintrag | ≤ 3,4e-7 (abgebildete Kontur per Differenzenquotient, adaptive Quadratur) |
| Phase 2: ASR mit RK (Lamellen) | Dicke, Kante, n, λ | ≤ 8e-6 relativ (3,7e-8 absolut; adaptive Schritte) |
| Phase 2: ASR, globale Parameter (Lamellen, Sinusprofil, 2D-Säulen in x und z) | λ, φ, ψ, δ, n Superstrat, n Substrat (Re, Im) | ≤ 2,5e-7 (λ der 2D-Säule), sonst ≤ 2e-8 |
| Phase 2: `profile2d` Eierkarton (Laurent, FFF_rig; anisotrop FFF_A, FFF_B), sinus_x anisotrop | Tiefe, n (Re, Im), λ, Hauptwert, Euler-Winkel | ≤ 3e-7 (Euler-Winkel, Gradient 1e-4), sonst ≤ 1,5e-8 |
| Phase 2: `profile2d` Verschiebungen (zwei Eierkartons übereinander; Laurent, FFF_rig, anisotrop FFF_A, FFF_B) | x0, z0 der oberen, x0 der unteren Schicht | ≤ 1,1e-8 (Gradienten beider Schichten entgegengesetzt gleich) |
| v2.13 Stufe A: Steg mit PML (1D, Laurent, Li, RK), Normaleinfall, P_oben, P_unten | Dicken (Steg, Film), Kanten, n (Re, Im), eps Film, λ, ψ, δ, n Superstrat, n Substrat | ≤ 3,3e-8 (RK 2,8e-8); Wert = Löser auf 1e-12 |
| v2.13 Stufe A: Kontaktloch (2D, PML in x und z, Laurent) | Radius, Mittelpunkt, n Film (Re, Im), Dicke | ≤ 2,8e-7 (Mittelpunkt, Gradient 4e-3), sonst ≤ 1,2e-8 |
| v2.13 Stufe B: Steg mit PML, schräg (φ = 25°, Laurent, Li) | Kanten, n Steg (Re, Im), eps Film, Dicken, λ, φ, ψ, δ, n Superstrat, n Substrat | ≤ 2,1e-8; Wert = Löser auf 1e-12 |
| v2.13 Stufe B: Kontaktloch 2D schräg in x und z; Steg mit `scattered_field` bei Normaleinfall | Radius, n Film (Re, Im), Dicke, λ, φ; Kante, n, eps, Dicke, λ | ≤ 1,4e-7 (φ, Gradient 1e-4), sonst ≤ 1e-7 |
| v2.13 Stufe B: Steg schräg mit RK (tol 1e-8) | Kante, eps Film, Dicke, λ | ≤ 2e-8 |
| v2.55 Stufe C1: Steg mit PML, Normaleinfall (N_x = 12), seitlicher Fluss / Absorption | Dicken, Kanten, n Steg (Re, Im), eps Film, λ, ψ, δ, n Superstrat, n Substrat | seitlich ≤ 6,6e-9, absorbiert ≤ 1,7e-7 (ψ 1e-6 relativ, 7e-10 absolut); Werte = Löser 1e-12 |
| v2.55 Stufe C1: Kontaktloch 2D (N = 3), seitlich + 0,5 absorbiert | Radius, Mittelpunkt, n Film (Re, Im), Dicke, λ | ≤ 5,2e-7 |
| v2.54: Profilhöcker mit PML (flach bei g0 = 0,06 µm), Normaleinfall und schräg (φ = 25°), RCWA und RK | n Profil (Re, Im), eps Film, Dicke, λ, φ | ≤ 2,2e-8; Wert = Löser auf 1e-12 |
| v2.13 Stufe A: modaler Einfall (Wellenleiter, Block, Wellenleiter als Substrat), \|r_0\|², \|t_0\|², \|r_1\|² | Dicke, Kante, n Block (Re, Im), n Kern (beide Halbräume), λ | ≤ 1e-8; Wert = Löser auf 1e-12 |

Aufwand (2D, N = 3, Laurent, 5 Parameter): Gradient = 1,3 × Vorwärtsrechnung (gemessen 1,3–2,3); die Kosten
hängen nicht von der Zahl der Parameter ab (Formableitungen kosten nur die Konturquadratur). Beispiele:
`optimize_littrow.py` (L-BFGS-B erreicht das Rasterscan-Optimum), `optimize_nurbs_pillar.py` (Adam über `fouriergrating.torch`, T(−1, 0) steigt in 25 Schritten monoton
von 0,178 auf 0,337; als Test in der Python-CI, die PyTorch, JAX und scipy installiert).
`checkpoint_every` = 2, 3, 7: Gradient bitgleich zu 1 (RCWA mit Scheiben und RK); gespeichert werden nur jedes k-te
Präfixprodukt, die Scheibendaten werden rückwärts je Abschnitt neu gerechnet (Aufwand: eine zusätzliche
Vorwärtsrechnung der Scheiben). Phase 2 (λ, Winkel, Superstrat/Substrat): Kotangenten von α, γ aus allen M und den Basen der Halbräume, von k₀ aus
den Dicken; Einfallsvektor und Effizienznormierung (wenige Skalare, O(n_f)) per lokalem zentralem
Differenzenquotienten. Anisotrope Scheiben: Rückwärtsregel der Faktorisierung (Laurent, FFF_A, FFF_B) bis zu den
Gittern der Koeffizientenfunktionen; deren Ableitung nach den Regionstensoren per zentralem Differenzenquotienten der
1D-Koeffizienten (holomorph, ohne O(n³)), Kanten analytisch, Profilverschiebung als Translation, Profiltiefe über das
mitbewegte Normalenfeld. Koordinatentransformationen (ASR): M mit A, G (Rückwärtsregel mit festen A, G),
numerische Basen der Halbräume und Einfall in X fest; die Wirkungsgrade aus den physikalischen Rayleigh-Amplituden
(lineare Abbildung, Kotangente durch ihre Adjungierte); Kanten und Profilkreuzungen in X mit dX/dx = 1/f'(X);
2D-Konturen über die abgebildeten Koeffizienten. Mit `"auto"`-Fixpunkten hängt die Abbildung von den Kanten und
Polygonecken ab: dafür Fehlermeldung (Fixpunkte explizit angeben). λ, Winkel und Superstrat/Substrat mit
Transformation: numerische Basen der Halbräume über die Eigenvektor-Rückwärtsregel (TE/TM-Paare entartet: nur der
Unterraum), Kotangenten von A, G aus allen M, Umrechnungsmatrizen nach k0 α0, k0 γ0 per zentralem
Differenzenquotienten der 1D-Matrizen. `profile2d`: isotrop als eps = eps_unten χ + eps_oben (1 − χ) mit dem
analytischen Normalenfeld, anisotrop über die Koeffizientengitter; die Tiefe hält die relative Höhe der Scheiben (gleiche
Kreuzungen und Quadraturlinien), es bewegt sich das Normalenfeld; Verschiebungen x0/z0 als Translation aller Gitter
(g_pq → −i (p Kx | q Kz) g_pq; ohne Transformation in dieser Richtung). Offen: Rückwärtsfunktionen über das
GPU-Backend, aperiodische Strukturen mit schrägem oder modalem Einfall. Aperiodisch, Normaleinfall (Plan v2.13,
Stufe A): drei Durchläufe desselben Treibers (Hintergrundstapel vorwärts, Stapel mit den Kotangenten der Streuflüsse,
Hintergrundstapel mit den negierten); Geometrieparameter wirken nur auf den Stapel, Dicken, Materialien und globale
Parameter auch auf den Hintergrund; S_inc und die numerischen Halbraumbasen mit PML über die
Eigenvektor-Rückwärtsregel. Modaler Einfall: Leistungen |b_k|² P_k / P_in, unabhängig von Phase und Skalierung der
Eigenvektoren; die Basen der strukturierten Halbräume über die Eigenvektor-Rückwärtsregel bis zu ihren Materialien.
Streufeldformulierung (Stufe B): Sternprodukt mit Quellvektoren rückwärts (`adjoint::starSources`), partikuläre Lösung
je Scheibe und Hintergrundexponent ((κ − M)⁻¹ σ), Quelle σ = (M − M_bg) F_bg nach M, M in der Blochbasis und M_bg;
Kotangenten der Hintergrunddaten (Amplituden, Exponenten, Modenvektoren, Schichtgrenzen) analytisch, ihre Ableitung
nach Dicken, Hintergrundmaterialien, λ und Winkeln per zentralem Differenzenquotienten der N = 0-Hintergrundrechnung
(ebenso Tsrc und I_0); Scheibengrenzen analytisch. RK: erweitertes System (F; 1) mit der Quelle σ(y) in der letzten
Spalte, Rückwärtsregel des RK-Schritts unverändert; der Löser teilt die RK-Scheiben zusätzlich an den Knoten der
seitlichen Flussquadratur (andere Schrittfolge), der Wert stimmt daher nur auf die RK-Toleranz (7e-10). Profilschichten
mit PML (v2.54): Scheiben, die die Profilhöhe der PML-Zone enthalten, werden dort geteilt (gleiches Q, je Stück die
Hintergrundschicht nach der Höhe), wie im Löser. Seitlicher Fluss und Absorption (v2.55, Stufe C1, Totalfeld): die
Kette wird an den Seitenknoten geschnitten, die Felder dort aus den Amplituden zwischen Präfix und Suffix
(`interfaceAmplitudes`); rückwärts zusätzlich eine Schleife über die Suffixe von unten nach oben. Rasterfeld von Algorithmus NF (v2.12): unterhalb
einer Rasterzelle springt die Vorwärtsrechnung um ~1e-7 (Kreismittelpunkt in 2-nm-Schritten), der volle
Differenzenquotient hängt dort von h ab (−0,82 / −0,086 / −0,0077 / −0,0012 für h = 1e-6 … 1e-3); eingefroren erzeugt
das stehenbleibende Feld eine Scheinableitung. Formparameter mit FFF bewegen das Feld deshalb mit
(Differenzenquotient der NF-Koeffizienten bzw. der FFF_B-Gitter über eine Rasterzelle, zwei NF-Rechnungen je
Parameter); der Indikatoranteil bleibt exakt. Anisotrope FFF_B-Querschnitte: Materialparameter als
Richtungsableitung der Gitter (zwei Auswertungen je Parameter).

## Zylinderkoordinaten (M20)

`tests/unit/test_cyl.cpp` (Bausteine), `tests/validation/test_v21_cyl.cpp` (`test_validation_v21`).

| Prüfung | Ergebnis | Toleranz |
|---|---|---|
| Bessel J_n(z), reell (n ≤ 40, x ≤ 120) gegen `std::cyl_bessel_j`; komplex und \|z\| ≈ 500 gegen scipy/mpmath | ≤ 1e-12 relativ / ≤ 1e-13 bezogen auf max(\|J\|, √(2/π\|z\|)) | 1e-13 |
| Lommel-Integrale (ν = 0, 1, 3, −2; ν = −ν′) gegen Gauß-Quadratur, reeller und deformierter Weg | ≤ 2,2e-14 | 5e-14 |
| Hankel-Paar r^ν e^(−r²) (ν = 0, 1, 3), diskrete Rücktransformation | ≤ 1e-12 | 1e-12 |
| Galerkin-Matrix wirkt als Multiplikation (Hintergrund analytisch, gleiche und gemischte Ordnungen) | ≤ 1e-12 | 1e-12 |
| Li in r mit ε ≡ ε∞; homogenes Medium: M² = k² − ε, Ψ Eigenvektoren; Fluss geschlossen = eq:cyl_flux | 1e-15 / 1e-13 / 1e-15 | – |
| Wurzelsubstitution ∫₀^{k_h} √(k_h² − k²) dk; deformierter Weg = reeller Weg (analytischer Integrand), 16 Knoten | ≤ 1e-13 | 1e-13 |

**V21b (Fasermoden).** Kern ε = 4, Radius k₀a = 2, Mantel ε = 1; Fehler \|β − β_exakt\| des nächstgelegenen
Eigenwerts von M (Panelbreite π/(2a), 8 Knoten):

| Mode | β exakt | Li in r, k_max = 8 / 16 | Laurent, k_max = 8 / 16 |
|---|---|---|---|
| TM01 (m = 0) | 1,197891848750 | 5,6e-4 / 7,7e-5 | 1,6e-2 / 8,1e-3 |
| HE11 (m = 1) | 1,711581418489 | 3,6e-4 / 1,9e-5 | 3,1e-3 / 1,6e-3 |
| m = 2 | 1,152981421773 | 4,9e-4 / 1,1e-4 | 1,1e-2 / 5,4e-3 |

Die Werte bei k_max = 16 stimmen mit dem Prototyp (Plan, Befund) überein; Laurent ist erster Ordnung, TE01 ist für
beide Faktorisierungen gleich (4e-5).

**V21a (Strahl durch Schichtstapel).** Luft / Film n = 1,45 (0,3 µm) / Film n = 1,38 bzw. 1,38 + 0,05i (0,2 µm) / Glas,
λ = 1 µm, Taille 0,5 µm bei y = 0,8 µm, κ_max = 2: R, T gegen die Dünnschichtlösung je Knoten (kartesischer Löser,
N = 0, sin φ = k_j, s und p gewichtet mit der einfallenden Leistung) für x-, radiale und azimutale Polarisation:
Abweichung ≤ 2e-15, verlustfrei R + T = 1 auf 1e-15. Ringschicht mit einem Material = homogene Schicht (1e-12).

Dielektrische Säule (n = 1,45, r = 0,4 µm, h = 0,3 µm, Luft, auf Glas), x-Polarisation, RCWA:

| Faktorisierung | κ_max = 2 (N_k = 96) | 4 (192) | 8 (384) | R + T − 1 (2 / 4 / 8) |
|---|---|---|---|---|
| Laurent | R = 0,049599 | 0,048994 | 0,048749 | 7e-15 / 9e-14 / −2e-13 |
| Li in r | R = 0,048787 | 0,048511 | 0,048496 | 7,3e-4 / −1,7e-5 / −1,5e-5 |

dk_factor = 2 ändert R nicht in den ersten 8 Stellen. m = −1 und m = +1 tragen gleiche Leistung (1e-12); RK = RCWA
(1e-8).

**Ebene Wellen (Streufeld).** Spalten der Galerkin-Matrizen am Zusatzknoten k′ = k_∥ (Lommel, gemischte Ordnungen per
Quadratur) gegen Quadratur ≤ 5e-14; Quellspektren nach Jacobi–Anger gegen numerische φ-Integration der ebenen Welle
≤ 1e-13; konsistente Quelle b̂ − Q â = (Q − ε∞) E_bg für k_∥ auf einem Knoten ≤ 1e-12. Säule wie oben (κ_max = 2):
Normaleinfall δ = 0 und 90° gleich (1e-12), schräger Einfall (φ = 30°) je Harmonischer unabhängig vom Azimut ψ
(1e-10), Ringe aus Hintergrundmaterial streuen nicht, RK mit Quellterm = RCWA (1e-8, Laurent und Li).

**V21c (Mie-Kugel).** n = 1,5, k₀a = 3 in Luft, λ = 1 µm, Treppe aus N Ringschichten (Radius in Scheibenmitte),
Li in r; σ_sca in µm² (Mie: 2,448002):

| N | κ_max = 4 | 8 | 12 | exakte Quelle, κ_max = 4 / 8 |
|---|---|---|---|---|
| 10 | 2,475229 | 2,480102 | 2,480778 | 2,536589 / 2,510417 |
| 20 | 2,448735 | 2,455710 | 2,457501 | 2,509851 / 2,485790 |
| 40 | 2,441324 | 2,447471 | – | 2,501089 / 2,476878 |
| 80 | 2,439120 | 2,444564 | – | 2,498970 / 2,473872 |

Extrapolation im Test (Ordnung in N aus κ_max = 4: p = 1,84; Richardson bei κ_max = 8; Rest κ^−2): 2,446596, also
−5,7e-4 (Toleranz 1e-3). Mit der exakten Quelle (P = (ε − ε∞) E_bg) konvergiert σ_sca nur erster Ordnung in κ_max.

**Nahfeld und Volumenintegrale (Plan v2.25).** `V21NearField` und `V21c.OpticalTheoremAndGoldSphere`:

| Fall | Referenz | Abweichung |
|---|---|---|
| Brennebene eines Gauß-Strahls in Luft (x, z, radial) | vorgegebenes Profil | 2e-10 (Apertur), radial 5e-9 |
| verlustbehaftete Säule (n = 1,5 + 0,1i), Laurent, κ_max = 4 | A_Ringe = 1 − R − T | 4e-14 |
| dasselbe, Li in r über G±, κ_max = 2 / 4 / 8 | 1 − R − T | 2,3e-3 / 5,8e-4 / 1,5e-4 (direkt: 7,5e-3 / 3,5e-3 / 1,9e-3) |
| Wand der Säule, D_r innen/außen (Abstand 2e-6 µm) | stetig | FFF 1e-5, direkt 56 % |
| dielektrische Kugel (V21c, N = 20, κ_max = 4) | optisches Theorem = σ_sca | 1,1e-4 |
| Goldkugel, Laurent (N = 10, κ_max = 16) | optisches Theorem = σ_sca + σ_abs | Rundung |

Goldkugel (a = 50 nm, λ = 550 nm, n = 0,43 + 2,45i), Li in r, σ in 1e-3 µm² (Mie: 9,76 / 10,93 / 20,70):

| N, κ_max | σ_sca | σ_abs | σ_ext | Theorem |
|---|---|---|---|---|
| 20, 16 | 6,82 | 13,45 | 20,27 | 18,00 |
| 20, 32 | 8,24 | 12,43 | 20,67 | 20,02 |
| 20, 64 | 10,03 | 12,71 | 22,74 | 22,62 |
| 40, 32 | 7,60 | 12,49 | 20,09 | 19,48 |

Die Extinktion trifft Mie auf wenige Prozent, die Aufteilung konvergiert langsam (Treppenecken der Metallkugel);
der Plan verlangt hier nur einen qualitativen Vergleich.

**Strahl als Streufeld (V21e, Plan v2.27).** Zwei Ringschichten zwischen Filmen, κ_max = 2: R, T, A_Ringe je
Harmonischer gleich dem Totalfeld auf 1e-13 (Laurent und Li in r, x und radial), Nahfeld 1e-12, RK gegen RCWA 5e-9.

**`rings_profile` und CAD-Rotationskörper (Plan v2.28, v2.29).** Treppe = explizite `rings`-Folge bitgleich (Kegel,
Kugel); Treppenfehler eines Kegels zweiter Ordnung (Verhältnis 3,6 bei K = 3, 6, 12). Zylinder + Kegelstumpf +
Kugelkappe aus STEP = Stapel aus `rings`/`rings_profile` auf 2e-13.

**Winkelspektrum und Fernfeld (Plan v2.30).** Σ Gewicht·(s + p) = Leistung je Harmonischer auf 6e-17; φ-Mittel des
Fernfeldes = Summe der Winkelspektren auf Rundung. Kleine Kugel (n = 1,5, k₀a = 1, 16 Scheiben) gegen Mie
(|S1|², |S2|² in den Ebenen φ_k = 0 und 90°), Abweichung bezogen auf das Maximum:

| κ_max | 4 | 8 |
|---|---|---|
| max \|dσ/dΩ − Mie\| / max Mie | 2,2e-2 | 4,7e-3 |

**V21f (Symmetrie und Weg, Plan v2.31).** Lochblende r = 0,25 µm in 0,1 µm Metall auf Glas, Strahl 633 nm, κ_max = 2:

| Fall | Größe | Wert |
|---|---|---|
| Silber, m = −1 / +1 | Leistung | gleich auf 1e-12 |
| Silber, δ = 0,05 → 0,1, dk_factor 1 | ΔR, ΔT | 9,4e-5, 4,6e-5 (Diskretisierungsfehler 1e-4, 2e-5) |
| Silber, δ = −0,05 (falsche Seite) | ΔT | 0,09 |
| n = 0,5 + 3i / 0,3 + 3,5i / 0,15 + 3,9i | \|ΔR\| + \|ΔT\| reell gegen deformiert | 2,3e-8 / 2,7e-6 / 1,2e-4 |

Werte seit Plan v2.33 (gestufte Panels um das Plasmon, V21h); vorher 1,6e-5 und 1,4e-4 für δ sowie 6,5e-8 / 2,3e-4 /
4,6e-3 reell gegen deformiert.

**V21h (Konvergenzstudie, Plan v2.32, v2.33).** `scripts/v21h_convergence.py` (nur Standardbibliothek, ruft `fgsolve`,
etwa 20 min auf einem Kern) schreibt `docs/validation_data/v21h_convergence.csv`. Dielektrische Kugel (n = 1,5,
k₀a = 2) in Luft, ebene Welle unter 40° (Harmonische automatisch bis m_cut = 8), `rings_profile`, σ_sca gegen Mie
(0,572454 µm²); Standard κ_max = 4, `dk_factor` 1, 32 Scheiben:

| Parameter | Werte | σ_sca/σ_Mie − 1 |
|---|---|---|
| κ_max | 2 / 3 / 4 / 6 / 8 | −3,4e-2 / −1,5e-2 / −1,2e-2 / −7,1e-3 / −4,8e-3 |
| `dk_factor` | 0,5 / 1 / 2 | −1,21e-2 (gleich auf 5e-8) |
| Scheiben | 8 / 16 / 32 / 64 | −5,1e-3 / −1,0e-2 / −1,2e-2 / −1,3e-2 |
| m_cut | 1 / 2 / 3 / ≥ 4 | −5,7e-2 / −1,2e-2 / −1,2e-2 / −1,21e-2 (ab 4 gleich auf 1e-9) |

κ_max bestimmt den Fehler; die Stützstellendichte ist bei dielektrischen Körpern unkritisch, die Harmonischen sind
weit unter dem automatischen m_cut gesättigt. Der Treppenfehler hebt bei κ_max = 4 einen Teil des Abschneidefehlers
auf (daher wächst die Abweichung mit der Scheibenzahl auf den Grenzwert dieses κ_max zu).

Lochblende r = 0,25 µm in 0,1 µm Silber (n = 0,135 + 3,99i) auf Glas, Strahl 633 nm, κ_max = 2 (Referenz jeweils der
feinste Lauf):

| Fall | `dk_factor` | R | T |
|---|---|---|---|
| Luft darüber | 0,5 / 1 / 2 / 3 | 0,351471 / 0,350550 / 0,3504533 / 0,3504524 | 0,493471 / 0,494094 / 0,4940715 / 0,4940711 |
| Luft darüber, vor v2.33 | 1 / 2 / 3 / 4 / 5 / 6 | 0,330 / 0,328 / 0,349 / 0,361 / 0,352 / 0,350 | 0,416 / 0,453 / 0,537 / 0,504 / 0,477 / 0,495 |
| Glas darüber | 1 / 2 / 3 / 4 | 0,474004 / 0,473997 / 0,4740000 / 0,4740000 | 0,418350 / 0,418152 / 0,4181481 / 0,4181480 |

Mit Luft darüber liegt das Plasmon der Grenzfläche Silber/Luft (k ≈ 1,033 + 0,0023i) zwischen den
Verzweigungspunkten 1 und 1,5; es strahlt durch den Film ins Glas ab (leckende Mode) und erzeugt auf dem reellen Weg
eine Resonanz der Breite ≈ 5e-3. Gleichmäßige Panels lösen sie erst bei sehr dichten Stützstellen auf (vor v2.33:
`dk_factor` 6, T noch ±0,02). Seit Plan v2.33 setzt das Gitter gestufte Panels um Re k_p (Breiten Im k_p · 2^j). Damit
ist `dk_factor` 1 auf 2e-5 in T konvergiert, bei 224 statt 144 Stützstellen (Test
`V21h.LeakyPlasmonBetweenTheBranchPoints`). Mit Glas auf beiden Seiten sind
beide Plasmonen gebunden und liegen auf dem deformierten Weg: `dk_factor` 3 ist auf 4e-8 konvergiert. In κ_max
konvergiert auch der symmetrische Fall langsam (R 0,532 / 0,474 / 0,487 bei 1,5 / 2 / 3), wie bei der Goldkugel.

**V21i (Gradienten, Plan v2.34).** Gauß-Strahl (633 nm, w₀ = 0,4 µm) auf zwei Ringe (Kern n = 1,7 + 0,02i,
r = 0,18 µm; Schale ε = 1,96, r = 0,32 µm; 0,15 µm hoch) über einem Film (n = 1,6, 0,08 µm) auf einem Substrat
n = 2, κ_max = 2. Zielfunktion R + ½T + 2T(m = 1), Adjungierte gegen zentrale Differenzen (Schritt
1e-5 · max(0,1, |p|)), Kriterium wie M19 1e-6 · max(1, |fd|):

| Fall | Parameter | max \|adjungiert − fd\| |
|---|---|---|
| Li in r | Dicken (Film, Ringe), beide Radien, Kern n (Re, Im), Schale ε, Film n | 2,9e-7 |
| Laurent | beide Radien, Schale ε | 7,0e-8 |
| Kern n = 0,5 + 3i, reeller Weg | Ringdicke, innerer Radius, Kern n (Re, Im) | 2,7e-7 |

Ebene Welle (Stufe 2, Plan v2.35): dieselben Ringe unter φ = 30°, ψ = 20°, Zielfunktion
σ↑ + ½σ↓ + 2σ↓(m = 0) in µm². Beide Radien, Kern n (Re, Im), Schale ε: max |adjungiert − fd| = 1,2e-7 (Li in r)
bzw. 6,6e-8 (Laurent). Beide Dicken und Film n, auch als Außenmaterial der Ringe (Hintergrund über zentrale
Differenzen der Einfallsordnung): ≤ 2e-10.

Halbraummaterialien (Plan v2.48), Differenzen bei festem Spektralgitter mit Richardson-Extrapolation (die
Verzweigungspunkte folgen den Halbräumen): Strahl, Substrat n = 2 + 0,05i (Re, Im) und Superstrat „air" (zugleich
Außenmaterial der Ringe) ≤ 1,8e-8; ebene Welle, Substrat ≤ 2e-9; Klassenpfad P = 2 (V22g-Fall, „air" zugleich
Superstrat, Hintergrund der Ellipse und Außenmaterial der Ringe) ≤ 9,4e-9. Superstrat der ebenen Welle (Plan v2.51,
über k∥ der Quellen, den Hintergrund und die Normierung; Richardson mit h = 4e-6/2e-6, darunter rauschen die
Differenzen): 1,1e-8, Klassenpfad 5e-10.

Strahl als Streufeld (Plan v2.49, Gradient über die Totalfeldkette, Differenzen der Streufeldrechnung): Li in r, Film
als Außenmaterial der Ringe (äußerer Radius, Im n Kern, Film n, Filmdicke) ≤ 3,2e-10; Klassenpfad P = 1 (Fünfeck,
Film als Hintergrund: Eckpunkt, Film n, Dicke) und P = 2 (Ellipse über Ringen: Halbachse, Ringradius) ≤ 1,7e-8.

`rings_profile`-Parameter (Plan v2.50, über die Treppe): Kugel innen, Tabelle außen, drei Scheiben unter einem Film;
Strahl (Dicke, Kugelradius und -mitte, Tabellenwert, Im n Kern) und ebene Welle (Dicke, Tabellenwert) ≤ 9,7e-10;
Klassenpfad P = 2 (Kegel in zwei Scheiben unter der Ellipse: Kegelradius, Dicke) ≤ 6e-9. Liegt eine Scheibenmitte auf
einem Tabellenknoten, ist die Ableitung nach der Dicke dort nur einseitig definiert (zentrale Differenzen 3,5e-7).

Die Rückwärtsregeln selbst (`CylAdjoint.*`, Fluss, Propagationsmatrix, Galerkin-Matrizen und Faktorisierung, Laurent
und Li, m = 0, 1, −2, 3, komplexe ε auf dem deformierten Weg) stimmen auf 1e-6 relativ mit Differenzenquotienten überein.

**V21g (Lochblende im Metallschirm, Plan v2.36).** Qualitativer Vergleich mit Popov/Bonod, phys. stat. sol. (b) 244,
3463 (2007). Ebene Welle 500 nm, senkrecht, E längs x; Silber/Aluminium nach Rakić (ε_Ag = −7,632 + 0,731i,
ε_Al = −34,23 + 8,98i; Annahme, die Arbeit nennt keine Werte), Schirm 200 nm (Abb. 4, 6) bzw. 50 nm (Abb. 7, Annahme).

| Abbildung | Aussage der Arbeit | hier (Test) |
|---|---|---|
| Abb. 4, R = 125 nm, \|E\|² 15 nm unter dem Schirm | Spitzen am Rand in Polarisationsrichtung, quer rascher Abfall | Spitze 1,89 bei \|x\| = R, Mitte 0,90; bei 250 nm längs x 0,37, längs z 0,014 (`Fig4EdgePeaksAlongThePolarization`, κ_max = 8) |
| Abb. 6, R = 50 nm, E_y längs x | \|E_y\| ∝ \|H₁(k_p x)\| (zylindrisches Plasmon) | \|E_y\|/\|H₁\| 0,89 … 1,00 über 1–5 µm (fgsolve, dk_factor 100); Test 0,5–1 µm ≤ 1,05 (`Fig6CylindricalPlasmonAlongTheSurface`) |
| Abb. 7, R = 10 nm, Fernfeld | in der E-Ebene flach bei idealem Leiter, gerichteter mit sinkender Leitfähigkeit; quer (H-Ebene) für alle Metalle gleich | dP/dΩ normiert auf θ = 0, bei 69°: E-Ebene Ag 0,38, Al 0,74; H-Ebene beide 0,13 = cos²θ (auf 4 %); k_max = 80 k₀ wie 160 und 320 k₀ (`Fig7FarFieldOfASmallHole`). 33 ε_Al: flacher als Al (69°: 0,82 … 0,99), aber in k_max nicht konvergiert (P(0) 2,9e-10 / 4,2e-14 / 5,6e-15 bei 80 / 160 / 320 k₀; Signal σ↓ ≈ 2e-14 µm²); 1000 ε_Al flach (0,98 … 1,02 bei 69° für alle k_max), P(0) ebenso nicht konvergiert (3,6e-10 / 2,6e-11 / 1,7e-13); beide daher nicht im Test |

Befund: Die Amplitude der 50-nm-Lochblende konvergiert erst ab κ_max ≈ 24 (σ_sca 0,0051 / 0,0012 / 0,0015 / 0,0021 µm²
für κ_max 4 / 6 / 8 / 24); die Form des Plasmons ist schon bei κ_max = 4 richtig. κ_max muss bei kleinen Löchern in
Metall in Einheiten von 1/R gewählt werden (k_max R ≈ 15). Auf feinen Gittern liegt der Knoten k ≈ 0 beim
Hintergrund-Exponenten (Abstand 2,1e-8 < 1e-8 max|λ|): früher RK-Rückfall (Abb. 6: über 90 min), jetzt die dividierte
Differenz; gegen ein Gitter ohne Resonanz σ auf 5e-14, Nahfeld auf 7,5e-12 (`V21g.BackgroundExponentNextToTheNodeKZero`).

**V21d (gegen kartesisch mit PML, Plan v2.38).** Normaleinfall, Li in r; kartesisch `contour2d` mit PML 0,3 µm in x und
z. Verglichen wird der gesamte Streuquerschnitt: M18 oben + unten + seitlich (nur die Summe ist vom Fenster
unabhängig), M20 σ_sca bzw. beim Film die Extinktion aus dem optischen Theorem (die geführten Moden fließen seitlich ab
und fehlen in oben/unten).

| Fall | M20 | M18 | Abweichung |
|---|---|---|---|
| Säule n = 1,45, r = 0,4 µm, 0,3 µm auf Glas, λ = 1 µm | 0,38965 µm² (κ_max = 8; 16: 0,38972) | Zelle 2 µm, N = 4/6/8: 0,38458 / 0,38723 / 0,38903; Zelle 2,4 µm, N = 10: 0,39006 | 6,2e-4 (N = 8) |
| Kontaktloch wie V20g (r = 0,25 µm, Film n = 2,1, 0,3 µm, 633 nm) | Extinktion 0,554838 (κ_max = 4), 0,554369 (κ_max = 8); oben + unten 0,4339 | Zelle 2 µm, N = 8/10/12: 0,5700 / 0,5686 / 0,5640 (Absorption bis 6e-3) | Zielsystem (`FG_SERVER_TESTS=1`, N = 14/18/22) |

Befund (Planfehler): Der Film führt TE₀, TM₀ und TE₁ (n_eff 1,9655 / 1,9026 / 1,5666). Mit dem Standardweg
(δ = 0,05 (k_b − k_a)) liegt TE₁ nur 0,010 unter dem Sinusweg, die Knoten haben bei R_s = 0,25 µm den Abstand ≈ 0,08:
Extinktion 0,6146 / 0,5485 / 0,5170 bei dk_factor 0,5 / 1 / 2. Seit v2.38 bestimmt der Löser die Moden des
verlustfreien planaren Hintergrunds und legt gestufte Panels um sie: 0,554838 für dk_factor 0,5, 1 und 2 gleich, wie ein
tiefer Weg mit doppelter Knotenzahl (δ = 0,3: 0,554877). N_k steigt dabei von 128 auf 312.

## P-zählige Drehsymmetrie (M21, V22)

`tests/validation/test_v22_cylp.cpp`; Säule n = 1,45 (0,3 µm) auf Glas, λ = 1 µm, κ_max = 3 (V22c: siehe Zeile).

Konvergenz V22c (Quadrat, lokal): Klassenpfad L = 1 bei κ_max = 4/6/8: 0,30063 / 0,29951 / 0,29919; κ_max = 6, L = 2:
0,29852. M18 N = 6/8: 0,30314 / 0,30206 (Absorption in der Box bis −2e-3). Beide nähern sich von oben.

| Test | Fall | Ergebnis |
|---|---|---|
| V22a | Ringschicht (Li in r) über den Klassenpfad, P = 1, 2, 3, L = 1, Strahl w₀ = 0,25 µm | R, T je Harmonischer = M20 auf 3e-13, dunkle Harmonische < 1e-14 |
| V22a | zentrierter Kreis als `contour2d` (P = 3, L = 2, Laurent) gegen `rings` (Laurent), φ = 30° | 3e-13 µm² bei σ = 0,24 µm² (L = 5: 4e-14) |
| V22b | exzentrischer Faserkern (P = 1) gegen zentrierten Kern, HE₁₁ | 8,4e-5 (x₀ = 0,3, L = 3) … 1,6e-10 (x₀ = 0,6, L = 7) |
| V22c | Quadrat 0,6 µm (P = 4), n = 1,45, Normaleinfall, gegen M18 mit PML (Zelle 2 µm), beide Laurent | lokal grob: Klassenpfad κ_max = 4, L = 1: 0,30063 µm², M18 N = 6: 0,30314 (2,5e-3); Abnahme auf dem Zielsystem (`FG_SERVER_TESTS=1`: κ_max = 8, L = 2/3 gegen N = 10/14; auch Sechseck und Silber) |
| V22d | Quadrat: Blöcke l − l' ∉ 4ℤ (eine Klasse von P = 1) | 6e-18 gegen 1,4e-2 |
| V22d | Quadrat, P = 4, L = 1, Einfallsazimut Ψ = 0° / 90° | σ gleich auf 1e-15, Fernfeld gedreht auf 4e-15 |
| V22e | Quadrat, P = 4, L = 2, Strahl zirkular + / − | R, T und Harmonische m ↔ −m auf 1e-13 |
| V22Volume | verlustbehaftetes Quadrat (n = 1,45 + 0,1i), P = 4, L = 1, Laurent | Strahl: Ringabsorption = 1 − R − T auf 1,3e-13; ebene Welle: Theorem = σ_sca + σ_abs auf 4e-14 relativ, σ_abs bei Ψ = 90° bitgleich |
| V22Volume | verlustbehaftete Ringschicht über P = 3 gegen M20; Kreis (P = 3, L = 2) gegen Ringschicht, Laurent | Ringabsorption 5e-15; σ_abs 5e-14, Theorem 9e-14 µm² |
| V22g | Gradienten im Klassenpfad (Strahl, κ_max = 2, L = 1, Kriterium von M19): P = 1, unregelmäßiges `contour2d`-Fünfeck (Eckpunkte, Materialien, Dicke); P = 2, Ellipse (Halbachse, Winkel) über Ringen mit Li in r (Radius, Materialien); Halbachse einer Kreis-Ellipse (Randterm, entartete Paare m, −m; Plan v2.52) | ≤ 2,9e-9 bzw. 6,5e-8, Kreis 2,7e-8; symmetriebrechender Parameter abgelehnt |
| V22g (Streufeld) | Gradienten des Strahls als Streufeld im Klassenpfad (über die Totalfeldkette, Plan v2.49): P = 1 Fünfeck mit Film als Hintergrund (Eckpunkt, Film n, Dicke), P = 2 Ellipse über Ringen (Halbachse, Ringradius) | ≤ 1,7e-8 |
| V22g (`rings_profile`) | Kegel in zwei Scheiben unter der Ellipse (P = 2): Kegelradius, Dicke | ≤ 6e-9 |
| V22g (ebene Welle) | Gradienten im Klassenpfad für eine schräge ebene Welle (φ = 30°, ψ = 20°, Harmonische −1..1, Zielfunktion σ↑ + ½σ↓ + 2σ↓,m=1): P = 1 Fünfeck (Eckpunkt, Kernmaterial, Filmmaterial, Dicke), P = 2 Ellipse über Ringen (Halbachse, Ringradius, Im n Kern, Dicke) | ≤ 4,6e-10 |
| CAD_V22f | Sechseck (P = 6, L = 1) aus STEP gegen `polygon2d`; Kegelstumpf (P = 4, L = 0) als `contour2d`-Treppe gegen `rings_profile`; gestörtes Sechseck, Quadrat mit P = 6 / 2, Dreieck (Achse im Schwerpunkt) | R, T 1,3e-14 bzw. 2,5e-14; `E_ROT_SYMMETRY` mit Ordnung 1 bzw. 4, `rotation_order` |
| V22NearField | Nahfeld im Klassenpfad: Ringschicht (Li in r, Strahl) gegen M20; Kreis gegen Ringschicht (ebene Welle); Materialkarte des Quadrats | 1,5e-12 bzw. 3e-10 relativ (Schnitte bis R_s) |

## CAD-Import für 1D-Gitter (M23, V24)

`tests/cad/test_cad.cpp` (`CAD_V24`, Label `validation`, nur mit Nurbforge), ohne Nurbforge
`tests/unit/test_line_section.cpp` und `tests/unit/test_sectioned1d.cpp`, Python `tests/python/test_cad.py`.

| Fall | Referenz | Abweichung |
|---|---|---|
| Linienschnitt Kreis (rationale Quadriken), bis 1e-9 an die Tangente | geschlossene Form | 2e-15 + 4e-16 r/h |
| Normalenfeld eq:cad1d_N | geschlossene Form der Koeffizienten | ≤ 1e-14 |
| V24a Rechtecklinie (senkrechte Wände) | `lamellar1d` (Li); Analytik-Provider: Q | 1e-12 / Q 1e-13 |
| V24b Trapez, 15°-Wände, Laurent | `profile1d` `points`: RCWA / RK; Q an 5 Höhen | 1e-12 / 1e-10; Q 1e-13 |
| V24b FFF_sym (anderes Normalenfeld im Inneren) | `profile1d` FFF_sym | 5,8e-5 (N = 5) → 1,3e-5 (N = 20) |
| V24c gerundetes Trapez der Meldung (CD 100 nm, h 100 nm, r 10 nm) | `profile1d` `Function` mit analytischem g(x), RK | 2,6e-13 |
| V24c Schnittpunkte / Normalen auf den Bögen | Analytik | ≤ 1e-13 / ≤ 1e-12 |
| V24d Überhang mit Deckel (zweites Material) und Hohlraum | RCWA (100 Scheiben) gegen RK; R + T; Reziprozität | 2,4e-10; 1e-14; δ = 0°: 8,5e-7 → 3,4e-9, δ = 90°: 2,2e-5 → 3,5e-7 (N = 8 → 24) |
| V24e z_cut, Extrusionslänge 0,1 … 1000 µm, Verschiebung über den Zellrand | Ausgangsfall | ≤ 1e-13 / 1e-12 |
| V24f Restschicht + Linie als ein Körper | zwei Körper; Selbstüberlappung mit `periodic_union` = homogene Schicht | 1e-12 / 1e-13 |
| V24g Befunde | `NOT_Z_INVARIANT`, `OUT_OF_RANGE`, `BODY_OVERLAP`, `priority` = `lamellar1d` | genau / 1e-12 |
| V24h Gradienten (Oberbreite über `StackParameter`, n) | zentrale Differenzen, RCWA und RK | ≤ 1,7e-8 / ≤ 2,1e-9 relativ |
| V24i Python, Beispiel `cad_line_grating_1d` | Golden File (neu) | bitgleich |

## NURBS-Profil in 1D-Gittern (M24, V25)

`profile1d` mit `"shape": "nurbs"` (Plan v2.37). `tests/unit/test_nurbs_profile.cpp` (V25a),
`tests/validation/test_v25_nurbs.cpp` (V25b–g, Label `validation`), `CAD_V25` im CI-Job `cad`, Python
`tests/python/test_profile1d_nurbs.py`. Periode 0,3 µm, 633 nm konisch (φ = 20°, ψ = 30°, δ = 35°) wie V24.

| Fall | Referenz | Abweichung |
|---|---|---|
| V25a Bézierstücke (A5.6) | NURBS nach de Boor | < 1e-14 |
| V25a y-Extrema kubisch / Kreisbögen (rational) | geschlossene Form | 1e-14 |
| V25a Linienschnitt (Höcker, Pilz mit Überhang) | `lineSection` des geschlossenen Körpers | Bruchstellen und Normalen < 1e-13 |
| V25b Grad 1 = Trapez (`points`), isotrop / anisotrop, Laurent | `profile1d` `points` | RCWA 1,4e-15 / 5,6e-16, RK 7,9e-14 / 3,1e-13 |
| V25b Grad 1 = Dreieck | `profile1d` `triangle` | RCWA 1,1e-15, RK 5,8e-14 |
| V25c gerundetes Trapez (exakte Bögen), Bandgrenzen an den Bogenenden | `profile1d` `Function` (RK, Laurent) | 2,1e-13 |
| V25c dieselbe Linie | `cad_line_grating_1d` (STEP, drei Bandschichten; RK, tol 1e-10) | 3,2e-10 (RK-Toleranz; CI-Job `cad`) |
| V25d Pilz (glatte Quadriken, Überhang) | R + T = 1 | RK 1,7e-14, RCWA ≤ 8e-15 |
| V25d RCWA → RK | 20 / 80 / 320 Scheiben | 2,5e-3 / 1,4e-5 / 2,3e-6 |
| V25d Kurvenanfang verschoben | Ausgangsfall | RCWA 4,4e-16, RK 1,1e-15 |
| V25d Reziprozität (δ = 90°) | N = 8 → 24 | 1,9e-6 → 2,0e-7 |
| V25d Abbruch, R₀ gegen FFF_sym N = 56 | FFF_sym N = 8; Laurent N = 8 / 24 | 2,0e-5; 6,4e-4 / 2,3e-4 |
| V25e Naht glatt / Knick, Anfang verschoben | Ausgangsfall | ≤ 1,3e-14 |
| V25f Tschebyschow (16 Knoten je Band) | exakte Schnitte | 4,8e-15 |
| V25f ASR (η = 0,5), N = 16 | ohne ASR; Energie | 1,6e-6; 2,0e-7 (~N⁻⁴, wie M17) |
| V25g Fehler mit Pointer und Kennung, Gradient nach `depth_um`/Punkt abgelehnt, Sweep über einen Punkt | Neuberechnung | < 1e-14 |
| V25h Beispiel `profile1d_nurbs` | Golden File (neu) | bitgleich |
