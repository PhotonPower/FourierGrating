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
Kreuzungen und Quadraturlinien), es bewegt sich das Normalenfeld. Offen: Rückwärtsfunktionen über das GPU-Backend,
aperiodische Strukturen (PML), Verschiebungen x0/z0 von `profile2d`. Rasterfeld von Algorithmus NF (v2.12): unterhalb
einer Rasterzelle springt die Vorwärtsrechnung um ~1e-7 (Kreismittelpunkt in 2-nm-Schritten), der volle
Differenzenquotient hängt dort von h ab (−0,82 / −0,086 / −0,0077 / −0,0012 für h = 1e-6 … 1e-3); eingefroren erzeugt
das stehenbleibende Feld eine Scheinableitung. Formparameter mit FFF bewegen das Feld deshalb mit
(Differenzenquotient der NF-Koeffizienten bzw. der FFF_B-Gitter über eine Rasterzelle, zwei NF-Rechnungen je
Parameter); der Indikatoranteil bleibt exakt. Anisotrope FFF_B-Querschnitte: Materialparameter als
Richtungsableitung der Gitter (zwei Auswertungen je Parameter).
