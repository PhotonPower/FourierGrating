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

Der Peak liegt stets oberhalb des SPP-Pols der Ag/Glas-Grenzfläche (1,3694 µm), wandert aber mit N nach Rot
und seine Höhe schwankt; die Buchangabe liegt im Streubereich, ist aber nicht bestätigt. Laurent schwankt
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
