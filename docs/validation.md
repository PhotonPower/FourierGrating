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
| V7 | Co-Sinus anisotrop: d = 0,6, a = 0,5 µm, λ = 0,6328 µm, 30°, TM-Einfall, ε_ii = −8,19+16,38i, ε_xz = −ε_zx = −0,495−0,106i | TM₀ 0,0567; TM₋₁ 0,1435; TE₀ 1,3e-7; TE₋₁ 2,37e-6 (intern, N = 20, 100 Scheiben) | 3 % |
| V8 | anisotroper Code, isotroper Tensor | = isotroper Code | 1e-10 |
| V9 | 2D-Code, z-invariante Struktur | = 1D konisch | 1e-10 |
| V10 | Ag-Lochgitter: Film 0,2 µm auf Glas, Quadratlöcher 0,25 µm, Periode 0,9 µm, ν_Ag = 0,1+8,94i, Normaleinfall | T(1,388 µm) ≈ 7,5 % (Buch); SPP-Pol eben: λ = 1,369406 µm | 10 % |
| V11 | Al-Sinus 40°, ν = 1,3+7,6i: Li-Lamellentreppe vs. FFF (Buch Abb. VI.10) | Treppe langsam, FFF schnell | qualitativ |
| V12 | Nahfeld: Stetigkeit tangentialer Komponenten; Feld im Superstrat = Rayleigh-Summe | Sprung < 1e-6 | 1e-6 |

Zusätzliche Kontrollen aus der Kapitelprüfung (Al-Sinus V4, TM, FFF-Treppe):
60 / 150 / 300 Scheiben → 0,853 / 0,865 / 0,869. Treppe mit Li-Lamellenformel bei N = 30:
M = 1 → 0,565, M = 20 → 0,804, M = 200 → 0,774 (keine monotone Konvergenz).
