# Buchfehler mit Relevanz für die Implementierung

Nevière/Popov, *Light Propagation in Periodic Media* (2003). Alle Punkte wurden in der Prüfung
der Kapitel numerisch bestätigt; der Plan verwendet jeweils die korrigierte Form.

| Stelle | Buch | korrekt |
|---|---|---|
| (II.19)/(II.20) | β₀⁽¹⁾ im einfallenden Term | β₀⁽²⁾ |
| TM-Matrix S. 46 | Indizes (m,n) vertauscht | M₁₂,ₙₘ = (k²)ₙ₋ₘ usw. |
| (III.24)/(III.26) | C und C′ vertauscht | Grenzfläche über Ψ_ref⁻¹ Ψ_X |
| S. 56 | A′⁽¹⁾ = A⁽¹⁾ | nur bei ν₁ = ν_M |
| (III.29) | T₁₂, T₂₁ vertauscht | Sternprodukt |
| (III.45) | S₂₁ statt S₂₂; T₁₂⁽q⁾ statt T₁₂⁽q+1⁾ | Sternprodukt (Plan eq:star) |
| (IV.39)/(IV.40) | δ_{n,p} fehlt | Einheitsstartvektoren |
| konisch | β ohne γ₀ | β = sqrt(k² − α² − γ₀²) |
| S. 102 (Text) | Spaltenreihenfolge von C_ε | (V.11) ist richtig |
| S. 108 | C_ε einmal berechnen | hängt über ε(y) von y ab |
| S. 111 | Faktor 1/(ω²μ₀) in Zwischenschritt | Endergebnis richtig |
| (VI.13) | M̂ᵢⱼ = Mᵢⱼ + L | M̂ᵢⱼ = Mᵢⱼ + δᵢⱼ L |
| (VII.72) | Vorfaktor √(μ₀/ε₀)/(ik_j) | 1/(iωε_j) |
| (VIII.2) | Vorzeichen von C₂₃ | −(N·ε)_z N_y / ξ |
| (VIII.3) | M₃₂ ohne −iα²/(ωμ₀) | mit |
| (IX.1) | β₀ = k₀ cos Ψ | k₀ cos Φ |
| (IX.14) | Int[(r−1)/(N+1)] | Int[(r−1)/(2N+1)] |
| (IX.16)–(IX.19) | Indizes m(s), α_{m(r)} | alle zu r, α_{n(r)} |
| (IX.20) | Transmission ohne k₂ → k₁ | mit |
| (IX.22) | e⁴_{g₂} a₂² | a₂⁴ |
| (IX.28) | d[E_y]/dy | d[E_z]/dy |
| S. 197 | ρ² = −129,0973, Faktor 0,00534 | −136,88 µm⁻², Faktor 0,096 |
| Anh. IV | 2D-Blockgröße (4N+2)² | 2(2N+1)² |
| Anh. IV (A.IV.7) | ε_j in q_h | ε_M |

Zusätzlich zu beachten: Die Formeln in Anh. IV dividieren durch k² − γ_m², das im 2D-Fall
verschwinden kann; der Plan verwendet stattdessen die Basis der tangentialen E-Komponenten
(nur β im Nenner).
