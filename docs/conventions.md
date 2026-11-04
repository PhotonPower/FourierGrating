# Konventionen (Kurzfassung des Plans)

* Zeitabhängigkeit `exp(-iωt)`.
* Periodizität in x (dx) und z (dz, nur 2D); Schichtung in y. Substrat y < 0 (Medium 1),
  Superstrat oben (Medium M), Licht fällt von oben ein.
* Normierte Einheiten im gesamten Kern: Längen · k0, `H̃ = Z0 H`, also ω = μ0 = ε0 = 1.
  Maxwell: `rot E = i μ_r H̃`, `rot H̃ = −i ε_r E`.
* α_n = α0 + n Kx, γ_m = γ0 + m Kz mit Kx = 2π/(k0 dx). Einfall: k = (α0, −β0, γ0),
  α0 = ν sinΦ cosΨ, β0 = ν cosΦ, γ0 = ν sinΦ sinΨ.
* β = sqrt(ε − α² − γ²), Im β ≥ 0, Re β ≥ 0 falls Im β = 0; |β| < 1e-10 → β + 1e-10 i.
* Zustandsvektor `F = ([Ex], [Ez], [Hx], [Hz])`, dF/dy = M F.
* Einzelindex r = (m + Nz) + (2Nz + 1)(n + Nx), 0-basiert.
* Amplituden in homogenen Medien: tangentiale E-Komponenten, `V = (a_x, a_z | b_x, b_z)`,
  a = abwärts, b = aufwärts, Phasenbezug lokal an der jeweiligen Grenzfläche.
* S-Matrix: `(b_oben; a_unten) = S (b_unten; a_oben)`.
* Moden-Klassifikation: Feld ∝ exp(λy); aufwärts ⇔ Re λ < −τ oder (|Re λ| ≤ τ und S_y > 0).
