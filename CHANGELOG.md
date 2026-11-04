# Changelog

## [0.4.0] – M5 (anisotrop, 1D)

### Hinzugefügt
- Anisotrope FFF: Variante A (Buch Kap. V/VIII, C_ε) und Variante B (Projektoren, Eigenherleitung);
  Fourierkoeffizienten stückweise glatter Funktionen per Gauß-Quadratur.
- Anisotrope Modenbasis homogener Medien (4×4 je Ordnung, Klassifikation); anisotrope Substrate,
  homogene Schichten, Lamellen (N = x konstant) und Profile; Transmission über den Modenfluss.
- Materialeingabe über Hauptachsen + Euler-Winkel; s/p-Anteile der Effizienzen (CSV-Spalten eff_s, eff_p).
- Validierung V7 (Kobalt), V8 (isotroper Tensor), einachsiges Substrat (analytisch), biaxiale Schicht.

### Behoben
- Eigen konjugiert bei `dot()`/`cross()` komplexer Vektoren: ξ und C_ε waren falsch, s/p-Amplituden konjugiert.

## [0.3.0] – M3, M4

### Hinzugefügt
- Profilschichten y = g(x): Sinus (exakte Schnittpunkte), stückweise linear (Dreieck, Echelette, Punkte;
  exakte Koeffizienten), Funktion (Nullstellensuche); Normalenfeld (NN und vorzeichenbehaftet).
- RK mit M(y) (exakte Koeffizienten an jeder Stufenhöhe), RCWA-Treppe mit glattem Normalenfeld.
- FFF-Varianten `FFF_sym` und `FFF_rig` (Plan eq:Qsym, eq:Qrig); Standard je nach Profilform;
  Warnung für die Buchvariante bei schwach absorbierenden Metallen.
- TE/TM-Entkopplung im in-plane-Fall (RK ≈ 6×, RCWA ≈ 2,5× schneller).
- Validierung V4 (Al-Sinus gegen C-Methode), V5, V6, konische Profile, V13 (Gold-Sinus gegen C-Methode),
  V14 (Gold-Dreieck, offen); Paarungs- und Profil-Unit-Tests.

### Behoben
- RK-Scheibenteilung: Spektralradius aus der abgeschnittenen Matrix statt Materialschätzung,
  plus Halbierung bei zu großem Wachstum (vorher Kontamination bei Gold).
- Modenklassifikation in metallischen Scheiben (Sortieren und Halbieren statt Schwellenwert).

### Befund
- Buch (IV.32) liefert für Gold unphysikalische Effizienzen (Σe > 1); FFF erhält die Energie nicht exakt.

## [0.2.0] – M1 abgeschlossen, M2

### Hinzugefügt
- Q-Matrizen: Laurent, Li-Lamellen, FFF isotrop (Plan eq:Qiso); Abtastung stetiger Funktionen.
- RCWA: reduziertes und allgemeines Eigenwertproblem mit Klassifikation, stabile Scheiben-S-Matrix.
- RK: Transfermatrix für konstantes M (Schrittverdopplung) und adaptiv für M(y); Scheiben-S aus T.
- Stapellöser für beide Verfahren: Substrat-/Referenzgrenzfläche, homogene und Lamellenschichten,
  Effizienzen aus dem Poynting-Fluss, s/p-Amplituden, Energiebilanz.
- JSON-Konfiguration (homogeneous, lamellar1d, Wellenlängen-/Winkel-Scans), `fgsolve` mit CSV-Ausgabe.
- Validierung V1 (Fresnel, konisch, verlustbehaftet), V2 (Airy, RK und RCWA), V3 (Energie,
  RK = RCWA, Reziprozität); 41 Tests.

### Geändert
- Plan v1.1: volles 4n_f-System auch in-plane (2n_f-Systeme als spätere Optimierung);
  Reziprozität im FMM nur bis auf Abschneidung.

## [0.1.0] – Grundgerüst

### Hinzugefügt
- Implementierungsplan (LaTeX/PDF) mit allen Formeln, Herleitungen, Tests und Meilensteinen.
- CMake-Projekt (Presets CPU/CUDA), Abhängigkeiten (Eigen3, nlohmann_json, GoogleTest), CI (CPU).
- Implementiert und getestet: `IndexMap` (korrigiert Buch IX.14), exakte Fourierkoeffizienten
  (Intervall, Rechteck, Polygon), Toeplitz-Matrizen (1D, 2D), S-Matrix (T→S, Sternprodukt,
  interne Amplituden), Butcher-Tableaus (Gauß 2/3, Radau IIA 3) und RK-Schrittpropagator,
  Propagationsmatrix `M` (μ = 1), isotrope Modenbasis, einfallende Welle (korrigiert Buch IX.1).
- Schnittstellen (noch `NotImplemented`): Q-Matrizen (FFF), anisotrope Moden, Löser RK/RCWA,
  Nahfeld, Konfiguration, CUDA-Backend.
