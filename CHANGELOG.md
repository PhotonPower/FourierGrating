# Changelog

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
