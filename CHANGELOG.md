# Changelog

## [Unreleased]

### Hinzugefügt
- M20/M21, Gradienten nach den Halbraummaterialien (Plan v2.48): Substrat (Re, Im) und Superstrat (Re) für Strahlen,
  Substrat für ebene Wellen, über die Basen der Halbräume (`fg::cyl::homogeneousModesBar`), den einfallenden Strahl
  und P_inc; Spektralgitter fest. C++-API `CylindricalSpec::fixedK`/`fixedW`/`fixedEdges` für Prüfungen bei festem
  Gitter. Tests `V21i.HalfSpaceMaterials`, `V22g.HalfSpaceMaterialsOfTheClassPath`, `CylAdjoint.HomogeneousModesBar`.

### Behoben
- M21: Ringschichten einer Klasse werden je Harmonischer in Moden zerlegt (`CylBasis::modes`). Die Harmonischen m
  und −m haben dieselben Eigenwerte; die Zerlegung der ganzen Klassenmatrix mischte die Blöcke, kostete etwa 1e-9 im
  Ergebnis und verfälschte den Gradienten nach dem Außenmaterial der Ringe (V22g-Fall: 2 %).
- M21: Formparameter, die ein Randstück auf einem Kreis um die Achse bewegen (Kreis, Kreisring), lieferten einen
  falschen Gradienten (0); sie werden jetzt mit einer Fehlermeldung abgelehnt.

### Hinzugefügt
- M21, Gradienten für ebene Wellen im Klassenpfad (Plan v2.47): die Quellenkette von V21i Stufe 2 läuft je Klasse
  (Quelle der Klasse aus `CylBasis::planeWaveSource`), Querschnittsschichten über die neue Rückwärtsregel
  `fg::cyl::coupledSourceBar` (Materialien und Form über `coupledQShape`), Ringschichten einer Klasse über
  `planeWaveSourceBar` je Harmonischer; Hintergrundgrößen (Dicken, Hintergrundmaterialien) wie Stufe 2b. Tests
  `V22g.PlaneWaveThroughTheClassPath` (≤ 4,6e-10), `CylCoupled.SourceReverseRuleAgainstDifferences`.
- M21, Strahl als Streufeld im Klassenpfad (Plan v2.46): `"scattered_field": true` mit `"symmetry": {"rotation": P}`;
  der Hintergrundstapel (Ringe → Außenmaterial, `polygon2d`/`contour2d` → `background`) ist diagonal in den
  Harmonischen und wird je Klasse auf dem Gitter gelöst, Quelle σ = (M − M_bg) F_bg wie V21e. Test
  `V22Scattered.BeamScatteredFieldEqualsTotalField` (R, T, A_Ringe je Harmonischer 6e-14, Nahfeld 4,6e-13 relativ).
  M20 unverändert.
- M21, Beschleunigung des Klassenpfads (Plan v2.45): Halbräume und homogene Schichten der Hankel-Basis gruppenweise
  (`groupedLayerSMatrix`, `groupedFromTransfer`, `groupedStar`, Marke `SMatrix::groups`, `star` mit einem gruppierten
  Faktor in O(n²) außer den zwei LU), Ψ⁻¹W der strukturierten Schichten gruppenweise, Quellen der ebenen Welle über das
  Schur-Komplement (2n_f statt 4n_f). Säule P = 3, L = 2: 29 s → 16 s; Quadrat zwischen homogenen Schichten 7,2 s → 2,6 s;
  volle lokale Testsuite 49 → 36 min. Test `CylGrouped.EqualTheDenseOperations`.
  Kartesische Ergebnisse bitgleich.
- M21 V22g, Gradienten im Klassenpfad (Plan v2.44): Gauß-Strahlen mit `"symmetry": {"rotation": P}` nach Dicken,
  Ringradien, Materialien und jeder Zahl einer `polygon2d`/`contour2d`-Geometrie, die die P-zählige Symmetrie erhält.
  Formableitung der Winkelkoeffizienten über die Bogenenden (`fg::cyl::angularCoefficientsShape`), Rückwärtsregel
  von `coupledQ` (`fg::cyl::coupledQBar`, `coupledQShape`); die M20-Kette läuft je Harmonischer unverändert. Tests
  `V22g` (Adjungierte gegen Differenzenquotienten ≤ 6,5e-8), `CylAngular.ShapeDerivativeAgainstDifferences`,
  `CylCoupled.ReverseRuleAgainstDifferences`.
- M21 V22f, CAD mit Drehsymmetrie (Plan v2.43; die LLM-Referenz nennt `E_ROT_SYMMETRY`): Schicht `cad` in einem Zylinderstapel mit `"symmetry": {"rotation": P}`
  (ohne `"axisymmetric"`); Achse ohne `cell_origin_um` im Schwerpunkt des untersten Schnitts, Bänder `contour2d`
  bzw. `contour2d_staircase`, Befund `E_ROT_SYMMETRY` mit der gefundenen Ordnung, Bericht `"dimension": "rotation"`
  mit `rotation_order`. Tests `CAD_V22f` (Sechseck aus STEP = `polygon2d` auf 1,3e-14, Kegelstumpf-Treppe =
  `rings_profile` auf 2,5e-14, Befunde).
- M21, Volumenintegrale und Nahfeld im Klassenpfad (Plan v2.42): σ_abs, optisches Theorem und die Ringabsorption
  eines Strahls mit `"symmetry": {"rotation": P}` (Querschnittsschichten über die Winkelkoeffizienten von Im ε und
  ε − ε_∞ auf der r-Regel von `coupledQ`, `fg::cyl::sectionVolume`; Ringschichten je Harmonischer wie M20), Nahfeld je
  Klasse in die Harmonischen zerlegt. Laurent erhält die Energie diskret (1 − R − T = Ringabsorption auf 1e-13).
  Tests `V22Volume`, `V22NearField`; V22c gibt für Silber zusätzlich σ_abs aus. M20 bitgleich.
- M21 V22c (Plan v2.41): quadratische und hexagonale Säule gegen den kartesischen Pfad mit PML (M18), beide Laurent;
  lokal ein grober Vergleich (`V22c.SquarePillarCoarse`), die Abnahme bei 1e-3 und die Silbersäule als
  Konvergenzstudie nur mit `FG_SERVER_TESTS=1` (`V22c.PillarsAgainstCartesianPml`).
- M21, Löser je Symmetrieklasse (Plan v2.40): `"cylindrical": {"symmetry": {"rotation": P, "L": n}}` mit Schichten
  `polygon2d` und `contour2d` um die Achse (Symmetrieprüfung, Fehlermeldung mit der gefundenen Ordnung). `CylBasis`
  rechnet eine Klasse m₀ + lP (Ringschichten blockdiagonal, Querschnitte über `coupledQ`), ebene Wellen mit der Quelle
  je Klasse (`fg::cyl::coupledSource`), Strahlen je Klasse; Ausgabe je Harmonischer, Python `rotation`, `class_L`.
  M20 bitgleich. Tests V22a (Klassenpfad = M20, 3e-13), V22d (Auswahlregeln, Drehung des Einfallsazimuts), V22e
  (Spiegelsymmetrie), `V22Config`.
- M21, erster Schritt (Plan v2.39): Winkelkoeffizienten χ_q(r) von Querschnitten um die Achse
  (`fg::cyl::arcsInside`, `angularCoefficients`, Plan eq:cylP_chi), kritische Radien (`criticalRadii`),
  Prüfung der P-zähligen Drehsymmetrie und Symmetrieordnung (`hasRotationSymmetry`, `rotationOrder`). Unit-Tests
  `test_cyl_angular` (n-Ecke und Kreise gegen exakte Bögen 1e-13). Gekoppelte Laurent-Blöcke einer Symmetrieklasse
  (`fg::cyl::coupledQ`, `HarmonicMap`, Plan eq:cylP_Q) und V22b (exzentrischer Faserkern = zentrierter Kern,
  1,6e-10 bei L = 7). Die Bernstein-Hilfen aus M24 liegen jetzt in `fg/geometry/bernstein.hpp` (M24 bitgleich).

### Geändert
- Plan v2.38 (Planfehler V21d): Gebundene Moden des verlustfreien planaren Hintergrunds nahe der Abschneidung lagen zu
  dicht am deformierten Weg (Sinusprofil, an k_a flach) und wurden bei kleinem R_s nicht aufgelöst; das Kontaktloch
  (Film n = 2,1) wanderte mit `dk_factor` um 10 %. Der Löser bestimmt jetzt die TE/TM-Moden des Hintergrunds
  (`cyl::planarBoundModes`, Schichtmatrizen, Bisektion) und legt gestufte Panels um sie (`SpectrumOptions::boundModes`,
  Breite = Abstand des Pols zum Weg). Ohne geführte Moden im Hintergrund ist das Gitter unverändert; Golden Files
  bitgleich. Tests `CylSpectrum.PlanarBoundModesOfAFilm`, `CylSpectrum.GradedPanelsAroundBoundModes`.

### Hinzugefügt
- M20 V21d (Plan v2.38): Säule und Kontaktloch gegen den kartesischen Pfad mit PML (M18), gesamter Streuquerschnitt
  (`V21d.*`); der Vergleich des Kontaktlochs braucht N ≥ 14 kartesisch und läuft nur mit `FG_SERVER_TESTS=1`
  (Zielsystem).
- `efficiencies.csv` enthält bei PML die Zeile `S_lat` (seitlicher Fluss des Streufeldes, bisher nur in C++/Python).

### Hinzugefügt
- M24 (Plan v2.37, Issue #30): `profile1d` mit `"shape": "nurbs"`: eine offene, periodische NURBS-Linie in 1D-Gittern
  (Überhänge erlaubt, kein Selbstschnitt), unter ihr `below`. Bleibt eine Schicht: Die kritischen Höhen (y-Extrema
  der rationalen Bézierstücke, Knoten) sind interne Bandgrenzen für RK und Tschebyschow (`StackLayer::yBreaks_um`,
  `Sectioned1DLayer::yBreaks_um`); JSON-Pointer späterer Schichten bleiben. `depth_um` optional. Fehler mit Pointer
  und Kennung `NURBS_INPUT`, `NURBS_NOT_PERIODIC`, `NURBS_PERIOD`, `NURBS_RANGE`, `NURBS_SELF_INTERSECTION`,
  `NURBS_STACK`. Tests `NurbsProfile.*` (V25a), `V25b`–`V25g`, `CAD_V25`, Python `test_profile1d_nurbs.py`; Beispiel
  `examples/profile1d_nurbs.json` mit neuem Golden File (alle anderen bitgleich).

### Behoben
- M23: `geometry::lineCrossings` verlor eine innere Nullstelle, wenn die Kurve die Gerade genau am Ende eines
  Bézierstücks erreichte (`refine` gab das Ende zurück); gefunden mit M24 (Schnitt genau auf der Höhe der Naht).
  Test `LineSection.RootAtTheEndKeepsTheInteriorRoot`; Golden Files unverändert.

### Geändert
- Nurbforge-Pin `FG_NURBFORGE_TAG` v0.6.0 → v0.12.1 (Plan v2.26): exakte Randkurven (Kreise statt kubischer
  Interpolanten), Boolesche Differenzen mit koplanaren Flächen (`overlap: "priority"` bei Körpern gleicher Höhe, Test
  `CAD.PriorityWithCoplanarLids`). Golden File `cad_line_grating_1d.csv` neu erzeugt (Grund: die Schnittkanten der
  gerundeten Linie liegen jetzt auf 8e-17 µm statt 3,9e-10 µm auf der exakten Geometrie; Effizienzen ändern sich um
  höchstens 8e-11); alle anderen Golden Files bitgleich.

### Behoben
- Der Wechsel deckte einen Fehler in Nurbforge auf: `critical_heights` verlor Extrema auf Abtastpunkten (exakte Kreise
  ab 0.7.0), der Hüllquader eines Zylinders fiel quer zur Achse zusammen und die Abflachung rekursierte mit Toleranz 0
  (`CAD.TwoMaterialsAndCavity` hing). Behoben in Nurbforge 0.12.1 (Issue #83, PR #84). `Contour2D::flatten` lehnt eine
  Toleranz ≤ 0 oder NaN jetzt mit einer Fehlermeldung ab (Test `Contour.FlattenRejectsANonPositiveTolerance`).

### Geändert
- Plan v2.31 (Planfehler V21f): Die Ergebnisse hängen von der Tiefe δ des deformierten Weges innerhalb des
  Diskretisierungsfehlers ab (die Knoten sind die Galerkin-Diskretisierung), nicht auf 1e-8; das Kriterium ist
  entsprechend angepasst. `cylindrical.spectrum.deform.delta` muss > 0 sein (ein negativer Wert wurde stillschweigend
  als reeller Weg gerechnet; in der C++-Schnittstelle bedeutet δ < 0 jetzt die falsche Seite, für V21f).

### Hinzugefügt
- M20 V21g (Plan v2.36): qualitativer Vergleich mit Popov/Bonod (phys. stat. sol. (b) 244, 3463, 2007; Quelle
  geprüft): Randspitzen in Polarisationsrichtung (Abb. 4), zylindrisches Plasmon |E_y| ∝ |H₁(k_p x)| längs der
  Oberfläche (Abb. 6), Fernfeld eines 10-nm-Lochs (Abb. 7: Silber gerichteter als Aluminium in der E-Ebene, H-Ebene cos²θ) (Tests `V21g.*`). Befund: Bei kleinen Löchern in Metall konvergiert die
  Amplitude erst bei k_max R ≈ 15 (50 nm: κ_max ≈ 24).

### Behoben
- M20 (Plan v2.36, Planfehler): Ebene Welle auf feinen Gittern: Lag der Knoten k ≈ 0 einer Ringschicht näher als
  1e-8 max|λ| am Exponenten des planaren Hintergrunds, integrierte der Löser die Scheibe mit RK (Minuten bis Stunden,
  im Nahfeld abgelehnt). Jetzt nimmt die dividierte Differenz des Strahls diese Terme; kartesisch bleibt der
  RK-Rückfall (Test `V21g.BackgroundExponentNextToTheNodeKZero`: gegen ein Gitter ohne Resonanz σ auf 5e-14).

- M20 V21i, Stufe 2 (Plan v2.35): Gradienten für ebene Wellen (Streufeld). Zielfunktion auf den Streuquerschnitten
  nach oben/unten (`("cylindrical", "R")`/`"T"`, je Harmonischer mit m); Parameter Ringradien, Ring- und
  Schichtmaterialien, Dicken (nicht die Halbräume). Die Hintergrundlösung der Einfallsordnung geht wie in M19 (B2)
  über zentrale Differenzen ein, alles Übrige analytisch. Gegen Differenzenquotienten auf ≤ 1,2e-7
  (`CylAdjoint.PlaneWaveSource`, `V21i.PlaneWave`).
- M20 V21i, Stufe 1 (Plan v2.34): Gradienten zylindrischer Strukturen (Gauß-Strahl, Totalfeld, RCWA) durch den
  diskreten Adjungierten. Zielfunktion R und T gesamt oder je Harmonischer (`("cylindrical", "R")`,
  `("cylindrical", "T", m)`); Parameter Ringradien `/layers/k/radii_um/i`, Dicken, Materialien (nicht die
  Halbräume). Neue Rückwärtsregeln für Fluss, Propagationsmatrix, Galerkin-Matrizen und die Faktorisierungen.
  Gegen Differenzenquotienten auf ≤ 2,9e-7 (Tests `CylAdjoint.*`, `V21i.*`).

### Behoben
- M20 (Plan v2.33): Ein Plasmon des planaren Hintergrunds zwischen den Verzweigungspunkten (z. B. Silber/Luft bei
  einem Silberfilm auf Glas) ist eine schmale leckende Resonanz auf dem reellen Weg. Mit gleichmäßigen Panels
  konvergierten R und T in `dk_factor` nicht (T zwischen 0,42 und 0,54 bis `dk_factor` 6). Jetzt setzt das Gitter
  gestufte Panels um Re k_p (`SpectrumOptions::resonances`), und `dk_factor` 1 ist auf 2e-5 konvergiert (Tests
  `CylSpectrum.GradedPanelsAroundAResonance`, `V21h.LeakyPlasmonBetweenTheBranchPoints`). Auch der reelle Weg
  (`deform: false`) wird genauer (V21f). Fälle ohne solches Plasmon sind bitgleich.

### Hinzugefügt
- M20 V21h Konvergenzstudie (Plan v2.32): `scripts/v21h_convergence.py` schreibt
  `docs/validation_data/v21h_convergence.csv`, Tabellen in `docs/validation.md`. Dielektrische Kugel: κ_max bestimmt
  den Fehler, `dk_factor` unkritisch, Harmonische ab m_cut = 4 gesättigt. Befund: eine Lochblende im Silberfilm mit
  Luft darüber konvergiert in `dk_factor` nur langsam (leckendes Plasmon Silber/Luft zwischen den
  Verzweigungspunkten, auf dem reellen Weg); mit Glas auf beiden Seiten auf 4e-8. Abhilfe offen (Plan v2.32).
- M20 V21f (Plan v2.31): m ↔ −m 1e-12, Lochblende im Silberfilm unabhängig von δ innerhalb des
  Diskretisierungsfehlers, reeller = deformierter Weg bei stark verlustbehaftetem Metall (6,5e-8), falsche
  Ausweichrichtung springt (`V21f.SymmetryAndIntegrationPath`).
- M20 Winkelspektrum und Fernfeld (Plan v2.30): je Harmonischer und Seite dP/dθ an den Knoten des Spektrums, exakt
  getrennt nach s = (e₊ + e₋)/2 und p = (e₊ − e₋)/2, mit Quadraturgewicht (Summe = Leistung der Harmonischen);
  Fernfeld dP/dΩ auf (θ-Knoten) × `cylindrical.far_field.n_phi` Azimute. Strahl bezogen auf P_inc, ebene Welle als
  differentieller Streuquerschnitt (µm²/sr). Python `cylindrical["angular_spectrum"]`, `["far_field"]`; `fgsolve`
  schreibt `angular_spectrum.csv` und `far_field.csv`. Tests `V21FarField.*`: Leistung je Harmonischer auf 6e-17,
  φ-Mittel = Winkelspektren, kleine Kugel gegen Mie (S1, S2) auf 4,7e-3 des Maximums bei κ_max = 8.
- M20 CAD-Rotationskörper (Plan v2.29): `cad` mit `"axisymmetric": true` im Zylinderstapel (dort Pflicht, sonst
  unzulässig). `cell_origin_um` ist die Achse (Standard: Mitte des Hüllquaders); jeder Schnitt muss aus Kreisen um
  die Achse bestehen (Scheiben oder Ringe), sonst neuer Befund `NOT_AXISYMMETRIC` mit Höhe und Radiusbereich. Bänder
  `homogeneous`, `rings` (Extrusion) oder `rings_profile` (Radien aus den Schnitten in den Scheibenmitten). Bericht mit
  `"dimension": "axisymmetric"` und `axis_um`. Tests `CAD_M20.*`: Zylinder + Kegelstumpf + Kugelkappe = Stapel aus
  `rings`/`rings_profile` auf 2e-13, Rohr, Befunde, Konfigurationsfehler.
- M20 Schichttyp `rings_profile` (Plan v2.28): Radien als Funktionen der Höhe (Zahl, `cone`, `sphere`, `table`;
  C++ auch Funktion), gerechnet als Treppe aus `slices` `rings`-Schichten mit den Radien in der Scheibenmitte (Ringe
  ohne Breite entfallen, kreuzende Radien sind ein Fehler). Tests `V21RingsProfile.*`: Treppe = explizite
  `rings`-Folge bitgleich (Kegel, Kugel von V21c), Kugelkappe, Treppenfehler zweiter Ordnung.
- M20 Strahl als Streufeld (Plan v2.27, V21e): `"scattered_field": true` gilt jetzt auch für Zylinderstapel mit
  Gauß-Strahl. Hintergrund = Stapel mit dem äußeren Material jeder Ringschicht, auf dem Hankel-Gitter gelöst und je
  Schicht als 4 N_k Exponentiallösungen dargestellt; Quelle σ = (M − M_bg) F_bg, partikuläre Lösungen über die Moden
  der Scheibe, fast resonante Terme (Knoten k ≈ 0) in der Form der dividierten Differenz. R, T, A_Ringe und Nahfeld
  gleich dem Totalfeld auf 1e-13 bzw. 1e-12, RK auf 5e-9 (`V21e.BeamScatteredFieldEqualsTotalField`). Der
  kartesische Pfad bleibt bitgleich.
- M20 Nahfeld und Querschnitte (Plan v2.25): `nearField` für Zylinderstapel (horizontaler Schnitt x–z, vertikaler
  Schnitt x–y bei `vertical_z_um`; `nearfield.extent_um`, Standard 1,5 R_s). Spektren je Harmonischer aus den
  Zustandsvektoren, Synthese mit Bessel-Tabellen, kartesische Komponenten nach eq:cyl_pm; bei Li in r
  FFF-konsistent über die stetigen G± (`CylQ::A`, `PlaneWaveSource::Apm`), D_r an Ringwänden stetig. Ebene Wellen:
  Hintergrund im Ortsraum, Total- oder Streufeld. Punkte jenseits r_max eines deformierten Weges sind NaN mit
  Warnung (`NearFieldResult::warnings`, `fgsolve` und Python `near_field(...)["warnings"]`).
- Volumenintegrale über die Ringe: Absorptionsquerschnitt σ_abs, σ_ext = σ_sca + σ_abs und das optische Theorem in
  Volumenform als Kontrolle (CSV `C_Sabs`, `C_Sext`, `C_Sext_ot`; Python `absorption_um2`, `extinction_um2`,
  `extinction_theorem_um2`); bei Strahlen auf verlustbehafteten Ringen `absorbed_rings`. Die Volumenintegrale laufen
  nur bei verlustbehafteten Ringen oder mit `"cylindrical": {"optical_theorem": true}` (etwa doppelte Rechenzeit);
  verlustfreie ebene Wellen rechnen so schnell wie bisher. Neue Bausteine
  `fg/cyl/field.hpp` (`synthesis`, `besselTable`, `radialRule`, `fieldsFromG`).
- Tests `CylField.*` (Gauß-Paar, radiale Regel, 𝔄 und â am Knoten, G±), `V21NearField.*` (Brennebene = Gauß-Profil
  2e-10; A_Ringe = 1 − R − T: Laurent 4e-14, Li zweiter Ordnung; D_r an der Wand stetig), `V21c.OpticalTheoremAndGoldSphere`.

### Geändert
- Zustandsvektoren an beliebigen Höhen (`stateSamples`) und Hintergrundfeld (`backgroundAt`) sind aus dem kartesischen
  Nahfeld ausgelagert und basisunabhängig; das LU von `G` wird je Element einmal statt je Höhe berechnet. Kartesische
  Ergebnisse und Golden Files bitgleich.
- `cylindricalSolve` teilt Gitter, Schichten und Anregung mit dem Nahfeld (`CylContext`); Strahl-Ergebnisse unverändert.

## [0.20.0] – M23 (CAD-Import für 1D-Gitter), Überlappungsprüfung in 2D, Agentenreferenz

### Hinzugefügt
- M23 CAD-Import für 1D-Gitter (Plan Teil 9, v2.22–v2.23): Eine `cad`-Schicht in einem 1D-Stapel schneidet
  z-invariante Körper einmal bei `z_cut_um` und bildet je Höhe Intervalle mit Materialien und Normalen
  (exakte Kreisbögen bleiben exakt). Neuer interner Schichttyp `Sectioned1DLayer` (RK, RCWA, Tschebyschow, ASR,
  Nahfeld, Gradienten), Normalenfeld mit linear interpoliertem Winkel (eq:cad1d_N), Bänder `homogeneous`,
  `lamellar1d`, `sectioned1d`, Befund `NOT_Z_INVARIANT`, `cell_origin_um` in 1D als Zahl oder `[x0]`,
  `fg.cad.preview(step, cfg)` ohne `y_um` für das Profil. Nähte (Restschicht, Körper über dem Zellrand) ohne
  Grenzfläche; `periodic_union` und `priority` als Intervall-Operationen. Formgradienten über `StackParameter`
  auch für `sectioned1d`. PML mit `sectioned1d` wird abgelehnt.
- Tests `test_line_section`, `test_sectioned1d`, `CAD_V24` (V24a–h), Python `test_line_grating_1d` (V24i).
- Beispiel `examples/cad_line_grating_1d.json` (gerundete Linie aus der Meldung zu v2.21) mit STEP aus
  `scripts/make_cad_examples.py`; neues Golden File `cad_line_grating_1d.csv`, bestehende unverändert.
- Plan v2.22 (Teil 9: Spezifikation M23), v2.23 (Umsetzung, Abweichungen: `FEATURE_SMALL` mit dem Rasterabstand
  wie in 2D, PML nicht unterstützt).
- Kompakte Referenz für Sprachmodell-Agenten `python/fouriergrating/LLM_REFERENCE.md` (englisch, ca. 4000 Tokens):
  Konventionen, alle Konfigurationsschlüssel, Python-API, Rezepte, CAD, Fallstricke. Wird mit dem Wheel
  ausgeliefert; Abruf mit `fg.llm_reference()` oder `python -m fouriergrating` (neu `__main__.py`). Test
  `tests/python/test_llm_reference.py`: jeder öffentliche Name, Schichttyp, vom Parser gelesene Schlüssel und
  CAD-Prüfcode steht darin, jeder `json`-/`python`-Block läuft, die Datei liegt im installierten Paket.
- Release 0.20.0: Version, README-Status, Handbuch (Programmstand, Plan v2.24).

### Geändert
- Die Tschebyschow-Interpolation in y ist eine gemeinsame Klasse `ChebyshevQ` (`SectionedModel` bitgleich).
- 1D-Stapel lehnen `cad`-Schichten nicht mehr ab (v2.21 → M23).

### Behoben
- Überlappungen in 2D (Plan v2.24): Körper bzw. Regionen mit kollinear verlaufenden Rändern (z. B. zwei Quader
  gleicher z-Ausdehnung) wurden weder vom Slicer noch von `contour2d`/`polygon2d`/CAD im Löser als überlappend
  erkannt; berührende Körper galten je nach Startpunkt der Schleife als überlappend. Jetzt: echte Kreuzung mit
  Rundungstoleranz (`geometry::properCrossing`) oder Probepunkt knapp innerhalb einer Randkante
  (`Region::interiorProbes`); eine Region, die ihr eigenes Bild entlang des Randes berührt, wird weiter abgelehnt
  (`geometry::sharedBoundary`).
- `periodic_union` in 2D (Planfehler, v2.24): Die Vereinigung eines Körpers mit seinem überlappenden Bild überlappt
  wieder ihr Bild; nach vier Runden bestand die Prüfung trotzdem, erst der Löser brach ab. Jetzt `PERIODIC_OVERLAP`
  mit Hinweis; `periodic_union` wirkt nur in 1D.
- `overlap: "priority"`: höchstens n(n−1)/2 + 2 Runden mit abschließender Prüfung (vorher vier ohne); scheitert die
  Boolesche Differenz in Nurbforge (koplanare Flächen, offene Ergebnishülle), meldet der Slicer `BODY_OVERLAP` statt
  einer Ausnahme.
- Tests `CAD.PeriodicUnionCannotResolveSelfOverlapIn2D`, `CAD.PriorityResolvesEveryOverlap`,
  `CAD.OverlapWithCollinearEdgesIsFound`.
- CAD-Slicer (M14, Plan v2.21), gemeldet beim Import eines z-invarianten Liniengitters aus Nurbforge:
  - Hüllquader aus den getrimmten Körpern (Extrema von `nf::critical_heights` entlang x, y, z) statt aus den
    Kontrollpunkten der Trägerflächen. Ein aus einem größeren Rechteck getrimmter Deckel verschob vorher
    `y_range_um` und `cell_origin_um` (zusätzliche Hintergrundbänder unter und über dem Körper).
  - 1D-Stapel: Die `cad`-Schicht wird mit klarer Meldung abgelehnt. Vorher war `dz = 0`, jeder Körper wurde mit sich
    selbst als „periodischem Bild“ verglichen (`PERIODIC_OVERLAP` je nach Startpunkt der Schleife), die Flächenanteile
    waren NaN, und `periodic_union` vereinigte identische Shells (hing).
  - Ein Körper, der sein eigenes periodisches Bild entlang eines Randstücks berührt (z. B. eine über genau `dz`
    extrudierte Linie), meldet `PERIODIC_OVERLAP` (auch mit `periodic_union`). Vorher bestand er die Prüfung, und
    erst der Löser brach mit „regions overlap its periodic image“ ab.
  - Tests `CAD.BoundsOfTrimmedLids`, `CAD.ZInvariantLine`.

## [0.19.0] – M22 (Dispersive Materialien), Teile von M20, Benutzerhandbuch

### Hinzugefügt
- Release: `wheels.yml` baut bei einem Tag `v*` zusätzlich das sdist, prüft Tag = Version und hängt Wheels und
  sdist an ein GitHub-Release; Upload zu PyPI (Trusted Publishing, Umgebungen `pypi` und `pypi-cuda12`)
  von Hand per `workflow_dispatch`. `fouriergrating-cuda12` verlangt genau dieselbe Version von `fouriergrating`
  (Test `tests/python/test_release_metadata.py`).
- Benutzerhandbuch in LaTeX (`docs/handbuch/`, ca. 100 Seiten): Grundlagen und Konventionen, Installation, Schnellstart,
  vollständige Konfigurationsreferenz, alle Schichttypen, Löserwahl/Faktorisierung/Konvergenz (mit gemessener
  Konvergenztabelle V19), Scans und Sweeps, Nahfeld, Ausgaben, ASR/PML/Zylinder, CAD, Python, Gradienten, Leistung,
  C++, GUI, Rezepte mit echten Ergebnissen, Fehlersuche, Anhänge. `scripts/build_docs.sh` baut es mit.
  Dazu dispersive Materialien (M22), eine vollständige Befehlsreferenz im Anhang (`fgsolve`, `fgbench`, Skripte,
  Python-API, Umgebungsvariablen, CMake-Optionen und Presets) und ein Stichwortverzeichnis (Konfigurationsschlüssel
  automatisch über `\key`, `makeindex` in `build_docs.sh`). Test `tests/python/test_handbook_reference.py`: jeder
  öffentliche Python-Name, jede Option von `fgsolve`/`fgbench`, jede Umgebungsvariable, CMake-Option und jedes
  Skript steht in der Referenz. Behoben: `"auto"` erschien als `äuto"` (babel-Kürzel), `--profile` als `-profile`.
- Doku korrigiert (`config_format.md`, `user_guide.md`): `triangle` hat `apex_x_um` (nicht `blaze_deg`/`apex_deg`),
  `points` ist ein stückweise lineares Polygon (kein Spline), CSV-Spalten heißen `phi_deg` usw., `diagnostics.json`
  wird nicht geschrieben, `nearfield.quantities` ist derzeit wirkungslos, `breakpoints_um` wird nicht gegen die
  Periode geprüft.
- M22 dispersive Materialien (Plan Teil 8, v2.19): Materialwerte unter `n`, `eps` und in `eps_principal` dürfen
  Modelle sein – `sellmeier` (optional `range_um`), `cauchy`, `drude_lorentz`, `table` (`file`/`data`, Einheit
  um/nm/eV, `linear`/`cubic` Not-a-knot, `extrapolate`); Modul `fg::Dispersion` (`include/fg/geometry/dispersion.hpp`),
  `Material::at`/`Stack::at`, `detail::resolveStack` an jedem Einstieg (`solve`, `solveScan`, `fieldCoefficients`,
  `nearField`, `orderSMatrix`, Gradienten, Zylinderpfad), nicht aufgelöste Materialien werden abgewiesen. Frühe
  Prüfung beim Laden (Bereich, Sellmeier-Pol, verlustfreies Superstrat, jeder Sweep-Wert) mit JSON-Pointer.
  `Simulation` baut bei einem Wechsel von λ nur Schichten mit dispersivem Material neu. Ausgaben
  `Result::diag.materials`, `fgsolve`-Zeilen `material …`, `"output": {"materials": "materials.csv"}`,
  `fg::io::materialAt`, Python `fg.material_eps`, `diagnostics["materials"]`. Gradienten nach Modellwerten
  (eq:disp_grad_p) und nach λ mit Kettenregel (eq:disp_grad_lambda). Python
  `fouriergrating.materials.from_refractiveindex` (YAML, optional PyYAML: `fouriergrating[materials]`),
  `load_config` löst Tabellendateien relativ auf. Tests `test_dispersion` (V23b, V23g), `test_validation_v23`
  (V23a, V23c–f), `tests/python/test_materials.py` (V23h); Beispiele `examples/dispersive_metal_grating.json`
  (neues Golden File, keine bestehenden geändert) und `examples/python/dispersive_scan.py`.
- Plan v2.19: Umsetzungsdetails zu M22 (Abschnitt „Umsetzung (v2.19)“): Sellmeier `range_um` für den Import aus
  refractiveindex.info, Gradienten über die Materialgradienten von M19 mit aufgelöster Konfiguration,
  CAD-Schichten im Zwischenspeicher, Inhalt von `materials.csv`.
- S-Matrix des Stapels in der s/p-Basis der Ordnungen (Plan v2.17, Abschnitt „S-Matrix in der s/p-Basis der
  Ordnungen“, eq:Ssp): `Solver::orderSMatrix`, `Simulation::orderSMatrix`, Python `Simulation.smatrix` mit
  `R_top`, `T_down`, `T_up`, `R_bottom`. Ein Eigenproblem liefert beide Polarisationen und Licht von oben und von
  unten; andere Spalten beschreiben die Bloch-Familie bis auf den Trunkierungsfehler (das Ordnungsfenster liegt um die
  Einfallsordnung); im 1D-in-plane-Fall beide Polarisationsblöcke. Nur periodisch, homogenes isotropes
  Superstrat und Substrat, ohne Koordinatentransformation. Tests `test_order_smatrix` (zentrale Spalte = `solve`
  auf 1e-12, Spiegelsymmetrie für Licht von unten, Simulation = Löser) und
  `test_smatrix_central_column_equals_solve`. Für RTT-Fusion L4 (#3).
- M20 (ebene Wellen): Einfall `plane_wave` für Zylinderstapel in der Streufeldformulierung (Hintergrund wie M18,
  Quellen nach Jacobi–Anger, `fg::cyl::besselProjection`/`planeWaveSource`), Harmonische `auto` nach m_cut,
  Streuquerschnitte `C_Sup`/`C_Sdown`/`C_Ssca` (µm²) und je Harmonischer; V21c (Mie-Kugel −5,7e-4 nach
  Extrapolation), Symmetrie- und RK-Tests ebener Wellen.
- Plan v2.16: konsistente Quelle ebener Wellen bei Li in r, p = b̂ − Q â (der Plan hielt sie wegen des Deltaspektrums
  für unmöglich; mit ε_bg = ε∞ heben sich die Deltaanteile auf). Mit der exakten Quelle konvergiert σ_sca der
  Mie-Kugel nur erster Ordnung in k_max (+1,5 % bei κ = 8), mit der konsistenten bleibt der Treppenfehler.
- M20 (Löser): Schnittstelle `detail::ModalBasis` der Stapelassemblierung (kartesisch bitgleich), zylindrische Basis
  je Harmonischer, Schichttyp `rings`, Konfigurationsblock `"cylindrical"` mit Strahlanregung (`gaussian_beam`),
  Ausgabe `Result::cylindrical` und CSV-Seiten `C_R`, `C_T`, `C_A`, `C_Rm`, `C_Tm`; V21a; Beispiel
  `examples/cyl_pillar_beam.json`.
- M20 (Bausteine): Modul `fg::cyl` – Bessel-Funktionen ganzzahliger Ordnung mit komplexem Argument, diskretes
  Hankel-Spektrum (`spectralGrid`), Hankel-Galerkin-Matrizen (Lommel-Integrale, Quadratur), Faktorisierung Laurent und
  Li in r (`cylQ`), Propagationsmatrix, Modenbasis homogener Medien und Poynting-Fluss je Harmonischer; Unit-Tests
  `test_cyl`, Validierung V21b (Fasermoden: Li in r ≤ 1,1e-4 bei k_max = 16, Laurent erster Ordnung; Werte des
  Prototyps reproduziert).
- Plan v2.18 (Teil 8): Spezifikation dispersiver Materialien – M22 (Tabellen, Sellmeier, Cauchy, Drude–Lorentz in
  `n`, `eps`, `eps_principal`; Auswertung je Scanpunkt, frühe Prüfung mit JSON-Pointer, Gradienten nach λ und
  Modellwerten, Import aus refractiveindex.info); Roadmap und AGENTS.md (M22 darf vor M20/M21) ergänzt.
- Plan v2.15: Umsetzungsdetails zu M20 – Wurzelsubstitution je Teilpanel (Panel zwischen zwei Verzweigungspunkten
  geteilt, die Sinus-Substitution über das ganze Panel ist mit 8 Knoten 10³-fach ungenauer), Bessel nach Miller mit
  Genauigkeitsmaß relativ zur Amplitude (scipy bei |z| ≈ 500 nur auf ≈ 1e-14), Lommel-Toleranz 5e-14.
- Plan v2.14 (Teil 7): Spezifikation der Zylinderkoordinaten – M20 rotationssymmetrische Strukturen
  (Fourier–Bessel-Spektrum, Galerkin-Matrizen, Li in r über stetige Hilfsgrößen, Strahlen, Streufeld, Fernfeld) und
  M21 P-zählige Drehsymmetrie; Prototyp `scripts/proto_cylindrical.py` (Fasermoden, exzentrischer Kern);
  Roadmap ergänzt.
- M19 Phase 2: Gradienten nach Wellenlänge, Einfallswinkeln, Polarisationswinkel und Brechzahl/Permittivität von
  Superstrat und Substrat (`/incidence/...`, `/materials/<superstrat|substrat>/...`); Rückwärtsregeln
  `adjoint::isotropicModes`, Propagationsmatrix nach α, γ, Scheiben-S-Matrix nach der Referenzbasis.
- `"checkpoint_every": k` (`SolverOptions::checkpointEvery`): Gradienten speichern nur jedes k-te Präfixprodukt
  und keine Scheibendaten; rückwärts je Abschnitt neu gerechnet, Ergebnis bitgleich.
- Gradienten mit anisotropen Materialien (`eps_tensor`, `eps_principal`, `euler_deg`): Rückwärtsregel
  `adjoint::qAnisotropic` (Laurent, FFF_A, FFF_B); homogene Schichten, Lamellen, Profile, 2D-Querschnitte (Laurent).
- Gradienten mit Koordinatentransformationen (ASR, x und z): Lamellen, Profile, 2D-Querschnitte, anisotrop, RK;
  `detail::transformedCoordinates` gemeinsam mit dem Löser.
- Gradienten anisotroper 2D-Querschnitte mit FFF_B.
- Gradienten nach λ, Winkeln und Superstrat/Substrat auch mit Koordinatentransformation; Rückwärtsregel der
  Propagationsmatrix nach dichten A, G.
- Gradienten in `profile2d`-Schichten (isotrop und anisotrop): Materialien, `depth_um`, `x0_um`, `z0_um`.
- Plan v2.13: Spezifikation der Gradienten aperiodischer Strukturen (Streuflüsse oben/unten, Modenleistungen;
  Stufen A–C).
- Gradienten aperiodischer Strukturen bei Normaleinfall (Stufe A): `Objective::scattering` (`ScatterRef` Up/Down),
  Python `{("scattering", "up"): w}`; bei modalem Einfall Leistungen der geführten Moden (`Objective::modes`,
  `ModeRef`; Python `{("mode", "R", j): w}`).
- Gradienten mit der Streufeldformulierung (schräger Einfall mit PML, Stufe B, RCWA und RK); Rückwärtsregel
  `adjoint::starSources` des Sternprodukts mit Quellvektoren.
- V20f: Reziprozität der Streuamplituden (Winkelspektrum des Streufeldes, 5e-2 bei L = 8 µm).
- Python-CI mit PyTorch, JAX und scipy; NURBS-Beispiel als Test.

### Geändert (Standard)
- `backend_threshold` (`backend = auto`): 160 statt 1000 (Messung auf der RTX 6000 Ada, `docs/benchmarks.md`); GPU = CPU im Rundungsbereich.

### Behoben
- `CylPropagationMatrix.HomogeneousMediumEigenvectors` auf arm64 (FMA-Kontraktion): Toleranz relativ zu
  max|k² − ε| statt absolut.
- Build mit Apple clang 21 (macOS, libc++): Die Skalarprodukte `ph.transpose() * X` im Nahfeldschnitt
  (`solvers/solver.cpp`) werden jetzt mit `.value()` ausgewertet; Eigen 3.4 liefert dort `operator const Scalar()`,
  das clang nicht an `push_back(T&&)` bzw. eine Zuweisung an `std::complex` bindet. Ergebnisse unverändert.
  Die Tests, die `std::cyl_bessel_j`/`std::cyl_bessel_k` als unabhängige Referenz nutzen (fehlen in libc++), werden
  dort mit Begründung übersprungen (`tests/support/std_bessel.hpp`); mit libstdc++ und MSVC laufen sie wie bisher.
- CUDA: `cusolverDnXgeev` lieferte bei gleichzeitigem Aufruf aus mehreren Threads `CUSOLVER_STATUS_INTERNAL_ERROR` (7)
  (RTX 6000 Ada, CUDA 12.9, `FG_NUM_THREADS` >= 2, `FG_GPU_STREAMS` >= 2). Der Aufruf ist jetzt mit einem Mutex
  serialisiert (`CudaBackend::eigT`); Test `Backend.CudaConcurrentEigenSolves`.
- `NormalFieldSamples.PolygonsBitIdentical*` auf GCC mit FMA (x86-64-v3, AlmaLinux 10): `-ffp-contract=off` für
  `geometry/normal_field.cpp`; dadurch liegt `V8_2D.ZInvariantAnisotropicEquals1D` (1D- gegen z-invarianten 2D-Code) bei
  1,6e-12 statt unter 1e-12 (Rundung), Schranke dort 5e-12.
- M10-Abnahme auf der RTX 6000 Ada: Messwerte in `docs/benchmarks.md`.

### Geändert
- Plan v2.20: Literatur zu M22 [D1]–[D4] bibliographisch geprüft (Seitenbereiche ergänzt; refractiveindex.info
  unter CC0 1.0). Drude–Lorentz-Werte für Ag und Al nach Rakić geprüft: `fg.material_eps` gibt die Rakić-Tabellen
  von refractiveindex.info auf die 5 Stellen der Tabelle wieder. Al-Beispiel in Plan und `docs/config_format.md` als
  gekürzt gekennzeichnet, Pfade zu refractiveindex.info an den heutigen Aufbau der Datenbank angepasst.
- Gradienten nach Formparametern mit FFF in 2D: das Rasterfeld von Algorithmus NF bewegt sich mit der Geometrie
  (Differenzenquotient über eine Rasterzelle). Eingefroren war z. B. die Verschiebung eines Kreises um 30–150 %
  falsch (nur der Radius war geprüft). Kosten: zwei NF-Rechnungen je Formparameter.

## [0.18.0] – M19 (Differenzierbarkeit)

### Hinzugefügt
- Rückwärtsregeln (`include/fg/adjoint/primitives.hpp`): Lösen, Inverse, Eigenzerlegung mit Entartung, Toeplitz,
  Sternprodukt, T → S, Scheiben-S-Matrix, Propagationsmatrix, Q-Faktorisierungen (Laurent, Li, FFF iso/sym/rig mit
  Normalenfeld-Kotangenten), RK-Schritt; Formableitung `shapeDerivative` (Plan eq:shape2d).
- `fg::gradient` (diskrete Adjungierte, RCWA und RK): Parameter per JSON-Pointer (Dicke, n/ε, Lamellenkanten,
  Polygon-/Contour2D-Geometrie, Sinusprofil) und `StackParameter` (C++, Schnittschichten und CAD-Körper).
- `Simulation::valueAndGrad` / `vjp`; Python `Simulation.value_and_grad` / `vjp`, `fouriergrating.torch`,
  `fouriergrating.jax`; Beispiele `optimize_littrow.py`, `optimize_nurbs_pillar.py`.
- RK-Transfers zeichnen auf Wunsch ihre Schritte auf (`RKSteps`); Schichtmodelle liefern `adjointSlice`.
- Tests `test_adjoint`, `test_validation_gradient`, `CAD_M19`, `tests/python/test_gradient.py`.

### Geändert
- Plan v2.11 (Umsetzung und Abnahme M19), Version 0.18.0.

## [0.17.0] – M18 (PML, isolierte Strukturen, Streufeld)

### Hinzugefügt
- `"aperiodic": {"pml": {"x"|"z": {"width_um", "beta"}}, "formulation", "excitation"}`: PML als komplexe Streckung
  (Plan eq:pml, `PmlMap`), Zonenprüfung (`checkPmlZones`, auch CAD und strukturierte Halbräume), Hintergrundstapel
  (`backgroundStack`).
- Streubreiten bzw. -querschnitte (oben, unten, seitlich) und Absorption aus der Flussbilanz der geschlossenen Box
  (`Result::scattering`, CSV `S_up`/`S_down`/`S_abs`, CLI, Python `result.scattering`).
- Streufeldformulierung (Plan eq:scat_source, eq:scat_ode) mit σ = (M − M_bg) F_bg: Hintergrundfeld (Airy),
  partikuläre Lösungen je RCWA-Scheibe, RK über das erweiterte System, erweitertes Sternprodukt (`SourcedSMatrix`),
  innere Felder und Nahfeld (gesamt oder `"field": "scattered"`); automatisch bei schrägem Einfall mit PML,
  `scattered_field` / `SolverOptions::scatteredField` für periodische Gitter.
- Strukturierte Superstrate/Substrate (`superstrate_structure`, `substrate_structure`) und Anregung durch eine geführte
  Mode (`"excitation": {"mode": k}` / `{"n_eff": n}`): Modenleistungen, Abstrahlung (`Result::modal`, CSV `MR`/`MT`/`Mrad`).
- Tests `test_validation_v20` (V20a–V20g) und `CAD_V20g` (Kontaktloch als CAD-Körper).

### Geändert
- Numerische Halbraumbasen (M17/M18): propagierende Moden werden bei |Re λ| ≤ 1e-6 |λ| nach dem Fluss klassifiziert.
- Plan v2.10 (Umsetzung und Abnahme M18), Version 0.17.0.

## [0.16.0] – M17 (Koordinatentransformationen, ASR)

### Hinzugefügt
- `"coordinate_transform": {"x"|"z": {"type": "asr", "eta", "breakpoints_um": "auto" | [...]}}` (Plan eq:asr):
  `CoordinateMap` (f, f′, f⁻¹, [[1/f′]], Umrechnungen eq:ct_rayleigh / eq:ct_inc), `Stack::transform`,
  `resolveTransform`.
- Propagationsmatrix, E_y, H_y mit dichten A = [[1/f′]] α, G = [[1/g′]] γ (Plan eq:ct_deriv, Ersetzungsregel).
- Schichtmodelle in X, Z: Lamellen und Profile über abgebildete Intervalle, Polygone/Konturen/CAD-Schnitte über das
  abgebildete Konturintegral, Normalenfelder (Profile per Quadratur, Algorithmus NF umgetastet), Profile 2D mit Linien
  bei z = g(Z).
- Löser: numerische Halbraumbasen (Klassifikation mit Fluss), Moden homogener isotroper Medien aus dem n_f-Block von
  P R (exakt entartete TE/TM-Paare), Einfall ohne evaneszente Anteile, und Effizienzen über die physikalischen
  Ordnungen, Nahfeld auf dem physikalischen Raster.
- Tests `test_coordinate_transform` (Unit) und `test_validation_v19` (η = 0 bitgleich, Matrixpfad = Diagonalpfad, V1–V3,
  V9, Nahfeld, V19 mit Konvergenztabelle in `docs/validation.md`).

### Geändert
- Plan v2.9: Quadratur statt DFT/Trapezregel, Abnahmekriterium für V3/V9 präzisiert (gemeinsamer Grenzwert statt 1e-8
  bei gleichem N). Ohne Transformation bitgleich (Golden Files unverändert).

## [0.15.0] – M16 (Wählbare Genauigkeit)

### Hinzugefügt
- `"precision": "auto" | "fp64" | "mixed" | "fp32"` (`FG_PRECISION`), `"precision_fallback"`, `"growth_limit"`;
  `linalg::Precision`, `setPrecision`, `gpuSelected`, `gpu()`, `product()`; `diagnostics.precision`,
  `precision_fallback`, `refinement_steps`, `refinement_fallbacks` (C++, Python, `fgsolve`).
- CPU-Backend in einfacher Genauigkeit (Eigen `MatrixXcf`) mit FP64-Nachiteration (Mixed); CUDA: Gerätematrizen in Z
  oder C, `Cgemm`/`Cgeam`/`Cdgmm`, `Xgetrf`/`Xgeev` mit `CUDA_C_32F`, `cusolverDnZCgesv` (Mixed), `CUBLAS_DEFAULT_MATH`
  (kein TF32); Host-API auf den residenten Operationen aufgebaut.
- `src/cuda/device_policy.cpp`: `auto` = FP32 auf einer GeForce RTX 3090, sonst FP64; CPU immer FP64.
- Angepasste Parameter: Wachstumsgrenze auto (FP64 1e8, sonst 4e3), `rk.tol` ≥ 1e-5 in FP32, Klassifikationstoleranz
  1e-4 bei C-Eigenwerten; Rückfall auf FP64 bei verlustfreien Strukturen mit |Σe − 1| > 1e-3.
- Tests `test_precision` (Backend-Operationen, Nachiteration bis cond 1e6, Solver-Modi, GPU), Studie
  `test_precision_study` (Tabelle in `docs/validation.md`), Laufzeiten in `docs/benchmarks.md`.

### Geändert
- Produkte Ψ⁻¹·(…) (Substrat, RK-Scheiben, Host-Pfad der Scheiben-S) über das Backend: GPU-Laufzeit halbiert.
- `CMAKE_CUDA_ARCHITECTURES` 86;89 (`cuda-release`), 86 (`win-cuda`); Konsistenzprüfung der Klassifikation 2τ.
- `fgsolve` gibt Warnungen und die Genauigkeit aus; GPU-Vergleichstests auf 1e-10 setzen `precision = fp64`.
- Ergebnisse mit Standardoptionen auf der CPU bitgleich (Golden Files unverändert). Plan v2.8.

## [0.14.0] – M15 (Profiling und Parallelisierung)

### Hinzugefügt
- `fg::prof` (`fg/core/prof.hpp`): `FG_PROF_SCOPE`, threadsichere Aggregation (Summe, Anzahl, Maximum), JSON;
  `fgsolve --profile` (`<output.dir>/profile.json`), Python `solve(cfg, profile=True)` / `solve_scan(..., profile=True)`
  (`diagnostics["timings"]["scopes"]`).
- Threadmodell: `"threads": {"scan", "slices", "blas"}` (`auto` = Scanebene zuerst, Rest auf die Scheiben),
  `FG_NUM_THREADS`, `fg::threadPlan`, `fg::parallelFor` (OpenMP, sonst `std::thread`), `innerThreads()`.
- Parallele Assemblierung: S-Matrizen aller Schichten/RCWA-Scheiben/RK-Stabilitätsscheiben parallel (blockweise,
  Speicher O(Threads)); Spektralradien der RK-Schichten parallel.
- `"reduction": "tree"`: Sternprodukt-Kette als Binärbaum fester Form (`fg::starTree`, `BasicStarChain`); bitgleich für
  jede Threadzahl, ≤ 1e-12 von der sequentiellen Kette.
- Parallele Geometrie: Raster von Algorithmus NF (zeilenweise), Gauß-Glättung, Konturkoeffizienten je Segment,
  Tschebyschow-Knoten (bitgleich).
- GPU-residente Kette (`fg/solvers/resident.hpp`): `linalg::DeviceMatrix` mit Speicherpool, residente Operationen im
  Backend (Host-Referenz in Eigen, CUDA über cuBLAS/cuSOLVER), Moden, Scheiben-S und Sternprodukte auf dem Gerät,
  `FG_GPU_STREAMS` (Standard 4) Streams.
- Tests: `test_prof`, `test_parallel`, `test_resident`, `StarChain.*`, `ThreadPlan.*`,
  `Regression.ThreadCountDoesNotChangeTheBits` (alle Beispiele, 1/2/4 Threads, sequentiell und Baum); CI-Job `tsan`
  (Preset `cpu-tsan`); `bench/` mit den Profilfällen; Leistungsbericht in `docs/benchmarks.md`.

### Geändert
- Einzelne Rechnungen nutzen alle Kerne (vorher nur Scans); Ergebnisse bitgleich zu v0.13 (Golden Files unverändert).
- `SectionedModel`: Schnitt-Cache mit Mutex und `shared_ptr`, Interpolant über `call_once` (threadsicher).
- Modenklassifikation als `modeOrder` ausgelagert (bitgleich); Eigen-Threads nur noch über `threads.blas`.
- Plan v2.7 (Umsetzung M15).

## [0.13.0] – M14 (CAD-Slicer mit Nurbforge)

### Hinzugefügt
- Nurbforge v0.6.0 (eigener PR, gemergt und getaggt): Schnitte von Shells (`nf::section`), kritische Höhen,
  Extrusionsbänder, STEP-Namen (`read_step_named`), Einheit µm, Thread-Sicherheit (TSAN-Job).
- `FG_ENABLE_NURBFORGE` (FetchContent, Tag `FG_NURBFORGE_TAG` = v0.6.0, statisch), Presets `cpu-cad`, `win-cad`.
- Modul `fg::cad` (`fg/cad/cad.hpp`; nur `src/cad` inkludiert Nurbforge): STEP → FG-Koordinaten (Einheit, Achsen ohne
  Spiegelung, Zellursprung), Geometrieprüfung mit festen Codes, Bandzerlegung zwischen kritischen Höhen, je Band
  homogen / `contour2d` (exakter Schnitt) / `SectionedLayer` (Nurbforge-Provider, 3D-Wandnormalen).
- Schichttyp `cad` in JSON (expandiert beim Laden; Pfade relativ zur Konfiguration, `parseConfig(text, baseDir)`),
  `fgsolve --cad-report`, `fg::io::cadReport`, `fg::io::cadPreview`.
- `SectionedLayer::chebyshevNodes`: Q(y) aus Tschebyschow-Interpolation mit Fehlerkontrolle (K gegen 2K).
- Python `fg.cad.check`, `fg.cad.preview` (SVG), `fg.cad.llm_tool_spec`, `fg.cad.available`; `load_config` löst
  `file` von cad-Schichten auf; Wheels mit Slicer.
- `tests/cad/test_cad.cpp` (V16b, V17, V18, zwei Materialien mit Hohlraum, fehlerhafte Dateien, Dateien von Platte),
  `tests/data/cad/`, `scripts/make_cad_examples.py`, Beispiel `examples/cad_pillar_array.json` mit Golden File,
  `docs/cad.md`, CI-Job `cad`.

### Geändert
- `fgsolve --info` meldet den Slicer; Regression und pytest überspringen cad-Beispiele ohne Slicer.
- Plan v2.6 (Präzisierungen M14).

## [0.12.0] – M13 (Gekrümmte Querschnitte, höhenabhängige Schichtmodelle)

### Hinzugefügt
- `geometry::BezierSegment`, `Contour2D` (Polygon, Kreis/Ellipse als rationale Quadratiken, NURBS → Bézier nach
  Piegl–Tiller A5.6), `Region` (äußere Kontur, Löcher); `quadrature::gaussLegendre`.
- `fourier::regionCoefficients`: exakte Koeffizienten aus Konturintegralen (Plan eq:contour, eq:contour0), adaptive
  Gauß-Legendre-Quadratur je Bézierstück, alle (p, q) über zwei Matrixprodukte.
- Algorithmus NF aus Randstichproben (`BoundarySample`, `BoundaryChain`, `normalFieldRaster(chains, …)`):
  Hüllkörperhierarchie, 3D-Normalen, 3×3-Projektorfeld, `rasterCoefficients`; auf Polygonen bitgleich zum bisherigen Pfad.
- Schichttyp `contour2d` (JSON: `circle`, `ellipse`, `polygon_um`, `nurbs`, `holes`; Disjunktheitsprüfung),
  Beispiel `examples/v15_circular_pillars.json` mit Golden File.
- `SectionedLayer` + `SectionProvider`/`CrossSection` (`fg/geometry/section.hpp`): Querschnitte, die sich mit y
  ändern; Modell je Höhe mit Cache. Grundlage für M14.
- Tests `Contour.*`, `NormalFieldSamples.*`, V15 (`test_v15_contour.cpp`), V16a (`test_v16a_sectioned.cpp`).

### Geändert
- `detail::LayerModel::normalAt(x, y, z)`, `projectorAt(x, y, z)`, `NN(y)` (bisher ohne y); alle Aufrufer übergeben
  die Höhe. `Polygon2DModel` → `CrossSectionModel` mit Geometrie-Schnittstelle (Polygon/Kontur), bitgleich.
- Gauß-Glättung mit vorberechneten periodischen Indizes (bitgleich, schneller).
- `tests/regression/examples.txt` löst bei Änderung den Neubau von `test_regression` aus.
- Plan v2.5 (Präzisierungen M13: BVH, Ketten, V15 zweistufige Richardson, V16a Laurent).

## [0.11.0] – M12 (Python-Anbindung)

### Hinzugefügt
- `fg::Simulation` (`fg/solvers/simulation.hpp`): Konfiguration einmal parsen, Schichtmodelle über viele Rechnungen
  behalten; `setParameter(pointer, value)` baut nur geänderte Schichten neu (bitgleich zum Neuaufbau). Intern
  `detail::solveWithModels`, `detail::buildLayerModel`, `detail::applyThreadSettings`.
- Python-Paket `fouriergrating` (pybind11, scikit-build-core, `pyproject.toml`): `load_config`, `solve`, `solve_scan`
  (`ScanResult`, `to_pandas`), `near_field`, `Simulation`, `set_backend`, `cuda_available`, `core_variant`,
  `FourierGratingError`; GIL während jeder Rechnung freigegeben; `py.typed`, `_core.pyi`.
- Distribution `fouriergrating-cuda12` (`python/cuda12/`), Kernwahl beim Import (`FG_CORE=cpu` erzwingt CPU),
  CUDA-DLL-Verzeichnisse über `os.add_dll_directory`.
- CMake: `FG_BUILD_PYTHON`, `FG_PYTHON_PACKAGE`; unter MSVC ignoriert `FGDependencies.cmake` MSYS2-Präfixe und die
  Paketregistry (statt im Preset).
- `fg::derivedFields`, `fg::nearFieldMeta` (gemeinsam für `NearFieldWriter` und Python).
- pytest `tests/python`, Beispiel `examples/python/scan_lamellar.py`, `docs/python.md`.
- CI: Jobs `regression`, `windows-msvc`, `python`; `.github/workflows/wheels.yml` (cibuildwheel: manylinux_2_28 und
  Windows für CPython 3.10–3.13, CUDA-Wheel für Windows).

### Geändert
- Plan v2.4: Toleranz der Golden Files über Plattformen 1e-11 absolut (statt 1e-13); gemessen Linux-GCC 13 gegen
  MinGW-GCC 16 bis 1,3e-12 (RK Gold-Echelette, V10), MSVC gegen MinGW ≤ 1,2e-14. Gleiche Maschine: bitgleich.
- Eigen per `FetchContent` nur als Header (Eigens CMake scheiterte auf GitHub-Windows-Runnern an gfortran).
- `wheels.yml`: `Jimver/cuda-toolkit@v0.2.36` (CUDA 12.6.3).

### Abnahme
- CI und Wheels grün; CI-Wheels auf dem Windows-Rechner mit RTX 3090: CUDA-Kern geladen (DLLs aus `nvidia-*-cu12`),
  V10 GPU = CPU auf 9e-15 (5,3 s statt 16,9 s), pytest 27/27.

## [Unreleased] (vor M12)

### Plan
- Plan v2.2: Optionen `split_polarizations` (Spezialfälle) und `fft_sampling` (Allgemeine Funktionen).
- Plan v2.0: Erweiterungsplan Version 2 (Teile 5–6, Meilensteine M12–M19): Python-Anbindung, Konturkoeffizienten und
  höhenabhängige Schichtmodelle, CAD-Slicer mit Nurbforge, Profiling/Parallelisierung, wählbare Genauigkeit
  (FP64/gemischt/FP32), Koordinatentransformationen und ASR, PML für isolierte Strukturen, diskrete Adjungierte.
  Roadmap, `AGENTS.md` und `CLAUDE.md` entsprechend ergänzt.
- Plan v2.1: Entscheidungen des Projektinhabers eingearbeitet (isolierte Strukturen mit Streufeldformulierung für
  schrägen Einfall, Windows-Wheels mit MSVC und CUDA-Variante, Geometrieprüfung für LLM-erzeugte Nurbforge-Geometrie,
  Einheit µm, FP32-Standard auf der RTX 3090).

### Hinzugefügt
- MSVC-Port (M12): Presets `win-msvc` und `win-msvc-cuda` (`86;89`), `scripts/msvc.cmd` (vcvars64), MSVC-Flags
  `/W4 /permissive- /utf-8 /bigobj`, externe Header `/external:W0`; MSYS2-Präfixe und Paketregistry im Preset
  ausgeschlossen. Alle Tests grün, `Backend.*` mit CUDA 13.4 auf der RTX 3090 grün. `docs/build_windows.md` Abschnitt 4.
- Golden Files `tests/regression/golden/` und Test `Regression.*` (Label `regression`): Effizienzen aller Beispiele mit
  Standardoptionen; 1e-13 bezogen auf die einfallende Leistung (Plan v2.3; rein relativ scheitern sehr kleine
  Ordnungen beim Compilerwechsel), `FG_REGRESSION_EXACT=1` bitgleich, `FG_REGRESSION_UPDATE=1` erzeugt neu.
- Option `"split_polarizations"` (`SolverOptions::splitPolarizations`, Standard aus): 1D in-plane isotrop
  (γ = 0, isotropes Substrat und isotrope Schichten) rechnet `solve()` TE und TM als getrennte 2n_f-Systeme
  (Plan eq:TE, eq:TM) – Grenzflächen, Schicht-S-Matrizen, RK-Transfer und Sternprodukte mit halber Größe;
  Einbettung blockdiagonal, nicht angeregte Polarisation (reines TE/TM) übersprungen. RCWA Faktor 1,9–5.
  Neu: `Polarization`, `polarizationBlocks`, `polarizationBasis`, `buildPropagationMatrix(Q, alpha, Polarization)`,
  `layerModesPolarization`.
- Option `"fft_sampling"` (`SolverOptions::fftSampling`, Standard aus): `fourier::fft` (beliebige Länge, Radix 2 /
  Bluestein) für `sampledCoefficients1D/2D` und `gridCoefficients2D`; Schalter `fourier::setFFTSampling` (thread-lokal).
- Tests `PolarizationSplit.*`, `Sampling.FFTEqualsDirectDFT`, `Sampling.FFTSwitchMatchesDFT`,
  `Sampling.SampledCoefficientsExactForTrigPolynomials`.
- Mit Standardoptionen sind alle `examples/` bitgleich zum Stand vor der Änderung (CSV-Vergleich).
- Umgebungsvariablen `FG_BACKEND` / `FG_BACKEND_THRESHOLD` überschreiben die Backend-Wahl (ganze Testsuite auf
  der GPU: `FG_BACKEND=cuda FG_BACKEND_THRESHOLD=0 ctest`); Test `Backend.EnvironmentOverride` (portabel, auch MinGW).
- V10 bei N = 11: Peak ≈ 1,3898 µm, ≈ 9,2 %.

### Geändert
- `layerModesAuto`: Klassifikation der TE/TM-Blöcke in eine gemeinsame Funktion ausgelagert (bitgleich).

## [0.10.0] – M11 (Scans, Sweeps, Anwenderdoku)

### Validierung (nachgetragen)
- V10 bei N = 9 und 10: Transmissionspeak bei 1,3876 / 1,3889 µm (Buch: 1,388 µm), Höhe 4,2 / 6,2 %
  (Buch 7,5 %); Tensor- und Vektorvariante von Algorithmus NF gleichwertig (N = 7, 8).

### Hinzugefügt
- Scans über λ, φ, ψ und δ (Zahl, Liste oder Bereich; kartesisches Produkt), parallel über `solveScan`.
- Parameter-Sweeps über beliebige Zahlenwerte der Konfiguration (`"sweep": {"parameter": <JSON-Pointer>,
  "values": ...}`), `fg::io::parseSweep`/`loadSweep`, CSV mit `sweep_parameter`, `sweep_value`.
- `scripts/plot_efficiencies.py` (Spektren und Sweeps aus `efficiencies.csv`).
- Beispiel `examples/sweep_lamellar_depth.json` (Smoke-Test), Tests `Scan.*`.
- Anwenderhandbuch `docs/user_guide.md` mit Laufzeittabelle; README-Status aktualisiert.

### Geändert
- `fgsolve` rechnet Sweeps (ein Block je Wert) und Scans parallel.
- Plan v1.9.

## [0.9.0] – M10 (CUDA-Backend, ohne GPU-Abnahme)

### Hinzugefügt
- `CudaBackend` (`src/cuda/cuda_backend.cpp`): cuBLAS `Zgemm`, cuSOLVER `Xgetrf`/`Xgetrs` (64 Bit) und
  `Xgeev` über die Host-API; RAII für Device- und Pinned-Speicher, Stream je Backend. Kein `nvcc` nötig.
- Dispatch (`linalg::setPreference`, `offload`, `device`): RCWA-Eigenprobleme, stabile Schicht-S-Matrix,
  Sternprodukt, T→S und RK-Schrittpropagator laufen ab der Schwelle auf der GPU; CPU-Pfad unverändert.
- `SolverOptions::backend` / `backendThreshold`, JSON `"backend"`, `"backend_threshold"`.
- `solveScan`: paralleler Wellenlängen-/Winkelscan (OpenMP, ein Löser und Stream je Thread); `fgsolve` nutzt ihn.
- `fgbench` (Benchmark-Tabellen), `docs/benchmarks.md` (Sandbox-Referenz, Zielsystem offen).
- `fgsolve --info` meldet „CUDA gebaut, kein Gerät“ statt abzubrechen.
- Tests `Backend.*` (GPU-Vergleiche werden ohne Gerät übersprungen).

### Geändert
- CMake: CUDA ohne `enable_language(CUDA)` (nur `find_package(CUDAToolkit 12.6)`), `cuda_backend.cu` → `.cpp`.
- Plan v1.8.

## [0.8.0] – M9 (Nahfeld 2D)

### Hinzugefügt
- `Solver::nearField` für 2D-Gitter: horizontale Karte u(x, y0, z) (`nz` Punkte je Periode, Achse `z`) und
  vertikaler Schnitt u(x, y, z0); Synthese mit separierten Phasen.
- FFF-konsistente Rekonstruktion in 2D: `LayerModel::projectorAt(x, z)`, `rank1Projector()`;
  isotrope Polygone mit dem Tensorfeld (E = e_t + d_n/ε), Variante B / anisotrop mit dem Rang-1-Feld,
  Profile2D mit eq:N2d; Punktwerte aus dem NF-Raster (bilinear, vorzeichenfrei).
- `epsAt(x, y, z)` und `normalAt(x, z)` für alle Schichtmodelle (Materialkarten der 2D-Schnitte).
- Validierung V12 (2D): z-invariant = 1D, Stetigkeit, Flusskonstanz, RK = RCWA, D_x-Sprung an Säulenwänden
  (isotrop und anisotrop), Schnitte und Dateien; Beispiel `examples/nearfield_pillars_2d.json` (Smoke-Test).

### Geändert
- Plan v1.7.

## [0.7.0] – M8 (2D anisotrop)

### Hinzugefügt
- Anisotrope Materialien in `Polygon2D` (Laurent exakt; `FFF_B`, Standard) und `Profile2D`
  (`FFF_A`, Standard; `FFF_B`; Laurent); anisotrope 1D-Profile in 2D-Stapeln.
- Variante B für Querschnitte: Rang-1-Normalenfeld aus Algorithmus NF (normierter geglätteter Gradient,
  Hauptachse des Tensors als Rückfall), Koeffizienten als glatter Anteil (FFT) plus Faltung der
  Materialdifferenz mit den exakten Polygonkoeffizienten (`fourier::convolveCoefficients`).
- `FFF_B` auch für isotrope `Polygon2D`, `FFF_A`/`FFF_B` auch für isotrope `Profile2D` (Vergleiche).
- Radix-2-FFT in `gridCoefficients2D` (komplexe Variante); NF-Raster auf Zweierpotenzen.
- Validierung V8 (2D), Spiegeltest mit anisotropem Tensor, Konvergenzstudie Co-Säulen und
  anisotroper Eierkarton.

### Geändert
- Anisotrope Faktorisierungen intern vereinheitlicht (`anisoEval`, `anisoQFromGrids`), 1D unverändert.
- Plan v1.6.

## [0.6.0] – M7 (2D isotrop)

### Hinzugefügt
- `Polygon2D`-Schichten: exakte Polygon-Koeffizienten von ε und 1/ε (Polygone dürfen den Zellrand kreuzen,
  Orientierung wird normalisiert); Faktorisierung `FFF_sym` (Standard), `FFF_iso`, Laurent.
- Algorithmus NF (Plan alg:nfield, v1.5): exakte Abstände zu den Materialgrenzen (periodische Bilder, innere
  Kanten entfernt, Mittelachsen gemittelt), periodische Gaußglättung, DFT; Tensorvariante (Standard) und
  Vektorvariante (v1.4) wählbar (`normal_field.method`).
- `Profile2D`-Schichten (`egg_crate`, `sinus_x`, Funktion): ε(y) zeilenweise exakt in x, in z Trapezregel bzw.
  Gauß–Legendre zwischen kritischen Linien; Normalenfeld nach eq:N2d (vorzeichenbehaftet, `FFF_rig` Standard).
  RK (M(y)) und RCWA-Treppe.
- Isotrope `profile1d`-Schichten in 2D-Stapeln (z-invariant).
- JSON: `polygon2d` (`vertices_um`, `rectangle_um`, `normal_field`), `profile2d`; Prüfung `period_um.z > 0` in 2D.
- `fourier::gridCoefficients2D` (separable DFT eines Rasters).
- Validierung V9 (z-invariante Strukturen = 1D konisch), V10 (Ag-Lochgitter, qualitativ), Spiegeltest x↔z,
  Eierkarton (Energie, RK = RCWA); Unit-Tests Algorithmus NF und 2D-Schichtmodelle.

### Geändert
- Plan v1.5: Algorithmus NF mit Tensorglättung (Befund: die Vektorglättung fällt an Mittelachsen auf ½ zurück
  und reproduziert Lamellen nicht exakt); Profile2D ohne Marching Squares.
- RK-Spektralradius für die Scheibenzerlegung: oberhalb n = 800 Potenziteration statt voller Eigenzerlegung.

### Befund
- V10 (Ag, ε = −79,9 + 1,79i): `FFF_iso` liefert R > 1 (N = 3, 6); der Transmissionspeak oberhalb des
  SPP-Pols wandert mit N nach Rot (1,383 → 1,386 → 1,391 → 1,404 µm für N = 5…8), Höhe 3–10 %;
  quantitativer Vergleich mit dem Buch (7,5 %) offen.

## [0.5.0] – M6 (Nahfeld 1D)

### Hinzugefügt
- Gemeinsame Schichtmodelle (`detail::LayerModel`) für Lösung und Feldrekonstruktion.
- Aufzeichnende Assemblierung: interne Amplituden (S_L/S_U), modale Scheiben, RK-Transfermatrizen
  an den angeforderten Höhen.
- `Solver::fieldCoefficients` und `Solver::nearField`: horizontaler und vertikaler Schnitt,
  direkte oder FFF-konsistente Rekonstruktion (Projektorform, iso und aniso).
- NPY/JSON-Ausgabe (`NearFieldWriter`), `fgsolve` schreibt Schnitte bei `nearfield.enabled`.
- Validierung V12 (Airy analytisch, Stetigkeit, Flusskonstanz, RK = RCWA, D_x-Stetigkeit).

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
