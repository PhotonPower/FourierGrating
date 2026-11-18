# CAD-Slicer: STEP → Schichten (M14)

Eine STEP-Datei beschreibt den Inhalt einer Einheitszelle (B-rep-Körper). FourierGrating zerlegt sie entlang
der Stapelrichtung y in **Bänder** zwischen den kritischen Höhen und erzeugt daraus homogene, y-invariante
(`contour2d`) und y-variable (`SectionedLayer`) Schichten. Gelesen und geschnitten wird mit
[Nurbforge](https://github.com/PhotonPower/Nurbforge) (MIT), eingebaut v0.12.1. STEP-Dateien können mit jeder
Nurbforge-Version ab 0.6.0 (oder einem anderen CAD-Programm) geschrieben werden.

```text
STEP --nf::read_step_named--> Körper (+ Namen) --Einheit/Achsen--> FG-Koordinaten (µm)
  --nf::critical_heights--> Bänder [y_k, y_k+1]
  --je Band: homogen | Extrusion (ein Schnitt, exakt) | y-variabel (Schnitte bei Bedarf)
  --nf::section--> Schleifen + Flächennormalen --> Region/Contour2D + Randstichproben (M13)
```

## Bauen

| Weg | Befehl |
|---|---|
| CMake (Linux) | `cmake --preset cpu-cad && cmake --build --preset cpu-cad` |
| CMake (Windows, MinGW) | `cmake --preset win-cad` |
| Python | `pip install fouriergrating` (Wheels mit Slicer) bzw. `pip install .` |

`FG_ENABLE_NURBFORGE=ON` lädt Nurbforge (Tag `FG_NURBFORGE_TAG`, Standard `v0.12.1`) per `FetchContent` und linkt
es statisch; für eine lokale Kopie `-DFETCHCONTENT_SOURCE_DIR_NURBFORGE=<Pfad>`. Ohne die Option meldet eine
`cad`-Schicht `NotImplemented`. `fgsolve --info` zeigt, ob der Slicer eingebaut ist. Nur `src/cad/` inkludiert
Nurbforge-Header.

## STEP-Dateien erzeugen

Quelle sind Modelle aus Nurbforge (z. B. vom Sprachmodell über die Ollama-Werkzeuge): in **Mikrometern**
konstruieren, mit `unit="um"` schreiben und **einen benannten Körper je Material** ausgeben.

```python
import nurbforge as nf
pillar = ...  # nf.Shell, Koordinaten in µm, CAD +Z = Stapelrichtung
nf.write_step("pillars.step", [pillar, lid], unit="um", names=["Pillar", "Lid"])
```

`scripts/make_cad_examples.py` erzeugt so die Beispiele (`examples/cad/`) und Testdaten (`tests/data/cad/`).

## Konfiguration

```json
{"type": "cad", "file": "cad/pillar_array.step", "unit": "auto", "axes": {"y": "+Z", "x": "+X"},
 "y_range_um": [0.0, 0.26], "cell_origin_um": [0.0, -0.5], "background": "superstrate",
 "bodies": [{"name": "Pillar", "material": "si"}, {"index": 1, "material": "sio2"}],
 "overlap": "error", "periodic_union": false,
 "slicing": {"critical_tol_um": 1e-6, "min_band_um": 1e-5, "rcwa_slices": 40,
             "y_interpolation": "none", "chebyshev_nodes": 16},
 "factorization": "auto", "normal_field": {"method": "tensor", "grid": 512, "sigma": 3}}
```

| Schlüssel | Bedeutung |
|---|---|
| `file` | STEP-Datei, relativ zur Konfigurationsdatei (Python: `load_config` löst auf); alternativ `step` mit dem Dateiinhalt |
| `unit` | `auto`: wie in der Datei deklariert; `mm` / `um`: Koordinaten sind in dieser Einheit (korrigiert falsche Deklarationen) |
| `axes` | CAD-Achsen für FG y (Stapelrichtung) und x; FG z = x × y (Standard: z = −CAD Y). Ein angegebenes `z` muss dazu passen: gespiegelte (linkshändige) Achsen sind verboten, sie verfälschen chirale Strukturen |
| `y_range_um` | y-Bereich (Standard: Hüllquader der Körper); Körper darüber hinaus → `OUT_OF_RANGE` |
| `cell_origin_um` | (x₀, z₀) der Zelle in FG-Koordinaten (Standard: Minimum des Hüllquaders); Körper dürfen über die Zelle ragen (periodische Bilder) |
| `background` | Material außerhalb der Körper: Name, `"superstrate"` oder `"substrate"` |
| `bodies` | Materialzuordnung je Körper über den STEP-Namen (`name`) oder die Position (`index`); nicht zugeordnete Körper sind ein Fehler |
| `overlap` | `error` (Standard) oder `priority`: frühere Körper gewinnen (Boolesche Differenz) |
| `periodic_union` | nur 1D-Gitter: Körper mit ihren überlappenden periodischen Bildern vereinigen. In 2D wirkungslos: Überlappt ein Körper sein Bild, umläuft die Vereinigung die Zelle (`PERIODIC_OVERLAP`) |
| `slicing.critical_tol_um` | kritische Höhen näher als dies werden zusammengelegt |
| `slicing.min_band_um` | dünnere Bänder werden dem Nachbarband zugeschlagen (Warnung `FEATURE_SMALL`) |
| `slicing.rcwa_slices` | RCWA-Treppenstufen je y-variablem Band |
| `slicing.y_interpolation` | `chebyshev`: Q(y) je Band aus `chebyshev_nodes` Tschebyschow-Knoten (verdoppelt bis K = 2K auf 1e-10, höchstens 256; sonst exakte Schnitte mit Warnung) |
| `factorization`, `normal_field` | wie bei `contour2d` (Standard FFF_sym / FFF_B) |

In einem 2D-Gitter berührt eine z-invariante Linie, die über genau eine z-Periode extrudiert ist, ihr periodisches
Bild (`PERIODIC_OVERLAP`, auch mit `periodic_union`): Liniengitter importiert man als **1D-Gitter** (nächster
Abschnitt). Hüllquader (Standard für `y_range_um` und `cell_origin_um`) sind die der getrimmten Körper, nicht die der
Trägerflächen.

Der `cad`-Eintrag wird beim Laden in seine Bänder expandiert (von unten nach oben); Bandgrenzen sind
Schichtgrenzen. Es wird nie auf einer kritischen Höhe geschnitten.

## 1D-Gitter (M23)

Mit `"dimension": "1D"` liest der Slicer z-invariante Körper, also Linien eines Liniengitters, die entlang FG z
(standardmäßig −CAD Y) extrudiert sind. Die Extrusionslänge ist beliebig. Jeder Körper wird **einmal** in der Ebene
`z = z_cut_um` geschnitten; das liefert sein Profil in der (x, y)-Ebene. Ab da rechnet FourierGrating ohne Nurbforge:
Je Höhe schneidet es die Profilkurven (exakte rationale Bézierstücke, also auch Kreisbögen) mit der Waagerechten und
erhält Intervalle mit Materialien und die Normalen an den Grenzen.

```json
{"dimension": "1D", "period_um": {"x": 0.3}, "truncation": {"Nx": 12},
 "layers": [{"type": "cad", "file": "cad/line_grating.step", "unit": "um", "cell_origin_um": 0.0,
             "z_cut_um": -0.15, "background": "superstrate", "bodies": [{"name": "Line", "material": "resist"}]}]}
```

| Schlüssel | Bedeutung in 1D |
|---|---|
| `cell_origin_um` | nur x₀: Zahl oder `[x0]` (Standard: Minimum des Hüllquaders in x) |
| `z_cut_um` | FG z der Schnittebene (Standard: Mitte des gemeinsamen z-Bereichs aller Körper); nur in 1D erlaubt |
| `normal_field` | nur `grid` für `FEATURE_SMALL`; das Normalenfeld ist analytisch (Plan eq:cad1d_N) |

* Bänder: `homogeneous` (ein Material über die Periode), `lamellar1d` (alle Wände senkrecht, exakt mit der Li-Regel)
  oder `sectioned1d` (y-variabel; RK an den Stufenknoten, RCWA-Treppe mit `rcwa_slices`, optional Tschebyschow).
  Der Bericht gibt Längenanteile statt Flächenanteilen, dazu `"dimension": "1D"` und `z_cut_um`.
* Periodizität: Ein Körper, der die ganze Periode füllt (Restschicht) oder über den Zellrand ragt, wird an der Naht
  ohne Grenzfläche verschmolzen. Überlappt er sein Bild, folgt `PERIODIC_OVERLAP`; `periodic_union` vereinigt.
  `overlap: "priority"` arbeitet mit Intervallen.
* Faktorisierung: Standard `FFF_sym` (isotrop) bzw. `FFF_B` (anisotrop); erlaubt sind `laurent`, `FFF_iso`,
  `FFF_sym` bzw. `laurent`, `FFF_B`. Das Normalenfeld interpoliert den Winkel zwischen den Grenzen linear.
* Gradienten: Materialwerte wie gewohnt; Formparameter über `fg::StackParameter` (C++), auch für `sectioned1d`.
* Nicht mit PML (`aperiodic`).
* Beispiel: `examples/cad_line_grating_1d.json` (Linie mit CD 100 nm auf halber Höhe, h = 100 nm, 15°-Wänden und
  Ecken mit r = 10 nm, Kopf konvex, Fuß konkav), STEP aus `scripts/make_cad_examples.py`.

## Bandzerlegung

1. Kritische Höhen aller Körper im y-Bereich (`nf::critical_heights`: Ecken, ebene Flächen senkrecht zu y,
   Höhenextrema auf Flächen und Kanten), zusammengelegt innerhalb `critical_tol_um`.
2. Je Band: kein Rand → `homogeneous` (Material am Zellmittelpunkt der Bandmitte); alle berührenden Flächen
   parallel zu y (`nf::is_extrusion_band`) → `extrusion`, ein exakter Schnitt in Bandmitte als `contour2d`;
   sonst `sectioned` (`SectionedLayer` mit Schnitten bei Bedarf, Cache je Höhe; Normalen der Wände in 3D).
3. Schnitte: Extrusionen und Rotationen um y (Zylinder, Kegel, Kugel, Torus) liefern exakte Isoparameterkurven
   (rationale Bézier-Stücke → exakte Koeffizienten, M13); andere Flächen werden auf 1e-9 der Modellgröße verfolgt.

`fgsolve --cad-report config.json` gibt Zerlegung und Befunde als JSON aus (Exit-Code 0 bei `ok`, sonst 4).

## Geometrieprüfung (für die Prüfschleife eines Sprachmodells)

| Code | Prüfung |
|---|---|
| `SHELL_OPEN` | Körper nicht geschlossen (freie Kanten), mit Position |
| `SHELL_ORIENTATION` | nicht konsistent orientiert oder negatives Volumen (korrigiert, Warnung) |
| `BODY_OVERLAP` | zwei Körper überlappen (Schnittprüfung an drei Höhen je Band; berührende Körper sind erlaubt); bei `priority` auch, wenn die Boolesche Differenz in Nurbforge scheitert (offene Ergebnishülle) |
| `PERIODIC_OVERLAP` | ein Körper überlappt sein periodisches Bild (bzw. das eines anderen) oder berührt sein eigenes entlang eines Randstücks (Körper umspannt die Periode) |
| `OUT_OF_RANGE` | Geometrie außerhalb des y-Bereichs; 1D: `z_cut_um` schneidet nicht alle Körper |
| `UNIT_SUSPECT` | Abmessungen passen nicht zu Periode / Wellenlänge (Faktor 10³ daneben) |
| `FEATURE_SMALL` | Merkmal kleiner als zwei Rasterabstände des Normalenfeldes bzw. dünnes Band (Warnung) |
| `BODY_UNASSIGNED` | Körper ohne Material bzw. Zuordnung ohne Körper |
| `STEP_UNSUPPORTED` | nicht lesbare Datei oder nicht unterstützte Entität (Name und Instanz) |
| `NOT_Z_INVARIANT` | 1D-Gitter: Körper ist nicht entlang z extrudiert |

Positionen in µm, FG-Koordinaten. Python:

```python
import fouriergrating as fg
report = fg.cad.check(step_text, cfg)   # {"ok", "errors": [{"code", "message", "where_um"}], "warnings", "bands"}
svg = fg.cad.preview(step_text, cfg, y_um=0.1)   # 1D: fg.cad.preview(step_text, cfg) zeigt das Profil
tool = fg.cad.llm_tool_spec()           # Werkzeug "fouriergrating_check_cad(step, config)" als JSON-Schema
```

`llm_tool_spec()` lässt sich in der Ollama-Prüfschleife von Nurbforge als Werkzeug registrieren; die Kopplung
bleibt einseitig (Nurbforge kennt FourierGrating nicht). C++: `fg::cad::check`, `fg::cad::slice`,
`fg::cad::previewSvg`, `fg::cad::llmToolSpec` (`fg/cad/cad.hpp`), `fg::io::cadReport`, `fg::io::cadPreview`.

## Validierung

`tests/cad/test_cad.cpp` (Label `validation`, nur mit `FG_ENABLE_NURBFORGE`); STEP-Dateien werden zur Laufzeit
mit Nurbforge geschrieben:

* V16b: extrudiertes L-Profil → ein Band `extrusion`; Koeffizienten = `polygon2d` (1e-13), Effizienzen RK und
  RCWA = `polygon2d` (1e-12); gespiegelte Achsen werden abgelehnt.
* V17: Zylinder = exakter Kreis, Koeffizienten = Plan eq:circle (1e-12).
* V18: Kegelstumpf und Kugelkappe (y-variabel) gegen `profile2d` mit g(x, z) (RK, tol 1e-9): 1e-6; Kegelstumpf
  gegen den analytischen Provider aus M13: 1e-9; RCWA-Treppe konvergiert; Energie 1e-10; Tschebyschow = exakt (1e-9).
* Zwei Materialien (Säule mit Deckel) und ein Block mit kugelförmigem Hohlraum: Bandfolge, Energie, RK ≈ RCWA.
* Fehlerhafte Dateien liefern genau ihren Code, korrekte `ok`.
* V24 (1D-Gitter, M23): siehe `docs/validation.md`; ohne Nurbforge `tests/unit/test_line_section.cpp` und
  `tests/unit/test_sectioned1d.cpp`.
