# Arbeitsregeln für KI-Agenten (und Menschen)

Diese Datei gilt für jeden, der an FourierGrating arbeitet. Lies sie vollständig, bevor du Code änderst.

## 1. Quelle der Wahrheit

1. **`docs/implementation_plan/main.tex`** (PDF daneben) ist die verbindliche Spezifikation.
   Jede implementierte Formel muss dort stehen; im Code wird sie per Label zitiert,
   z. B. `// Plan eq:Mmu1b`.
2. Bei Widerspruch zwischen Plan und Buch (Nevière/Popov) gilt der Plan. Die bekannten
   Buchfehler stehen in `docs/book_errata.md` und im Plan, Abschnitt „Buchfehler“.
3. Findest du einen Fehler im Plan: erst einen Test schreiben, der ihn belegt, dann Plan **und**
   Code im selben Commit korrigieren und den Eintrag in `CHANGELOG.md` ergänzen.

## 2. Reihenfolge

* Arbeite die Meilensteine in `docs/roadmap.md` **der Reihe nach** ab (M0 → M11, dann Version 2: M12 → M19,
  Plan Teile 5–6). Die offene GPU-Abnahme von M10 blockiert nur die GPU-Teile von M15 und M16.
  Ein Meilenstein ist erst fertig, wenn alle seine Abnahmetests grün sind.
* Status in `docs/roadmap.md` nach jedem Meilenstein aktualisieren (`[x]`).
* Keine Funktionen „auf Vorrat“ implementieren, die der aktuelle Meilenstein nicht braucht.

## 3. Tests zuerst

* Jede neue Formel bekommt einen Unit-Test (`tests/unit/test_<modul>.cpp`), jede physikalische
  Funktion mindestens einen Referenzfall aus `docs/validation.md` (`tests/validation/`).
* Bewährte Prüfmuster: Vergleich mit direkter Elimination (Zufallsmatrizen), analytische
  Grenzfälle (Fresnel, Airy), Energieerhaltung bei verlustfreien Medien, Reziprozität,
  RK ↔ RCWA auf Lamellen, 2D-Code mit z-invarianter Struktur ↔ 1D-konisch,
  anisotroper Code mit isotropem Tensor ↔ isotroper Code.
* `ctest --preset cpu-release` muss vor jedem Commit grün sein. Validierungstests tragen das
  CTest-Label `validation` und dürfen länger laufen.

## 4. Numerik-Regeln

* Kern ausschließlich in `std::complex<double>`; normierte Einheiten (`k0 = 1`, `H̃ = Z0 H`).
  Einfache bzw. gemischte Genauigkeit nur innerhalb des Linear-Algebra-Backends (ab M16, Plan „Wählbare Genauigkeit“).
  Einheiten (µm, Grad) nur an der IO-Grenze umrechnen.
* Keine expliziten Inversen, wenn ein LU-Solve genügt.
* Wurzelzweig für β, Moden-Klassifikation und Rayleigh-Anomalie exakt wie im Plan.
* Konvention der S-Matrix: Plan eq:Sdef. Verkettung nur über `fg::star`.
* **Eigen-Falle:** `a.dot(b)` konjugiert bei komplexen Vektoren das erste Argument, und `a.cross(b)`
  liefert für komplexe Vektoren das **konjugierte** Kreuzprodukt. Für bilineare Ausdrücke
  (ξ = Nᵀ ε N, v × T, E·s) immer `(a.array() * b.array()).sum()` bzw. ein eigenes Kreuzprodukt
  verwenden. (Dieser Fehler hat in M5 die anisotrope FFF still verfälscht.)
* Unstetige Funktionen nie per FFT abtasten, wenn eine exakte Formel existiert
  (Intervall, Rechteck, Polygon).

## 5. Code-Stil

* C++20, `clang-format` (Datei im Repo), Warnungen aus (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion`).
* Namensraum `fg`, Unterräume nach Verzeichnis (`fg::fourier`, `fg::linalg`, …).
* Kein CUDA-Code außerhalb von `src/cuda/` und dem Backend; der Physikcode spricht nur mit
  `fg::linalg::Backend`.
* Öffentliche Funktionen bekommen einen Kurzkommentar mit Formelverweis.
* Fehler über `fg::Error` / `FG_REQUIRE`; Unfertiges über `fg::NotImplemented("… Mx")`.

## 6. Git

* Kleine, thematisch abgeschlossene Commits; Präfix mit Meilenstein, z. B. `M2: RCWA lamellar TE/TM`.
* **Niemals Zugangsdaten (Tokens, Passwörter) ins Repository schreiben**, auch nicht in Remote-URLs
  in Skripten.
* Erzeugte Dateien (Build, Ausgaben) gehören nicht ins Repo (`.gitignore`). Ausnahme: das PDF des
  Plans (`docs/implementation_plan/main.pdf`) wird bei Änderungen des Plans neu erzeugt und eingecheckt.

## 7. Version 2: Abhängigkeiten und Regressionsschutz

* Mit ausgeschalteten neuen Optionen bleiben alle Ergebnisse bitgleich zu v0.10 (Golden Files in
  `tests/regression/golden/`, Test `Regression.*`). Golden Files nur mit Begründung im CHANGELOG neu erzeugen.
* Nurbforge (`PhotonPower/Nurbforge`) wird nur über `FG_ENABLE_NURBFORGE` und ein festes Tag eingebunden; nur
  `src/cad/` inkludiert Nurbforge-Header. Fehlende CAD-Funktionen werden in Nurbforge selbst ergänzt (eigener PR mit
  Tests, neues Tag), nicht in FourierGrating nachgebaut.
* Python-Bindungen (`python/`) enthalten keine Physik, nur Umwandlung und Aufrufe.
* Gradienten (M19) werden nur über die dokumentierten Rückwärtsfunktionen in `src/adjoint/` gebildet; jede bekommt
  einen Test gegen Differenzenquotienten.

## 8. Zielsystem

AlmaLinux 10 (x86_64, GCC 14), CMake ≥ 3.24, CUDA ≥ 12.6 Update 2 (für `cusolverDnXgeev`),
GPU NVIDIA RTX 6000 Ada (sm_89, 48 GB); zusätzlich Windows mit RTX 3090 (sm_86, 24 GB) für M15/M16. Die FP64-Leistung dieser GPU ist begrenzt; der CPU-Pfad
bleibt Referenz und muss immer funktionieren (CI baut nur CPU).
