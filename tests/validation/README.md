# Validierungstests

Ab M2 entstehen hier die Referenzfälle V1–V12 aus `docs/validation.md`, je als GoogleTest
(`test_v<k>_*.cpp`) mit CTest-Label `validation`. Jeder Test

1. lädt die passende Konfiguration aus `examples/` (Feld `reference`),
2. rechnet mit RK **und** RCWA (sofern sinnvoll),
3. prüft Sollwert und Toleranz,
4. schreibt eine Konvergenztabelle (N, Scheiben, Ergebnis, Laufzeit) nach
   `${CMAKE_BINARY_DIR}/validation/<id>.csv`.
