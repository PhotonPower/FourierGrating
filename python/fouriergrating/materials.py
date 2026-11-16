"""Material models from other sources (Plan section "Import aus refractiveindex.info", M22).

Format conversion only (AGENTS.md section 7): the result is a material object of the configuration
(``docs/config_format.md``, "Materialwerte"), evaluated by the core like a handwritten one.

    import fouriergrating as fg
    cfg["materials"]["au"] = fg.materials.from_refractiveindex("data/main/Au/Johnson.yml")
"""

from __future__ import annotations

import os
from typing import Any, Union

__all__ = ["from_refractiveindex"]


def _rows(text: str) -> list:
    rows = []
    for line in str(text).splitlines():
        parts = line.replace(",", " ").split()
        if parts:
            rows.append([float(p) for p in parts])
    return rows


def _interp(x: float, xs: list, ys: list) -> float:
    """Linear interpolation inside [xs[0], xs[-1]] (xs increasing)."""
    if x == xs[-1]:
        return float(ys[-1])
    for i in range(len(xs) - 1):
        if xs[i] <= x < xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return float(ys[i] + t * (ys[i + 1] - ys[i]))
    raise ValueError(f"{x} outside [{xs[0]}, {xs[-1]}]")


def _sellmeier(entry: dict, formula: int) -> dict:
    c = [float(v) for v in str(entry["coefficients"]).split()]
    if len(c) < 3 or len(c) % 2 == 0:
        raise ValueError(f"formula {formula}: coefficients C1, then pairs (B_k, C_k)")
    b = c[1::2]
    cc = c[2::2]
    model: dict[str, Any] = {"A": 1.0 + c[0], "B": b, "C_um2": [v * v for v in cc] if formula == 1 else cc}
    if "wavelength_range" in entry:
        lo, hi = (float(v) for v in str(entry["wavelength_range"]).split())
        model["range_um"] = [lo, hi]
    return {"n": {"sellmeier": model}}


def from_refractiveindex(path_or_text: Union[str, "os.PathLike[str]"]) -> dict:
    """Material object from a YAML file of the refractiveindex.info database (Polyanskiy, Sci. Data 11, 94 (2024);
    check the licence of the data before using them in examples). Supported: ``tabulated nk``, ``tabulated n``,
    ``tabulated n`` + ``tabulated k`` (common nodes: union inside the overlap, linear interpolation) as a table
    under ``"n"``; ``formula 1`` (Sellmeier with C_k^2) and ``formula 2`` (Sellmeier with C_k) with the validity
    range as ``range_um``. Other formulas raise ``NotImplementedError``. Needs PyYAML
    (``pip install fouriergrating[materials]``)."""
    try:
        import yaml
    except ImportError as e:  # pragma: no cover - optional dependency
        raise ImportError("from_refractiveindex needs PyYAML: pip install fouriergrating[materials]") from e
    s = os.fspath(path_or_text)
    if "\n" not in s and os.path.isfile(s):
        with open(s, encoding="utf-8") as f:
            s = f.read()
    doc = yaml.safe_load(s)
    if not isinstance(doc, dict) or not isinstance(doc.get("DATA"), list) or not doc["DATA"]:
        raise ValueError("not a refractiveindex.info material file (DATA missing)")
    data = doc["DATA"]
    types = [str(d.get("type", "")).strip() for d in data]
    tab_n = next((d for d, t in zip(data, types) if t == "tabulated n"), None)
    tab_k = next((d for d, t in zip(data, types) if t == "tabulated k"), None)
    for d, t in zip(data, types):
        if t == "tabulated nk":
            return {"n": {"table": {"data": [[r[0], r[1], r[2]] for r in _rows(d["data"])]}}}
        if t.startswith("formula"):
            num = int(t.split()[1])
            if num in (1, 2):
                return _sellmeier(d, num)
            raise NotImplementedError(f"refractiveindex.info formula {num} (supported: 1, 2)")
    if tab_n is not None and tab_k is not None:
        n, k = _rows(tab_n["data"]), _rows(tab_k["data"])
        xn, yn = [r[0] for r in n], [r[1] for r in n]
        xk, yk = [r[0] for r in k], [r[1] for r in k]
        lo, hi = max(xn[0], xk[0]), min(xn[-1], xk[-1])
        if lo >= hi:
            raise ValueError("tabulated n and tabulated k do not overlap")
        xs = sorted({x for x in xn + xk if lo <= x <= hi})
        return {"n": {"table": {"data": [[x, _interp(x, xn, yn), _interp(x, xk, yk)] for x in xs]}}}
    if tab_n is not None:
        return {"n": {"table": {"data": [[r[0], r[1]] for r in _rows(tab_n["data"])]}}}
    raise NotImplementedError(f"refractiveindex.info data types {types} (supported: tabulated nk, tabulated n[, k], "
                              "formula 1, formula 2)")
