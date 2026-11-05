"""FourierGrating: Fourier modal method for 1D/2D periodic gratings (RK differential method and RCWA).

The configuration is the JSON format of ``docs/config_format.md`` as a ``dict`` (units µm and degrees);
CLI (``fgsolve``) and Python share one parser and one semantics. See ``docs/python.md``.

    import fouriergrating as fg
    cfg = fg.load_config("examples/v3_lamellar_dielectric.json")
    res = fg.solve(cfg)                     # first scan point
    scan = fg.solve_scan(cfg)               # all sweep values and scan points
    sim = fg.Simulation(cfg)                # layer models built once ...
    sim.solve(lambda_um=0.6328)             # ... many computations
    sim.set_parameter("/layers/0/thickness_um", 0.25)
"""

from __future__ import annotations

import copy
import json
import os
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Union

import numpy as np

from ._loader import load_core

_core = load_core()

FourierGratingError = _core.FourierGratingError
__version__: str = _core.version()

__all__ = [
    "FourierGratingError",
    "ORDER_DTYPE",
    "Result",
    "ScanResult",
    "Simulation",
    "core_variant",
    "cuda_available",
    "get_backend",
    "load_config",
    "near_field",
    "set_backend",
    "solve",
    "solve_scan",
]

#: dtype of ``Result.reflected`` / ``Result.transmitted`` (one row per propagating order)
ORDER_DTYPE = np.dtype(
    [("n", "i4"), ("m", "i4"), ("eff", "f8"), ("eff_s", "f8"), ("eff_p", "f8"), ("amp_s", "c16"), ("amp_p", "c16")]
)

Config = Union[Mapping[str, Any], str, "os.PathLike[str]"]

_backend: dict[str, Any] = {}


# ----------------------------------------------------------------------------- configuration
def load_config(path: Union[str, "os.PathLike[str]"]) -> dict:
    """Reads a JSON configuration; ``output.dir`` is resolved relative to the file."""
    path = os.fspath(path)
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    out = cfg.get("output")
    if isinstance(out, dict) and isinstance(out.get("dir"), str) and not os.path.isabs(out["dir"]):
        out["dir"] = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), out["dir"]))
    return cfg


def set_backend(preference: str = "auto", threshold: int = 1000) -> None:
    """Default for configurations without ``"backend"``: ``"auto"`` | ``"cpu"`` | ``"cuda"``."""
    if preference not in ("auto", "cpu", "cuda"):
        raise ValueError('backend must be "auto", "cpu" or "cuda"')
    _backend.update(backend=preference, backend_threshold=int(threshold))


def get_backend() -> dict:
    return {"backend": "auto", "backend_threshold": 1000, **_backend}


def cuda_available() -> bool:
    """True if the loaded core contains the CUDA backend and a GPU is usable."""
    return bool(_core.cuda_available())


def core_variant() -> str:
    """``"cpu"`` (package fouriergrating) or ``"cuda12"`` (package fouriergrating-cuda12)."""
    return str(_core.variant())


def _to_json(cfg: Config) -> str:
    if isinstance(cfg, (str, os.PathLike)):
        cfg = load_config(cfg)
    c = copy.deepcopy(dict(cfg))
    for k, v in _backend.items():
        c.setdefault(k, v)
    return json.dumps(c)


# ----------------------------------------------------------------------------- results
@dataclass
class Result:
    """Efficiencies of one incidence. ``reflected``/``transmitted``: structured arrays (ORDER_DTYPE)."""

    reflected: np.ndarray
    transmitted: np.ndarray
    absorption: float
    diagnostics: dict
    incidence: Optional[dict] = None

    def efficiency(self, side: str, n: int, m: int = 0) -> float:
        """Efficiency of order (n, m) on side ``"R"`` or ``"T"`` (0 if not propagating)."""
        a = self.reflected if side.upper() == "R" else self.transmitted
        hit = a[(a["n"] == n) & (a["m"] == m)]
        return float(hit["eff"][0]) if len(hit) else 0.0


def _orders(rows: list) -> np.ndarray:
    return np.array([tuple(r) for r in rows], dtype=ORDER_DTYPE)


def _result(d: dict) -> Result:
    return Result(_orders(d["reflected"]), _orders(d["transmitted"]), float(d["absorption"]), d["diagnostics"],
                  d.get("incidence"))


@dataclass
class ScanResult:
    """All sweep values and scan points. ``columns``: one row per order (as efficiencies.csv)."""

    points: list
    sweep_parameter: str = ""
    columns: dict = field(default_factory=dict)

    def to_pandas(self):  # pragma: no cover - optional dependency
        import pandas as pd

        return pd.DataFrame(self.columns)


def _columns(points: list) -> dict:
    keys = ["sweep_value", "lambda_um", "phi_deg", "psi_deg", "delta_deg", "side", "n", "m", "efficiency",
            "eff_s", "eff_p", "amp_s", "amp_p"]
    cols: dict = {k: [] for k in keys}
    for sv, r in points:
        inc = r.incidence
        for side, arr in (("R", r.reflected), ("T", r.transmitted)):
            for o in arr:
                for k, v in (("sweep_value", sv), ("lambda_um", inc["lambda_um"]), ("phi_deg", inc["phi_deg"]),
                             ("psi_deg", inc["psi_deg"]), ("delta_deg", inc["delta_deg"]), ("side", side),
                             ("n", o["n"]), ("m", o["m"]), ("efficiency", o["eff"]), ("eff_s", o["eff_s"]),
                             ("eff_p", o["eff_p"]), ("amp_s", o["amp_s"]), ("amp_p", o["amp_p"])):
                    cols[k].append(v)
    dtypes = {"side": "U1", "n": "i4", "m": "i4", "amp_s": "c16", "amp_p": "c16"}
    return {k: np.asarray(v, dtype=dtypes.get(k, "f8")) for k, v in cols.items()}


# ----------------------------------------------------------------------------- computations
def solve(cfg: Config) -> Result:
    """Efficiencies for the first scan point of the configuration ("sweep" is ignored)."""
    return _result(_core.solve(_to_json(cfg)))


def solve_scan(cfg: Config) -> ScanResult:
    """All sweep values and scan points, as ``fgsolve`` (scan points run in parallel, OpenMP)."""
    runs = _core.solve_scan(_to_json(cfg))
    points = []
    for run in runs:
        for p in run["points"]:
            points.append((run["sweep_value"], _result(p)))
    sweep = runs[0]["sweep_parameter"] if runs else ""
    return ScanResult(points=[r for _, r in points], sweep_parameter=sweep, columns=_columns(points))


def near_field(cfg: Config) -> dict:
    """Near-field cuts (block "nearfield") for the first scan point: ``{"horizontal", "vertical"}``,
    each a dict with the axes, Ex..Hz (complex), absE2, absH2, Sy, material and meta (as meta.json)."""
    d = _core.near_field(_to_json(cfg))
    for k in ("horizontal", "vertical"):
        d[k]["meta"] = json.loads(d[k]["meta"])
    return d


class Simulation:
    """Prepared computation: the layer models are built once and kept between ``solve`` calls;
    ``set_parameter`` rebuilds only what depends on the changed value. One object per thread."""

    def __init__(self, cfg: Config):
        self._sim = _core.Simulation(_to_json(cfg))

    def solve(self, lambda_um: Optional[float] = None, phi_deg: Optional[float] = None,
              psi_deg: Optional[float] = None, delta_deg: Optional[float] = None) -> Result:
        """Efficiencies; omitted incidence values come from the first scan point of the configuration."""
        return _result(self._sim.solve(lambda_um, phi_deg, psi_deg, delta_deg))

    def set_parameter(self, pointer: str, value: float) -> None:
        """Sets a number of the configuration by JSON pointer (as in sweeps), e.g. "/layers/0/thickness_um"."""
        self._sim.set_parameter(pointer, float(value))

    @property
    def config(self) -> dict:
        return json.loads(self._sim.config_json())

    @property
    def model_builds(self) -> int:
        """Number of layer models built so far."""
        return int(self._sim.model_builds)
