"""CAD slicer (M14): STEP bodies -> layers. Geometry check for an LLM tool loop, section preview.

The STEP text normally comes from Nurbforge (``nf.to_step(shells, unit="um", names=[...])``); the
configuration is the usual FourierGrating dict with a layer ``{"type": "cad", ...}`` (docs/cad.md).

    import fouriergrating as fg
    report = fg.cad.check(step_text, cfg)       # {"ok", "errors", "warnings", "bands", ...}
    svg = fg.cad.preview(step_text, cfg, y_um=0.1)
    tool = fg.cad.llm_tool_spec()               # JSON schema for function calling
"""

from __future__ import annotations

import json
import os
from typing import Any, Mapping, Optional, Union

from . import _core, _to_json

Step = Union[str, "os.PathLike[str]", None]

__all__ = ["available", "check", "llm_tool_spec", "preview"]


def available() -> bool:
    """True if the core was built with the CAD slicer (FG_ENABLE_NURBFORGE)."""
    return bool(_core.cad_available())


def _step_text(step: Step) -> str:
    if step is None:
        return ""
    if isinstance(step, os.PathLike) or (isinstance(step, str) and not step.lstrip().startswith("ISO-10303-21")):
        with open(os.fspath(step), encoding="utf-8", errors="replace") as f:
            return f.read()
    return str(step)


def _base_dir(cfg: Any) -> str:
    return os.path.dirname(os.path.abspath(os.fspath(cfg))) if isinstance(cfg, (str, os.PathLike)) else ""


def check(step: Step, cfg: Union[Mapping[str, Any], str, "os.PathLike[str]"]) -> dict:
    """Geometry check and band decomposition: dict with ok, errors, warnings (code, message, where_um)
    and bands. `step`: STEP text, a STEP file or None (the file of the cad layer)."""
    return json.loads(_core.cad_report(_to_json(cfg), _base_dir(cfg), _step_text(step)))


def preview(step: Step, cfg: Union[Mapping[str, Any], str, "os.PathLike[str]"], y_um: float) -> str:
    """SVG of the section at height y_um (FG coordinates, one cell, materials coloured)."""
    return str(_core.cad_preview(_to_json(cfg), _base_dir(cfg), _step_text(step), float(y_um)))


def llm_tool_spec() -> dict:
    """Tool description (name, description, JSON-schema parameters) of `check` for an LLM tool loop,
    e.g. Nurbforge's Ollama agent (fouriergrating_check_cad(step, config))."""
    return json.loads(_core.cad_llm_tool_spec())


def _tool_call(arguments: Optional[Mapping[str, Any]]) -> dict:
    """Runs the tool for a function call with the arguments of llm_tool_spec()."""
    a = dict(arguments or {})
    return check(a.get("step"), a.get("config", {}))
