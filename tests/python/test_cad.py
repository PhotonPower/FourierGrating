"""fg.cad (M14): check, preview, tool description, cad layers from Python."""

import copy

import pytest

import fouriergrating as fg
from conftest import ROOT

pytestmark = pytest.mark.skipif(not fg.cad.available(), reason="core built without FG_ENABLE_NURBFORGE")

CFG = {
    "dimension": "2D", "solver": "rcwa", "truncation": {"Nx": 2, "Nz": 2},
    "period_um": {"x": 0.5, "z": 0.5},
    "incidence": {"lambda_um": 0.633, "phi_deg": 10, "psi_deg": 0, "delta_deg": 0},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "pillar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "layers": [{"type": "cad", "file": str(ROOT / "tests" / "data" / "cad" / "cylinder.step"),
                "cell_origin_um": [0.0, -0.5], "bodies": [{"name": "Pillar", "material": "pillar"}],
                "normal_field": {"grid": 64}}],
}


def test_check_and_solve():
    r = fg.cad.check(None, CFG)
    assert r["ok"] and r["errors"] == []
    assert [b["kind"] for b in r["bands"]] == ["extrusion"]
    res = fg.solve(CFG)
    assert abs(res.diagnostics["energy_balance"] - 1.0) < 1e-10
    # the STEP text instead of the file gives the same check
    text = (ROOT / "tests" / "data" / "cad" / "cylinder.step").read_text()
    assert fg.cad.check(text, CFG)["bands"] == r["bands"]


def test_findings_and_tools():
    bad = copy.deepcopy(CFG)
    bad["layers"][0]["bodies"] = [{"name": "Other", "material": "pillar"}]
    r = fg.cad.check(None, bad)
    assert not r["ok"]
    assert {e["code"] for e in r["errors"]} == {"BODY_UNASSIGNED"}
    with pytest.raises(fg.FourierGratingError):
        fg.solve(bad)
    svg = fg.cad.preview(None, CFG, 0.1)
    assert svg.startswith("<svg") and "pillar" in svg
    spec = fg.cad.llm_tool_spec()
    assert spec["name"] == "fouriergrating_check_cad" and set(spec["parameters"]["required"]) == {"step", "config"}


def test_load_config_resolves_step_file():
    cfg = fg.load_config(ROOT / "examples" / "cad_pillar_array.json")
    assert cfg["layers"][0]["file"].endswith("pillar_array.step")
    assert fg.cad.check(None, cfg)["ok"]
