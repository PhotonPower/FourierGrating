"""Near fields from Python = fgsolve output files (small examples only, Plan section "M12: Tests")."""

import json
import subprocess

import numpy as np
import pytest

import fouriergrating as fg
from conftest import ROOT

SMALL = ["nearfield_lamellar_1d", "nearfield_pillars_2d"]


@pytest.mark.parametrize("name", SMALL)
def test_structure(name):
    nf = fg.near_field(fg.load_config(ROOT / "examples" / f"{name}.json"))
    for cut in ("horizontal", "vertical"):
        d = nf[cut]
        shape = d["Ex"].shape
        assert d["absE2"].shape == shape and d["Sy"].shape == shape and d["Hz"].dtype == np.complex128
        assert np.allclose(d["absE2"], abs(d["Ex"]) ** 2 + abs(d["Ey"]) ** 2 + abs(d["Ez"]) ** 2, rtol=1e-14)
        assert d["meta"]["cut"] == cut and "absE2" in d["meta"]["quantities"]
        assert shape[-1] == len(d["axis1_um"])


@pytest.mark.parametrize("name", SMALL)
def test_equals_fgsolve(name, fgsolve, tmp_path):
    cfg = fg.load_config(ROOT / "examples" / f"{name}.json")
    cfg["output"] = {"dir": str(tmp_path)}
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(cfg))
    subprocess.run([fgsolve, str(path)], check=True, capture_output=True)
    nf = fg.near_field(cfg)
    for cut in ("horizontal", "vertical"):
        ref_dir = tmp_path / f"nearfield_{cut}"
        for q in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "absE2", "Sy"):
            ref = np.load(ref_dir / f"{q}.npy")
            scale = np.max(np.abs(ref)) or 1.0
            assert np.max(np.abs(nf[cut][q] - ref)) <= 1e-12 * scale, (cut, q)
        meta = json.loads((ref_dir / "meta.json").read_text())
        assert meta == nf[cut]["meta"]
