"""Python results = fgsolve for all examples (Plan section "M12: Tests und Abnahme").

Reference: the golden files (tests/regression/golden, computed exactly as fgsolve with 17 digits).
Tolerance 1e-11 of the incident power (absolute on the efficiencies, as Regression.*): the core may be
built with another compiler / C library than the golden files (Plan v2.4). Same build: identical.
"""

import csv
import math

import numpy as np
import pytest

import fouriergrating as fg
from conftest import EXAMPLES, ROOT


def golden(name):
    with open(ROOT / "tests" / "regression" / "golden" / f"{name}.csv", newline="") as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("name", EXAMPLES)
def test_scan_equals_golden(name):
    scan = fg.solve_scan(fg.load_config(ROOT / "examples" / f"{name}.json"))
    rows = [r for r in golden(name) if r["side"] != "A"]
    absorption = [float(r["efficiency"]) for r in golden(name) if r["side"] == "A"]
    c = scan.columns
    assert len(c["efficiency"]) == len(rows)
    assert len(scan.points) == len(absorption)
    for k, r in enumerate(rows):
        assert c["side"][k] == r["side"] and c["n"][k] == int(r["n"]) and c["m"][k] == int(r["m"])
        assert math.isclose(c["lambda_um"][k], float(r["lambda_um"]), rel_tol=1e-15)
        assert math.isclose(c["sweep_value"][k], float(r["sweep_value"]), rel_tol=1e-15, abs_tol=0)
        for col in ("efficiency", "eff_s", "eff_p"):
            assert abs(c[col][k] - float(r[col])) <= 1e-11, (name, k, col)
    for p, a in zip(scan.points, absorption):
        assert abs(p.absorption - a) <= 1e-11


@pytest.mark.parametrize("name", ["v3_lamellar_dielectric", "v7_cobalt_anisotropic", "sweep_lamellar_depth"])
def test_solve_is_first_scan_point(name):
    cfg = fg.load_config(ROOT / "examples" / f"{name}.json")
    first = fg.solve_scan(cfg).points[0]
    if "sweep" in cfg:  # solve() ignores the sweep: the configured value
        cfg = dict(cfg)
        cfg.pop("sweep")
        first = fg.solve_scan(cfg).points[0]
    r = fg.solve(cfg)
    assert np.array_equal(r.reflected, first.reflected)
    assert np.array_equal(r.transmitted, first.transmitted)
    assert r.absorption == first.absorption
    assert r.diagnostics["timings"]["solve_s"] >= 0


def test_to_pandas():
    pd = pytest.importorskip("pandas")
    df = fg.solve_scan(fg.load_config(ROOT / "examples" / "v3_lamellar_rk_scan.json")).to_pandas()
    assert isinstance(df, pd.DataFrame)
    assert {"lambda_um", "side", "n", "efficiency", "amp_p"} <= set(df.columns)


def test_example_script():
    import subprocess
    import sys

    r = subprocess.run([sys.executable, str(ROOT / "examples" / "python" / "scan_lamellar.py")],
                       capture_output=True, text=True, check=True)
    assert "model builds: 4" in r.stdout
