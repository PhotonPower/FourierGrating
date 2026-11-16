"""Dispersive materials from Python (Plan section "Dispersive Materialien (M22)"): material_eps, diagnostics,
scan equivalence through the Python API and V23h (import from refractiveindex.info YAML files)."""

import copy
import json

import numpy as np
import pytest

import fouriergrating as fg

SIO2 = {"sellmeier": {"A": 1.0, "B": [0.6961663, 0.4079426, 0.8974794],
                      "C_um2": [0.00467914826, 0.0135120631, 97.9340025]}}
AG = {"drude_lorentz": {"eps_inf": 1.0, "omega_p_eV": 9.01, "f0": 0.845, "gamma0_eV": 0.048,
                        "oscillators": [{"f": 0.065, "omega_eV": 0.816, "gamma_eV": 3.886},
                                        {"f": 0.124, "omega_eV": 4.481, "gamma_eV": 0.452},
                                        {"f": 0.011, "omega_eV": 8.185, "gamma_eV": 0.065},
                                        {"f": 0.840, "omega_eV": 9.083, "gamma_eV": 0.916},
                                        {"f": 5.646, "omega_eV": 20.29, "gamma_eV": 2.419}]}}

CFG = {
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 6},
    "period_um": {"x": 0.6},
    "incidence": {"lambda_um": [0.45, 0.6, 0.75], "phi_deg": 12, "psi_deg": 0, "delta_deg": 30},
    "materials": {"air": {"n": 1.0}, "sio2": {"n": SIO2}, "ag": {"eps": AG},
                  "lnb": {"eps_principal": [SIO2, 4.98, {"cauchy": {"A": 2.1, "B_um2": 0.02}}],
                          "euler_deg": [10, 30, 0]}},
    "superstrate": "air", "substrate": "sio2",
    "layers": [{"type": "lamellar1d", "thickness_um": 0.05, "breakpoints_um": [0.0, 0.3, 0.6],
                "materials": ["ag", "air"]}],
}


def identical(a, b):
    return (np.array_equal(a.reflected, b.reflected) and np.array_equal(a.transmitted, b.transmitted)
            and a.absorption == b.absorption)


def fixed(cfg, lam):
    """The configuration at lam with the dispersive values in the layers written as numbers (from the core)."""
    c = copy.deepcopy(cfg)
    c["incidence"]["lambda_um"] = lam
    for name in ("sio2", "ag"):  # both models give eps (Sellmeier, Drude-Lorentz)
        eps = fg.material_eps(cfg, name, lam)
        c["materials"][name] = {"eps": [eps.real, eps.imag]}
    return c


def test_material_eps():
    lam = np.array([0.4, 0.5876, 1.2])
    eps = fg.material_eps(CFG, "sio2", lam)
    assert eps.shape == (3,) and eps.dtype == complex
    assert abs(np.sqrt(eps[1].real) - 1.4585) < 1e-4  # Malitson n_d
    assert isinstance(fg.material_eps(CFG, "sio2", 0.5), complex)
    ag = fg.material_eps(CFG, "ag", np.linspace(0.3, 1.5, 7))
    assert np.all(ag.imag > 0) and ag[-1].real < -50
    p = fg.material_eps(CFG, "lnb", [0.5, 0.6])
    assert p.shape == (2, 3)
    assert p[0, 1] == 4.98 and p[0, 0] == fg.material_eps(CFG, "sio2", 0.5)
    assert fg.material_eps(CFG, "air", [0.5, 0.6]).tolist() == [1.0, 1.0]  # fixed materials too


def test_scan_equals_fixed_values_and_diagnostics():
    scan = fg.solve_scan(CFG)
    sim = fg.Simulation(CFG)
    for lam, r in zip([0.45, 0.6, 0.75], scan.points):
        ref = fg.solve(fixed(CFG, lam))
        assert identical(r, ref)
        assert identical(sim.solve(lambda_um=lam), ref)
        mats = {(m["name"], m["component"]): m for m in r.diagnostics["materials"]}
        assert set(mats) == {("sio2", 0), ("ag", 0), ("lnb", 1), ("lnb", 2), ("lnb", 3)}
        assert mats[("ag", 0)]["eps"] == fg.material_eps(CFG, "ag", lam)
    assert identical(sim.solve(lambda_um=0.45), scan.points[0])  # back to the first wavelength
    assert sim.model_builds == 4  # the lamellar layer again at 0.6, 0.75 and 0.45


def test_gradient_chain_rule():
    sim = fg.Simulation({**CFG, "incidence": {**CFG["incidence"], "lambda_um": 0.6}})
    obj = {("R", 0, 0): 1.0, ("T", -1, 0): 0.5}
    params = ["/incidence/lambda_um", "/materials/ag/eps/drude_lorentz/omega_p_eV", "/materials/sio2/n/sellmeier/B/0"]
    _, g = sim.value_and_grad(obj, params)
    for k, p in enumerate(params):
        h = 1e-6
        vals = []
        for s in (h, -h):
            c = json.loads(json.dumps(sim.config))
            node = c
            keys = p.strip("/").split("/")
            for key in keys[:-1]:
                node = node[int(key)] if isinstance(node, list) else node[key]
            last = int(keys[-1]) if isinstance(node, list) else keys[-1]
            node[last] = node[last] + s
            vals.append(fg.Simulation(c).value_and_grad(obj, [])[0])
        assert g[k] == pytest.approx((vals[0] - vals[1]) / (2 * h), rel=1e-6, abs=1e-8), p


def test_errors_name_the_pointer():
    c = copy.deepcopy(CFG)
    c["materials"]["ag"] = {"n": {"table": {"data": [[0.5, 0.1, 3.0], [0.7, 0.14, 4.6]]}}}
    with pytest.raises(fg.FourierGratingError, match=r"/materials/ag/n/table: lambda 0\.45 um outside"):
        fg.solve_scan(c)


def test_load_config_resolves_table_files(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "ag.csv").write_text("lambda_um,n,k\n0.4,0.17,1.95\n0.8,0.14,5.4\n", encoding="utf-8")
    c = copy.deepcopy(CFG)
    c["incidence"]["lambda_um"] = 0.6
    c["materials"]["ag"] = {"n": {"table": {"file": "data/ag.csv"}}}
    (tmp_path / "cfg.json").write_text(json.dumps(c), encoding="utf-8")
    cfg = fg.load_config(tmp_path / "cfg.json")
    assert fg.material_eps(cfg, "ag", 0.6) == pytest.approx(complex(0.155, 3.675) ** 2, rel=1e-14)


# ----------------------------------------------------------------------------- V23h: refractiveindex.info import

NK = """REFERENCES: "test data"
DATA:
  - type: tabulated nk
    data: |
        0.40 1.47 1.95
        0.50 0.97 1.87
        0.60 0.25 2.98
"""
N_ONLY = """DATA:
  - type: tabulated n
    data: |
        0.40 1.52
        0.60 1.50
"""
N_AND_K = """DATA:
  - type: tabulated n
    data: |
        0.40 1.60
        0.50 1.55
        0.80 1.50
  - type: tabulated k
    data: |
        0.45 0.020
        0.70 0.010
        0.90 0.005
"""
F1 = """DATA:
  - type: formula 1
    wavelength_range: 0.21 6.7
    coefficients: 0 0.6961663 0.0684043 0.4079426 0.1162414 0.8974794 9.896161
"""
F2 = """DATA:
  - type: formula 2
    wavelength_range: 0.3 2.5
    coefficients: 0.5 1.0 0.01
"""


def test_import_refractiveindex(tmp_path):
    pytest.importorskip("yaml")  # optional dependency (fouriergrating[materials])
    imp = fg.materials.from_refractiveindex
    assert imp(NK) == {"n": {"table": {"data": [[0.4, 1.47, 1.95], [0.5, 0.97, 1.87], [0.6, 0.25, 2.98]]}}}
    assert imp(N_ONLY) == {"n": {"table": {"data": [[0.4, 1.52], [0.6, 1.5]]}}}
    nk = imp(N_AND_K)["n"]["table"]["data"]  # union of the nodes inside the overlap [0.45, 0.8]
    assert [r[0] for r in nk] == [0.45, 0.5, 0.7, 0.8]
    assert nk[0][1] == pytest.approx(1.575) and nk[1][2] == pytest.approx(0.018)
    f1 = imp(F1)
    assert f1 == {"n": {"sellmeier": {"A": 1.0, "B": [0.6961663, 0.4079426, 0.8974794],
                                      "C_um2": [0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2],
                                      "range_um": [0.21, 6.7]}}}
    assert imp(F2) == {"n": {"sellmeier": {"A": 1.5, "B": [1.0], "C_um2": [0.01], "range_um": [0.3, 2.5]}}}
    # evaluated by the core like the handwritten model; from a file path too
    (tmp_path / "sio2.yml").write_text(F1, encoding="utf-8")
    c = copy.deepcopy(CFG)
    c["materials"]["imp"] = imp(tmp_path / "sio2.yml")
    c["materials"]["hand"] = {"n": {"sellmeier": {"A": 1.0, "B": [0.6961663, 0.4079426, 0.8974794],
                                                  "C_um2": [0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2],
                                                  "range_um": [0.21, 6.7]}}}
    lam = np.linspace(0.25, 6.5, 9)
    assert np.array_equal(fg.material_eps(c, "imp", lam), fg.material_eps(c, "hand", lam))
    assert abs(np.sqrt(fg.material_eps(c, "imp", 0.5876).real) - 1.4585) < 1e-4
    with pytest.raises(fg.FourierGratingError, match="range_um"):
        fg.material_eps(c, "imp", 0.2)
    with pytest.raises(NotImplementedError, match="formula 3"):
        imp("DATA:\n  - type: formula 3\n    coefficients: 1 2 3\n")
