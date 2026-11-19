"""profile1d "nurbs" through Python (Plan section "M24: Tests und Abnahme", V25h): the JSON goes through unchanged,
a lossless overhang conserves energy, Simulation.set_parameter moves a control point, errors carry the pointer and
the code."""

import copy

import numpy as np
import pytest

import fouriergrating as fg

# mushroom: stem and an overhanging cap of quadratic pieces (interior knots of multiplicity 2)
PIECES = [[[0.0, 0.02], [0.065, 0.02], [0.13, 0.02]], [[0.13, 0.02], [0.14, 0.06], [0.12, 0.09]],
          [[0.12, 0.09], [0.05, 0.10], [0.08, 0.13]], [[0.08, 0.13], [0.15, 0.19], [0.22, 0.13]],
          [[0.22, 0.13], [0.25, 0.10], [0.18, 0.09]], [[0.18, 0.09], [0.16, 0.06], [0.17, 0.02]],
          [[0.17, 0.02], [0.235, 0.02], [0.3, 0.02]]]
POINTS = [PIECES[0][0]] + [q for p in PIECES for q in p[1:]]
KNOTS = [0, 0, 0] + [i for i in range(1, 7) for _ in range(2)] + [7, 7, 7]

CFG = {
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 8},
    "period_um": {"x": 0.3},
    "incidence": {"lambda_um": 0.633, "phi_deg": 20, "psi_deg": 30, "delta_deg": 35},
    "rcwa": {"slices": 40},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "pillar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "layers": [
        {"type": "homogeneous", "thickness_um": 0.05, "material": "glass"},
        {"type": "profile1d", "shape": "nurbs", "below": "pillar", "above": "air",
         "nurbs": {"degree": 2, "knots": KNOTS, "points_um": POINTS}},
    ],
}


def total(r):
    return r.reflected["eff"].sum() + r.transmitted["eff"].sum()


def test_overhang_is_lossless():
    r = fg.solve(CFG)
    assert abs(total(r) - 1) < 1e-10
    # the depth defaults to the top of the curve (0.16 um); an explicit larger depth adds air above it
    cfg = copy.deepcopy(CFG)
    cfg["layers"][1]["depth_um"] = 0.2
    assert abs(total(fg.solve(cfg)) - 1) < 1e-10
    assert not np.allclose(fg.solve(cfg).reflected["eff"], r.reflected["eff"])


def test_set_parameter_moves_a_control_point():
    sim = fg.Simulation(CFG)
    sim.set_parameter("/layers/1/nurbs/points_um/7/1", 0.2)  # the top of the cap
    moved = copy.deepcopy(CFG)
    moved["layers"][1]["nurbs"]["points_um"][7][1] = 0.2
    a, b = sim.solve(), fg.solve(moved)
    assert np.array_equal(a.reflected["eff"], b.reflected["eff"])


@pytest.mark.parametrize("change, code", [
    (lambda c: c["nurbs"]["points_um"][-1].__setitem__(1, 0.03), "NURBS_NOT_PERIODIC"),
    (lambda c: c["nurbs"]["points_um"][-1].__setitem__(0, 0.31), "NURBS_PERIOD"),
    (lambda c: c.__setitem__("depth_um", 0.15), "NURBS_RANGE"),
    (lambda c: c["nurbs"]["knots"].pop(), "NURBS_INPUT"),
])
def test_errors_carry_pointer_and_code(change, code):
    cfg = copy.deepcopy(CFG)
    change(cfg["layers"][1])
    with pytest.raises(fg.FourierGratingError, match=rf"/layers/1: \[{code}\]"):
        fg.solve(cfg)
