"""Simulation, set_parameter, GIL release, errors (Plan section "M12: Tests und Abnahme")."""

import copy
import threading
import time

import numpy as np
import pytest

import fouriergrating as fg

CFG = {
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 10},
    "period_um": {"x": 0.55},
    "incidence": {"lambda_um": 0.63, "phi_deg": 17, "psi_deg": 10, "delta_deg": 35},
    "rcwa": {"slices": 12},
    "materials": {"air": {"n": 1.0}, "glass": {"n": [1.5, 0.01]}, "al": {"n": [1.2, 7.0]}},
    "superstrate": "air", "substrate": "glass",
    "layers": [
        {"type": "homogeneous", "thickness_um": 0.07, "material": "glass"},
        {"type": "profile1d", "shape": "sinus", "depth_um": 0.12, "below": "al", "above": "air"},
        {"type": "lamellar1d", "thickness_um": 0.15, "breakpoints_um": [0.0, 0.2, 0.55],
         "materials": ["glass", "air"]},
    ],
}


def identical(a, b):
    return (np.array_equal(a.reflected, b.reflected) and np.array_equal(a.transmitted, b.transmitted)
            and a.absorption == b.absorption)


def test_simulation_equals_solve():
    sim = fg.Simulation(CFG)
    assert identical(sim.solve(), fg.solve(CFG))
    c = copy.deepcopy(CFG)
    c["incidence"].update(lambda_um=0.7, phi_deg=25.0)
    assert identical(sim.solve(lambda_um=0.7, phi_deg=25.0), fg.solve(c))
    assert sim.model_builds == 3


def test_set_parameter_equals_rebuild():
    sim = fg.Simulation(CFG)
    sim.set_parameter("/layers/1/depth_um", 0.09)
    assert sim.model_builds == 4  # only the profile layer
    c = copy.deepcopy(CFG)
    c["layers"][1]["depth_um"] = 0.09
    assert identical(sim.solve(), fg.Simulation(c).solve())  # bit-identical
    assert sim.config["layers"][1]["depth_um"] == 0.09
    sim.set_parameter("/materials/al/n/1", 6.5)  # outside the layers: everything rebuilt
    assert sim.model_builds == 7
    c["materials"]["al"]["n"][1] = 6.5
    assert identical(sim.solve(), fg.Simulation(c).solve())


def test_errors():
    with pytest.raises(fg.FourierGratingError):
        fg.Simulation(CFG).set_parameter("/layers/9/thickness_um", 1.0)
    with pytest.raises(fg.FourierGratingError):
        fg.solve({**CFG, "solver": "fdtd"})
    with pytest.raises(fg.FourierGratingError):
        fg.solve({k: v for k, v in CFG.items() if k != "period_um"})  # missing key (JSON error)
    assert issubclass(fg.FourierGratingError, RuntimeError)


def test_two_threads_equal_sequential():
    params = [0.09, 0.14]
    seq = []
    for d in params:
        s = fg.Simulation(CFG)
        s.set_parameter("/layers/1/depth_um", d)
        seq.append(s.solve())
    out = [None, None]

    def work(i):
        s = fg.Simulation(CFG)
        s.set_parameter("/layers/1/depth_um", params[i])
        out[i] = s.solve()

    th = [threading.Thread(target=work, args=(i,)) for i in range(2)]
    for t in th:
        t.start()
    for t in th:
        t.join()
    assert all(identical(a, b) for a, b in zip(out, seq))


def test_gil_released_during_solve():
    c = copy.deepcopy(CFG)
    c["truncation"]["Nx"] = 40
    c["rcwa"]["slices"] = 40
    sim = fg.Simulation(c)
    t0 = time.perf_counter()
    sim.solve()
    single = time.perf_counter() - t0
    ticks = 0
    th = threading.Thread(target=sim.solve)
    th.start()
    while th.is_alive():  # runs only if the worker does not hold the GIL
        ticks += 1
        time.sleep(0.001)
    th.join()
    assert single > 0.05
    assert ticks >= 10


def test_backend_default():
    before = fg.get_backend()
    fg.set_backend("cpu", threshold=500)
    try:
        assert fg.get_backend() == {"backend": "cpu", "backend_threshold": 500}
        assert identical(fg.solve(CFG), fg.solve({**CFG, "backend": "cpu"}))
        with pytest.raises(ValueError):
            fg.set_backend("gpu")
    finally:
        fg.set_backend(before["backend"], before["backend_threshold"])
    assert isinstance(fg.cuda_available(), bool)
    assert fg.core_variant() in ("cpu", "cuda12")
    assert fg.__version__.count(".") == 2


def test_profile_and_threads():
    """M15: scope times from fg::prof; thread settings and the tree reduction do not change the result
    beyond 1e-12 (bit-identical for the sequential chain)."""
    r = fg.solve(CFG, profile=True)
    scopes = r.diagnostics["timings"]["scopes"]
    assert {"solve", "solve.assemble", "linalg.eig"} <= set(scopes)
    assert scopes["solve"]["count"] == 1 and scopes["solve"]["total_s"] > 0
    assert "scopes" not in fg.solve(CFG).diagnostics["timings"]
    ref = fg.solve(CFG)
    for threads in ({"scan": 1, "slices": 1}, {"slices": 3}):
        c = copy.deepcopy(CFG)
        c["threads"] = threads
        assert identical(fg.solve(c), ref)
        c["reduction"] = "tree"
        t = fg.solve(c)
        assert np.allclose(t.reflected["eff"], ref.reflected["eff"], rtol=0, atol=1e-12)
    scan = fg.solve_scan(CFG, profile=True)
    assert "scan" in scan.points[0].diagnostics["timings"]["scopes"]
