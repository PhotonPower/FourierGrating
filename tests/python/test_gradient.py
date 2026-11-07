"""Gradients (Plan section "M19: Schnittstellen"): Simulation.value_and_grad / vjp against difference
quotients, and the optional PyTorch / JAX adapters."""

import copy

import numpy as np
import pytest

import fouriergrating as fg

CFG = {
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 8},
    "period_um": {"x": 0.5},
    "incidence": {"lambda_um": 0.6328, "phi_deg": 39.26, "delta_deg": 0},
    "materials": {"air": {"n": 1.0}, "al": {"n": [1.38, 7.62]}},
    "superstrate": "air", "substrate": "al",
    "layers": [{"type": "lamellar1d", "thickness_um": 0.12, "breakpoints_um": [0.0, 0.22, 0.5],
                "materials": ["al", "air"]}],
}
PARAMS = ["/layers/0/thickness_um", "/layers/0/breakpoints_um/1"]
OBJ = {("R", -1, 0): 1.0, ("R", 0, 0): -0.25}


def objective_value(cfg, values):
    sim = fg.Simulation(cfg)
    for p, v in zip(PARAMS, values):
        sim.set_parameter(p, v)
    r = sim.solve()
    return sum(w * r.efficiency(s, n, m) for (s, n, m), w in OBJ.items())


def test_value_and_grad_matches_difference_quotient():
    sim = fg.Simulation(CFG)
    value, grad = sim.value_and_grad(OBJ, PARAMS)
    x0 = [0.12, 0.22]
    assert value == pytest.approx(objective_value(CFG, x0), abs=1e-13)
    for k in range(len(PARAMS)):
        h = 1e-6
        xp, xm = list(x0), list(x0)
        xp[k] += h
        xm[k] -= h
        num = (objective_value(CFG, xp) - objective_value(CFG, xm)) / (2 * h)
        assert grad[k] == pytest.approx(num, rel=1e-6, abs=1e-9)


def test_vjp_equals_gradient_and_list_form():
    sim = fg.Simulation(CFG)
    _, grad = sim.value_and_grad(OBJ, PARAMS)
    v = sim.vjp([("R", -1, 0, 1.0), ("r", 0, 0, -0.25)], PARAMS)
    np.testing.assert_allclose(v, grad, rtol=1e-12, atol=1e-15)
    assert sim.last_result.efficiency("R", -1) > 0


PML_CFG = {
    "dimension": "1D", "truncation": {"Nx": 16}, "period_um": {"x": 3.0},
    "incidence": {"lambda_um": 0.6328, "phi_deg": 0, "delta_deg": 40},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "bar": {"n": [1.8, 0.05]}},
    "superstrate": "air", "substrate": "glass",
    "layers": [{"type": "lamellar1d", "thickness_um": 0.3, "breakpoints_um": [0.0, 1.2, 1.75, 3.0],
                "materials": ["air", "bar", "air"]}],
    "aperiodic": {"pml": {"x": {"width_um": 0.5, "beta": 1.5}}},
}


def test_aperiodic_scattering_objective():
    """Plan v2.13, Stufe A: scattered powers of an aperiodic ridge (PML, normal incidence)."""
    obj = {("scattering", "up"): 1.0, ("scattering", "down"): 0.5}
    params = ["/layers/0/thickness_um", "/materials/bar/n/0"]
    sim = fg.Simulation(PML_CFG)
    value, grad = sim.value_and_grad(obj, params)
    sc = sim.solve().scattering
    assert value == pytest.approx(sc["up"] + 0.5 * sc["down"], rel=1e-12)
    for k, p in enumerate(params):
        h = 1e-6
        vals = []
        for sgn in (1, -1):
            s = fg.Simulation(PML_CFG)
            s.set_parameter(p, (0.3 if k == 0 else 1.8) + sgn * h)
            r = s.solve().scattering
            vals.append(r["up"] + 0.5 * r["down"])
        assert grad[k] == pytest.approx((vals[0] - vals[1]) / (2 * h), rel=1e-6)
    v = sim.vjp([("scattering", "up", 1.0), ("scattering", "down", 0.5)], params)
    np.testing.assert_allclose(v, grad, rtol=1e-12)


def test_unsupported_parameter_raises():
    sim = fg.Simulation(CFG)
    with pytest.raises(Exception):
        sim.value_and_grad(OBJ, ["/period_um/x"])  # the period is no parameter


def test_incidence_and_substrate_parameters():
    """M19 Phase 2: angle, wavelength and the (lossy) substrate index."""
    params = ["/incidence/phi_deg", "/incidence/lambda_um", "/materials/al/n/0", "/materials/al/n/1"]
    sim = fg.Simulation(CFG)
    _, grad = sim.value_and_grad(OBJ, params)
    x0 = [39.26, 0.6328, 1.38, 7.62]

    def value(x):
        s = fg.Simulation(CFG)
        for p, v in zip(params, x):
            s.set_parameter(p, v)
        r = s.solve()
        return sum(w * r.efficiency(side, n, m) for (side, n, m), w in OBJ.items())

    for k in range(len(params)):
        # small gradients (~1e-3): h = 1e-5 keeps the rounding noise of the quotient below 1e-9
        h = 1e-5 * max(1.0, abs(x0[k]))
        xp, xm = list(x0), list(x0)
        xp[k] += h
        xm[k] -= h
        assert grad[k] == pytest.approx((value(xp) - value(xm)) / (2 * h), rel=1e-6, abs=1e-8)


def test_torch_adapter():
    torch = pytest.importorskip("torch")
    from fouriergrating.torch import efficiencies

    sim = fg.Simulation(copy.deepcopy(CFG))
    x = torch.tensor([0.12, 0.22], dtype=torch.float64, requires_grad=True)
    e = efficiencies(sim, PARAMS, x, [("R", -1, 0), ("R", 0, 0)])
    (e[0] - 0.25 * e[1]).backward()
    _, grad = fg.Simulation(CFG).value_and_grad(OBJ, PARAMS)
    np.testing.assert_allclose(x.grad.numpy(), grad, rtol=1e-10)


def test_jax_adapter():
    jax = pytest.importorskip("jax")
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from fouriergrating.jax import efficiency_function

    f = efficiency_function(fg.Simulation(copy.deepcopy(CFG)), PARAMS, [("R", -1, 0), ("R", 0, 0)])
    g = jax.grad(lambda x: f(x)[0] - 0.25 * f(x)[1])(jnp.array([0.12, 0.22]))
    _, grad = fg.Simulation(CFG).value_and_grad(OBJ, PARAMS)
    np.testing.assert_allclose(np.asarray(g), grad, rtol=1e-10)


def test_littrow_example_reaches_raster_optimum(root):
    pytest.importorskip("scipy")
    import runpy

    mod = runpy.run_path(str(root / "examples" / "python" / "optimize_littrow.py"))
    assert mod["main"]() == 0


def test_nurbs_pillar_example_falls_monotonically(root):
    pytest.importorskip("torch")
    import runpy

    mod = runpy.run_path(str(root / "examples" / "python" / "optimize_nurbs_pillar.py"))
    assert mod["main"]() == 0
