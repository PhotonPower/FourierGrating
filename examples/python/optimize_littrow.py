"""Littrow grating optimization (Plan section "M19", example 1).

Binary aluminium grating (period 0.5 um, lambda 0.6328 um, TM) in the -1st-order Littrow mount
(sin phi = lambda / (2 d)): the depth and the ridge width are optimized for the -1st reflected order with
L-BFGS-B from a fixed start, using the gradient of the discrete adjoint (Simulation.value_and_grad). The
result is checked against a raster scan: the optimum must reach the best raster value within 1e-3.

    python examples/python/optimize_littrow.py
"""

import math
import sys

import numpy as np
from scipy.optimize import minimize

import fouriergrating as fg

LAMBDA, D = 0.6328, 0.5
CFG = {
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 12},
    "period_um": {"x": D},
    "incidence": {"lambda_um": LAMBDA, "phi_deg": math.degrees(math.asin(LAMBDA / (2 * D))), "delta_deg": 90},
    "materials": {"air": {"n": 1.0}, "al": {"n": [1.38, 7.62]}},
    "superstrate": "air", "substrate": "al",
    "layers": [{"type": "lamellar1d", "thickness_um": 0.15, "breakpoints_um": [0.0, 0.25, D],
                "materials": ["al", "air"]}],
}
PARAMS = ["/layers/0/thickness_um", "/layers/0/breakpoints_um/1"]
BOUNDS = [(0.03, 0.40), (0.05, 0.45)]
ORDER = {("R", -1, 0): 1.0}


def efficiency(sim, x):
    for p, v in zip(PARAMS, x):
        sim.set_parameter(p, float(v))
    return sim.solve().efficiency("R", -1)


def main() -> int:
    sim = fg.Simulation(CFG)
    # raster scan
    hs = np.linspace(*BOUNDS[0], 36)
    ws = np.linspace(*BOUNDS[1], 41)
    grid = np.array([[efficiency(sim, (h, w)) for w in ws] for h in hs])
    i, j = np.unravel_index(np.argmax(grid), grid.shape)
    print(f"raster scan ({grid.size} points): e_-1 = {grid[i, j]:.6f} at h = {hs[i]:.4f}, w = {ws[j]:.4f}")

    # gradient-based optimization from a fixed start
    calls = 0

    def fun(x):
        nonlocal calls
        calls += 1
        for p, v in zip(PARAMS, x):
            sim.set_parameter(p, float(v))
        value, grad = sim.value_and_grad(ORDER, PARAMS)
        return -value, -grad

    res = minimize(fun, x0=[0.2, 0.3], jac=True, method="L-BFGS-B", bounds=BOUNDS, options={"ftol": 1e-12})
    e_opt = -res.fun
    print(f"L-BFGS-B ({calls} value+gradient evaluations): e_-1 = {e_opt:.6f} at h = {res.x[0]:.5f}, "
          f"w = {res.x[1]:.5f}")
    ok = e_opt >= grid[i, j] - 1e-3
    print("optimum >= raster optimum - 1e-3:", "yes" if ok else "NO")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
