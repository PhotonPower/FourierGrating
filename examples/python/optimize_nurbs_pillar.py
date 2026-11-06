"""NURBS pillar optimization with PyTorch (Plan section "M19", example 2; demonstration).

A dielectric pillar (n = 2) on glass whose cross-section is a closed quadratic NURBS curve (contour2d,
initially an exact circle). The inner control points are optimized with Adam to maximize the transmitted
(-1, 0) order; the gradient comes from the discrete adjoint through fouriergrating.torch. The objective
should fall monotonically (no reference value). Requires torch.

    python examples/python/optimize_nurbs_pillar.py
"""

import math
import sys

import fouriergrating as fg

try:
    import torch
except ImportError:  # pragma: no cover
    print("this example needs PyTorch (pip install torch)")
    sys.exit(0)

from fouriergrating.torch import efficiencies

R, CX, CZ = 0.15, 0.25, 0.25
s2 = math.sqrt(0.5)
# rational quadratic circle: 9 points (first = last), weights 1, 1/sqrt 2, ...
POINTS = [[CX + R, CZ], [CX + R, CZ + R], [CX, CZ + R], [CX - R, CZ + R], [CX - R, CZ], [CX - R, CZ - R],
          [CX, CZ - R], [CX + R, CZ - R], [CX + R, CZ]]
CFG = {
    "dimension": "2D", "solver": "rcwa", "truncation": {"Nx": 3, "Nz": 3},
    "period_um": {"x": 0.5, "z": 0.5},
    "incidence": {"lambda_um": 0.6328, "phi_deg": 10, "psi_deg": 0, "delta_deg": 0},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "pillar": {"n": 2.0}},
    "superstrate": "air", "substrate": "glass",
    "layers": [{"type": "contour2d", "thickness_um": 0.3, "background": "air", "factorization": "laurent",
                "regions": [{"material": "pillar", "nurbs": {
                    "degree": 2, "knots": [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4],
                    "points_um": POINTS, "weights": [1, s2, 1, s2, 1, s2, 1, s2, 1]}}]}],
}
# inner control points 1..7 (x and z); point 0 = point 8 stays (closed curve)
PARAMS = [f"/layers/0/regions/0/nurbs/points_um/{k}/{c}" for k in range(1, 8) for c in (0, 1)]


def main() -> int:
    sim = fg.Simulation(CFG)
    x = torch.tensor([POINTS[k][c] for k in range(1, 8) for c in (0, 1)], dtype=torch.float64, requires_grad=True)
    opt = torch.optim.Adam([x], lr=2e-3)
    history = []
    for it in range(25):
        opt.zero_grad()
        e = efficiencies(sim, PARAMS, x, [("T", -1, 0)])
        loss = -e[0]
        loss.backward()
        opt.step()
        history.append(loss.item())
        print(f"iteration {it:2d}: T(-1, 0) = {-history[-1]:.6f}")
    monotone = all(b <= a + 1e-9 for a, b in zip(history, history[1:]))
    print("objective falls monotonically:", "yes" if monotone else "no")
    return 0 if monotone else 1


if __name__ == "__main__":
    sys.exit(main())
