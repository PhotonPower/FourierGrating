"""Spectrum of a silver grating with a dispersive (Drude-Lorentz) material model (M22).

    python examples/python/dispersive_scan.py            # table
    python examples/python/dispersive_scan.py --plot     # plus matplotlib figure

The material values are functions of the wavelength and are evaluated by the core at every scan point
(Plan section "Dispersive Materialien"). A Simulation rebuilds only the layer that contains silver when the
wavelength changes; fg.material_eps evaluates a model for plots.
"""

import argparse
import pathlib

import numpy as np

import fouriergrating as fg

HERE = pathlib.Path(__file__).resolve().parent


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    cfg = fg.load_config(HERE.parent / "dispersive_metal_grating.json")
    scan = fg.solve_scan(cfg)  # 0.4 .. 1.0 um, parallel over the scan points
    lambdas = np.array([p.incidence["lambda_um"] for p in scan.points])
    r0 = np.array([p.efficiency("R", 0) for p in scan.points])
    absorbed = np.array([p.absorption for p in scan.points])
    eps = fg.material_eps(cfg, "ag", lambdas)
    print("lambda_um        R0   absorbed   eps(Ag)")
    for lam, r, a, e in zip(lambdas, r0, absorbed, eps):
        print(f"{lam:9.3f}  {r:8.5f}  {a:9.5f}   {e.real:9.3f}{e.imag:+8.3f}i")

    # the same with a Simulation: one rebuild per wavelength change (only the lamellar layer)
    sim = fg.Simulation(cfg)
    for lam in lambdas[:5]:
        res = sim.solve(lambda_um=float(lam))
        assert abs(res.efficiency("R", 0) - r0[list(lambdas).index(lam)]) < 1e-12
    print(f"model builds: {sim.model_builds} (initial build + 4 wavelength changes)")
    # n and eps at the scan point are in the diagnostics
    print("diagnostics:", scan.points[0].diagnostics["materials"])

    if args.plot:
        import matplotlib.pyplot as plt

        fig, (a1, a2) = plt.subplots(2, 1, sharex=True)
        a1.plot(lambdas, r0, label="R₀")
        a1.plot(lambdas, absorbed, label="absorbed")
        a1.legend()
        a2.plot(lambdas, eps.real, label="Re ε (Ag)")
        a2.plot(lambdas, eps.imag, label="Im ε (Ag)")
        a2.set_xlabel("λ (µm)")
        a2.legend()
        plt.show()


if __name__ == "__main__":
    main()
