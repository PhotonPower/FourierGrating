"""Wavelength scan and depth sweep of a glass lamellar grating with the Python API (M12).

    python examples/python/scan_lamellar.py            # table
    python examples/python/scan_lamellar.py --plot     # plus matplotlib figure

The Simulation object keeps the layer models: changing the wavelength costs no rebuild, changing the
depth rebuilds only the lamellar layer.
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

    cfg = fg.load_config(HERE.parent / "v3_lamellar_dielectric.json")
    sim = fg.Simulation(cfg)
    lambdas = np.linspace(0.55, 0.75, 21)
    depths = [0.2, 0.3, 0.4]
    t0 = np.empty((len(depths), len(lambdas)))
    for i, d in enumerate(depths):
        sim.set_parameter("/layers/0/thickness_um", d)
        for k, lam in enumerate(lambdas):
            t0[i, k] = sim.solve(lambda_um=lam).efficiency("T", 0)
    print(f"model builds: {sim.model_builds} (3 layer rebuilds + initial build)")
    print("lambda_um  " + "  ".join(f"T0(d={d:.2f})" for d in depths))
    for k, lam in enumerate(lambdas):
        print(f"{lam:9.4f}  " + "  ".join(f"{t0[i, k]:11.6f}" for i in range(len(depths))))

    # the same scan via the configuration (parallel over the scan points, as fgsolve)
    scan_cfg = dict(cfg, incidence=dict(cfg["incidence"], lambda_um={"from": 0.55, "to": 0.75, "steps": 21}))
    scan = fg.solve_scan(scan_cfg)
    assert abs(scan.points[-1].efficiency("T", 0) - t0[-1, -1]) < 1e-12

    if args.plot:
        import matplotlib.pyplot as plt

        for i, d in enumerate(depths):
            plt.plot(lambdas, t0[i], label=f"d = {d} µm")
        plt.xlabel("λ (µm)")
        plt.ylabel("T₀")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
