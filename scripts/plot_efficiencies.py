#!/usr/bin/env python3
"""Plot diffraction efficiencies from fgsolve's efficiencies.csv (scans and sweeps).

Usage:
  plot_efficiencies.py out/v3/efficiencies.csv                       # x axis chosen automatically
  plot_efficiencies.py out/sweep/efficiencies.csv --x sweep_value --orders R:0 R:-1 T:0 --delta 90
  plot_efficiencies.py file.csv --x lambda_um --out spectrum.png

The x axis is the first of sweep_value, lambda_um, phi_deg, psi_deg, delta_deg that varies (or --x).
Each curve is one order (side:n or side:n:m) at fixed values of the remaining scan columns
(filter them with --lambda/--phi/--psi/--delta).
"""
import argparse

import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--x", help="x column (sweep_value, lambda_um, phi_deg, psi_deg, delta_deg)")
    ap.add_argument("--orders", nargs="*", help="orders as side:n or side:n:m, e.g. R:0 T:-1 (default: all)")
    ap.add_argument("--quantity", default="efficiency", help="efficiency | eff_s | eff_p")
    for c in ("lambda", "phi", "psi", "delta"):
        ap.add_argument(f"--{c}", type=float, help=f"keep rows with this {c} value")
    ap.add_argument("--out", help="write the figure instead of showing it")
    a = ap.parse_args()

    df = pd.read_csv(a.csv)
    for c, col in (("lambda", "lambda_um"), ("phi", "phi_deg"), ("psi", "psi_deg"), ("delta", "delta_deg")):
        v = getattr(a, c)
        if v is not None:
            df = df[(df[col] - v).abs() < 1e-9]
    cands = [c for c in ("sweep_value", "lambda_um", "phi_deg", "psi_deg", "delta_deg") if c in df.columns]
    x = a.x or next((c for c in cands if df[c].nunique() > 1), "lambda_um")
    others = [c for c in cands if c != x and df[c].nunique() > 1]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    groups = df.groupby(others) if others else [((), df)]
    for key, g in groups:
        for (side, n, m), o in g.groupby(["side", "n", "m"]):
            name = f"{side}:{n}" + (f":{m}" if m else "")
            if a.orders and name not in a.orders and f"{side}:{n}:{m}" not in a.orders:
                continue
            o = o.sort_values(x)
            label = name + ("" if not others else " (" + ", ".join(f"{c}={v:g}" for c, v in
                                                                 zip(others, key if isinstance(key, tuple) else (key,))) + ")")
            ax.plot(o[x], o[a.quantity], marker=".", label=label)
    ax.set_xlabel(x)
    ax.set_ylabel(a.quantity)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    if a.out:
        fig.savefig(a.out, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
