#!/usr/bin/env python3
"""Plot near-field cuts written by FourierGrating (Plan, section "Nahfeld").

Each cut directory contains meta.json and one .npy file per quantity, e.g.
    nearfield_vertical/meta.json, Ex.npy, Ey.npy, ..., absE2.npy, material.npy

meta.json (written by fg::NearFieldWriter):
    {"cut": "vertical" | "horizontal", "dimension": "1D" | "2D",
     "axes": {"axis1": {"name": "x", "values_um": [...]},
              "axis2": {"name": "y", "values_um": [...]}},     # axis2 missing for a 1D line
     "lambda_um": 0.6328, "incidence": {...}, "quantities": ["Ex", ...]}
Arrays are C-ordered with shape (len(axis2), len(axis1)) for maps, (len(axis1),) for lines.

Usage:
    plot_nearfield.py <cut_dir> [--quantity absE2] [--part abs|re|im|phase] [--out file.png]
    plot_nearfield.py --demo <dir>      # writes synthetic test data and plots it
"""
import argparse
import json
import os
import sys

import numpy as np


def load_cut(path):
    with open(os.path.join(path, "meta.json")) as f:
        meta = json.load(f)
    data = {}
    for q in meta.get("quantities", []):
        fn = os.path.join(path, q + ".npy")
        if os.path.exists(fn):
            data[q] = np.load(fn)
    mat = os.path.join(path, "material.npy")
    if os.path.exists(mat):
        data["material"] = np.load(mat)
    return meta, data


def select_part(arr, part):
    if np.isrealobj(arr) or part == "re":
        return np.real(arr)
    return {"abs": np.abs, "im": np.imag, "phase": np.angle}[part](arr)


def plot_cut(meta, data, quantity, part, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if quantity not in data:
        sys.exit(f"quantity '{quantity}' not found; available: {sorted(data)}")
    a1 = meta["axes"]["axis1"]
    x = np.asarray(a1["values_um"])
    val = select_part(data[quantity], part)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    title = f"{quantity} ({part}), {meta.get('cut', '')} cut, λ = {meta.get('lambda_um', '?')} µm"
    if "axis2" in meta["axes"]:
        a2 = meta["axes"]["axis2"]
        y = np.asarray(a2["values_um"])
        im = ax.pcolormesh(x, y, val, shading="auto", cmap="inferno" if part == "abs" else "RdBu_r")
        fig.colorbar(im, ax=ax)
        if "material" in data:
            ax.contour(x, y, data["material"], levels=np.unique(data["material"])[:-1] + 0.5,
                       colors="w", linewidths=0.6)
        ax.set_ylabel(f"{a2['name']} (µm)")
        ax.set_aspect("equal")
    else:
        ax.plot(x, val)
        ax.set_ylabel(quantity)
    ax.set_xlabel(f"{a1['name']} (µm)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("written:", out)


def write_demo(path):
    """Synthetic vertical cut: standing wave above a lamellar region (only for testing the plot)."""
    os.makedirs(path, exist_ok=True)
    x = np.linspace(0, 1.0, 200)
    y = np.linspace(-0.3, 0.6, 180)
    X, Y = np.meshgrid(x, y)
    mat = ((Y > 0) & (Y < 0.25) & ((X % 0.5) < 0.2)).astype(float)
    Ez = np.exp(1j * 2 * np.pi * Y / 0.6328) + 0.5 * np.exp(-1j * 2 * np.pi * Y / 0.6328) * (1 + 0.3 * mat)
    np.save(os.path.join(path, "Ez.npy"), Ez)
    np.save(os.path.join(path, "absE2.npy"), np.abs(Ez) ** 2)
    np.save(os.path.join(path, "material.npy"), mat)
    meta = {"cut": "vertical", "dimension": "1D", "lambda_um": 0.6328,
            "axes": {"axis1": {"name": "x", "values_um": x.tolist()},
                     "axis2": {"name": "y", "values_um": y.tolist()}},
            "quantities": ["Ez", "absE2"]}
    with open(os.path.join(path, "meta.json"), "w") as f:
        json.dump(meta, f)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cut_dir", nargs="?")
    p.add_argument("--quantity", default="absE2")
    p.add_argument("--part", default="abs", choices=["abs", "re", "im", "phase"])
    p.add_argument("--out", default=None)
    p.add_argument("--demo", metavar="DIR")
    a = p.parse_args()
    if a.demo:
        a.cut_dir = write_demo(a.demo)
    if not a.cut_dir:
        p.error("cut_dir or --demo required")
    meta, data = load_cut(a.cut_dir)
    out = a.out or os.path.join(a.cut_dir, f"{a.quantity}_{a.part}.png")
    plot_cut(meta, data, a.quantity, a.part, out)


if __name__ == "__main__":
    main()
