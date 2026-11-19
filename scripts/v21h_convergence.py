#!/usr/bin/env python3
"""V21h (Plan "Tests und Abnahme (M20)"): convergence study of the cylindrical solver.

Runs fgsolve on generated configurations and writes docs/validation_data/v21h_convergence.csv:

* sphere: dielectric sphere (n = 1.5, k0 a = 2) in air, plane wave at phi = 40 deg (many harmonics), staircase of a
  rings_profile; sigma_sca against Mie (BHMIE below) while one of kappa_max, dk_factor, slices, m_cut changes;
* aperture: hole (r = 0.25 um) in a 0.1 um silver film on glass, Gaussian beam at 633 nm, R and T against the
  finest run: in air (aperture_air, dk_factor 1 ... 6) the plasmon of the silver/air interface leaks into the glass
  and lies between the branch points (real path), so the convergence in dk_factor is slow; with glass on both
  sides (aperture_sym) both plasmons are bound and on the deformed path, and dk_factor 2 is converged.

Usage: python3 scripts/v21h_convergence.py [path/to/fgsolve]   (default build/cpu-release/apps/fgsolve)
Only the standard library; the runs take about 20 minutes (one core, 8 GB suffice).
"""

import csv
import json
import math
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
FGSOLVE = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "build" / "cpu-release" / "apps" / "fgsolve"
OUT = ROOT / "docs" / "validation_data" / "v21h_convergence.csv"


def mie_qsca(x: float, m: float) -> float:
    """Scattering efficiency of a sphere (Bohren-Huffman BHMIE, real relative index m)."""
    nstop = int(x + 4 * x ** (1 / 3) + 2)
    nmx = int(max(nstop, m * x)) + 16
    d = [0.0] * (nmx + 1)
    for n in range(nmx, 0, -1):
        d[n - 1] = n / (m * x) - 1 / (d[n] + n / (m * x))
    psi0, psi1, chi0, chi1 = math.cos(x), math.sin(x), -math.sin(x), math.cos(x)
    xi1 = complex(psi1, -chi1)
    q = 0.0
    for n in range(1, nstop + 1):
        psi = (2 * n - 1) * psi1 / x - psi0
        chi = (2 * n - 1) * chi1 / x - chi0
        xi = complex(psi, -chi)
        da, db = d[n] / m + n / x, m * d[n] + n / x
        an = (da * psi - psi1) / (da * xi - xi1)
        bn = (db * psi - psi1) / (db * xi - xi1)
        q += (2 * n + 1) * (abs(an) ** 2 + abs(bn) ** 2)
        psi0, psi1, chi0, chi1, xi1 = psi1, psi, chi1, chi, xi
    return 2 * q / x**2


def run(cfg: dict, work: pathlib.Path, name: str) -> dict:
    """fgsolve on cfg; returns {side: value} of efficiencies.csv (n = 0 rows)."""
    out = work / name
    cfg = dict(cfg, output={"dir": str(out)})
    path = work / f"{name}.json"
    path.write_text(json.dumps(cfg))
    subprocess.run([str(FGSOLVE), str(path)], check=True, capture_output=True, cwd=work)
    vals = {}
    with open(out / "efficiencies.csv") as f:
        for row in csv.DictReader(f):
            if row["n"] == "0":
                vals[row["side"]] = float(row["efficiency"])
    return vals


LAM, KA, NS = 1.0, 2.0, 1.5
A = KA * LAM / (2 * math.pi)


def sphere(kappa=4.0, dk=1.0, slices=32, mcut=None) -> dict:
    cyl = {"spectrum": {"kappa_max": kappa, "dk_factor": dk}, "factorization": "li_r"}
    if mcut is not None:
        cyl["harmonics"] = list(range(-mcut, mcut + 1))
    return {
        "cylindrical": cyl,
        "materials": {"air": {"n": 1.0}, "sph": {"n": NS}},
        "superstrate": "air",
        "substrate": "air",
        "incidence": {"type": "plane_wave", "lambda_um": LAM, "phi_deg": 40, "psi_deg": 0, "delta_deg": 0},
        "layers": [{"type": "rings_profile", "thickness_um": 2 * A, "slices": slices,
                    "radii": [{"sphere": {"center_y_um": A, "radius_um": A}}], "materials": ["sph", "air"]}],
    }


def aperture(kappa=2.0, dk=1.0, sup="air") -> dict:
    return {
        "cylindrical": {"spectrum": {"kappa_max": kappa, "dk_factor": dk}, "factorization": "li_r"},
        "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "ag": {"n": [0.135, 3.99]}},
        "superstrate": sup,
        "substrate": "glass",
        "incidence": {"type": "gaussian_beam", "lambda_um": 0.633, "waist_um": 0.4, "focus_y_um": 0.1,
                      "polarization": "x"},
        "layers": [{"type": "rings", "thickness_um": 0.1, "radii_um": [0.25], "materials": ["air", "ag"]}],
    }


def main() -> None:
    rows = []
    sigma_mie = mie_qsca(KA, NS) * math.pi * A * A
    m_auto = math.ceil(KA * math.sin(math.radians(40)) + 4 * (KA * math.sin(math.radians(40))) ** (1 / 3) + 2)
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        studies = [("kappa_max", [2, 3, 4, 6, 8], lambda v: sphere(kappa=v)),
                   ("dk_factor", [0.5, 1, 2], lambda v: sphere(dk=v)),
                   ("slices", [8, 16, 32, 64], lambda v: sphere(slices=v)),
                   ("m_cut", list(range(1, m_auto + 2)), lambda v: sphere(mcut=v))]
        for param, values, make in studies:
            for v in values:
                s = run(make(v), work, f"sphere_{param}_{v}")["C_Ssca"]
                rows.append(["sphere", param, v, "sigma_sca_um2", s, sigma_mie, s / sigma_mie - 1])
                print(f"sphere {param} = {v}: sigma = {s:.8f} um^2, {s / sigma_mie - 1:+.2e} vs Mie", flush=True)
        aperture_studies = [("aperture_air", "dk_factor", [1, 2, 3, 4, 5, 6], lambda v: aperture(dk=v)),
                            ("aperture_sym", "dk_factor", [1, 2, 3, 4], lambda v: aperture(dk=v, sup="glass")),
                            ("aperture_sym", "kappa_max", [1.5, 2, 3], lambda v: aperture(kappa=v, dk=2, sup="glass"))]
        for case, param, values, make in aperture_studies:
            res = [(v, run(make(v), work, f"{case}_{param}_{v}")) for v in values]
            ref = res[-1][1]  # the finest run
            for v, r in res:
                for q, side in [("R", "C_R"), ("T", "C_T")]:
                    rows.append([case, param, v, q, r[side], ref[side], r[side] - ref[side]])
                print(f"{case} {param} = {v}: R = {r['C_R']:.10f} T = {r['C_T']:.10f}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "parameter", "value", "quantity", "result", "reference", "error"])
        w.writerows(rows)
    print(f"written: {OUT} (Mie sigma_sca = {sigma_mie:.8f} um^2, m_cut auto = {m_auto})")


if __name__ == "__main__":
    main()
