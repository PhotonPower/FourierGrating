"""Writes the STEP files of the CAD examples and test data with the Nurbforge Python package (M14).

    pip install nurbforge>=0.6.0
    python scripts/make_cad_examples.py          # examples/cad/*.step, tests/data/cad/*.step

Conventions for FourierGrating (docs/cad.md): models in micrometres, written with unit="um", one named
solid per body (names -> materials in the "bodies" block), CAD +Z = stacking direction (FG y).
"""

from __future__ import annotations

import math
import pathlib

import nurbforge as nf

ROOT = pathlib.Path(__file__).resolve().parents[1]
STAMP = "2026-11-05T00:00:00"  # reproducible files


def disk(cx: float, cy: float, z: float, r: float) -> nf.TrimmedSurface:
    rect = nf.make_rectangle((cx - r, cy - r, z), (2 * r, 0, 0), (0, 2 * r, 0))
    return nf.TrimmedSurface(rect, nf.circle_loop(0.5, 0.5, 0.5))


def cylinder(cx: float, cy: float, r: float, z0: float, z1: float) -> nf.Shell:
    side = nf.extrude(nf.make_circle((cx, cy, z0), (0, 0, 1), r), (0, 0, z1 - z0))
    return nf.Shell.sew([nf.TrimmedSurface(side), disk(cx, cy, z0, r), disk(cx, cy, z1, r)])


def frustum(cx: float, cy: float, r0: float, r1: float, z0: float, h: float) -> nf.Shell:
    side = nf.revolve(nf.make_line((cx + r0, cy, z0), (cx + r1, cy, z0 + h)), (cx, cy, 0), (0, 0, 1), 2 * math.pi)
    return nf.Shell.sew([nf.TrimmedSurface(side), disk(cx, cy, z0, r0), disk(cx, cy, z0 + h, r1)])


def prism(xy: list[tuple[float, float]], z0: float, z1: float) -> nf.Shell:
    faces = []
    for (ax, ay), (bx, by) in zip(xy, xy[1:] + xy[:1]):
        faces.append(nf.TrimmedSurface(nf.extrude(nf.make_line((ax, ay, z0), (bx, by, z0)), (0, 0, z1 - z0))))
    x0, x1 = min(p[0] for p in xy), max(p[0] for p in xy)
    y0, y1 = min(p[1] for p in xy), max(p[1] for p in xy)
    uv = [((x - x0) / (x1 - x0), (y - y0) / (y1 - y0)) for x, y in xy]
    for z in (z0, z1):
        rect = nf.make_rectangle((x0, y0, z), (x1 - x0, 0, 0), (0, y1 - y0, 0))
        faces.append(nf.TrimmedSurface(rect, nf.polygon_loop(uv)))
    return nf.Shell.sew(faces)


def write(path: pathlib.Path, shells: list[nf.Shell], names: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    nf.write_step(path, shells, unit="um", names=names, name=path.stem, timestamp=STAMP)
    print("written", path.relative_to(ROOT))


def main() -> None:
    # example: tapered pillar (frustum) with a cylindrical lid, period 0.5 um (examples/cad_pillar_array.json)
    write(ROOT / "examples" / "cad" / "pillar_array.step",
          [frustum(0.25, 0.25, 0.16, 0.11, 0.0, 0.2), cylinder(0.25, 0.25, 0.11, 0.2, 0.26)], ["Pillar", "Lid"])
    # test data: cylinder (V17) and L prism (V16b)
    write(ROOT / "tests" / "data" / "cad" / "cylinder.step", [cylinder(0.25, 0.25, 0.15, 0.0, 0.2)], ["Pillar"])
    L = [(0.08, 0.06), (0.38, 0.06), (0.38, 0.18), (0.20, 0.18), (0.20, 0.42), (0.08, 0.42)]
    write(ROOT / "tests" / "data" / "cad" / "l_prism.step", [prism(L, 0.0, 0.2)], ["L"])


if __name__ == "__main__":
    main()
