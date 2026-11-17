"""The compact reference for language-model agents (python/fouriergrating/LLM_REFERENCE.md) is complete and runs.

Static checks (no compiled core): every public Python name and class member, every layer type, every JSON key the
parser reads and every CAD finding code appears in the reference. With the core: every ``json`` block is a complete
configuration that solves, every ``python`` block runs, and the text ships with the package.
"""

import json
import pathlib
import re
import subprocess
import sys

import pytest

from test_handbook_reference import PKG, ROOT, class_members, module_all

REF = PKG / "LLM_REFERENCE.md"


def text() -> str:
    return REF.read_text(encoding="utf-8")


def blocks(lang: str) -> list:
    return re.findall(r"```" + lang + r"\n(.*?)```", text(), flags=re.S)


def missing(names) -> list:
    t = text()
    return sorted(n for n in set(names) if n not in t)


def test_python_names():
    names = list(module_all(PKG / "__init__.py")) + ["__version__"]
    for mod in ("cad", "materials", "torch", "jax"):
        names += module_all(PKG / f"{mod}.py")
    for members in class_members(PKG / "__init__.py").values():
        names += members
    assert not missing(names), "add to python/fouriergrating/LLM_REFERENCE.md"


def test_layer_types_keys_and_codes():
    cfg = (ROOT / "src" / "io" / "config.cpp").read_text(encoding="utf-8")
    types = set(re.findall(r'type == "([a-z0-9]+)"', cfg)) - {"asr", "none", "plane_wave", "gaussian_beam"}
    keys = set(re.findall(r'\bj(?:\.value|\.at|\.contains)\("([A-Za-z_0-9]+)"', cfg))
    keys |= set(re.findall(r'\bj\["([A-Za-z_0-9]+)"\]', cfg))
    common = (ROOT / "src" / "cad" / "cad_common.cpp").read_text(encoding="utf-8")
    codes = set(re.findall(r"\b([A-Z]+(?:_[A-Z]+)+)\b", common.split('"description"', 1)[1].split("positions", 1)[0]))
    assert len(types) >= 7 and len(keys) > 50 and len(codes) == 10
    assert not missing(f"`{t}`" for t in types | {"cad"})
    assert not missing(keys)
    assert not missing(codes)


def test_compact():
    assert len(text()) < 20000  # about 5000 tokens: fits next to a task in a model's context


@pytest.fixture
def fg():
    return pytest.importorskip("fouriergrating")


def test_ships_with_the_package(fg):
    assert fg.llm_reference() == text()
    out = subprocess.run([sys.executable, "-m", "fouriergrating"], capture_output=True, text=True, check=True)
    assert out.stdout == text()


@pytest.mark.parametrize("k", range(len(blocks("json"))))
def test_json_blocks_solve(k, fg):
    res = fg.solve(json.loads(blocks("json")[k]))
    assert len(res.reflected) > 0 and -1e-9 < res.absorption < 1


@pytest.mark.parametrize("k", range(len(blocks("python"))))
def test_python_blocks_run(k, fg, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    exec(compile(blocks("python")[k], f"LLM_REFERENCE.md python block {k}", "exec"), {"__name__": "__ref__"})
