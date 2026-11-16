"""The command reference of the handbook (docs/handbuch/anh_befehle.tex) lists every public name.

Static checks of the sources (no compiled core needed): public Python names and class members, the options
of fgsolve and fgbench, the environment variables read by the core and the Python package, the CMake options
and the scripts. A name that is added to the code but not to the reference makes this test fail.
"""

import ast
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
PKG = ROOT / "python" / "fouriergrating"


def reference_text() -> str:
    tex = (ROOT / "docs" / "handbuch" / "anh_befehle.tex").read_text(encoding="utf-8")
    return tex.replace("\\_", "_").replace("-{}-", "--")


def module_all(path: pathlib.Path) -> list:
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "__all__" for t in node.targets):
            return [ast.literal_eval(e) for e in node.value.elts]
    return []


def class_members(path: pathlib.Path) -> dict:
    """Public methods, properties, dataclass fields and attributes assigned to self, per public class."""
    out = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if not isinstance(node, ast.ClassDef) or node.name.startswith("_"):
            continue
        names = set()
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                if not item.name.startswith("_"):
                    names.add(item.name)
                for sub in ast.walk(item):
                    if (isinstance(sub, ast.Attribute) and isinstance(sub.ctx, ast.Store)
                            and getattr(sub.value, "id", "") == "self" and not sub.attr.startswith("_")):
                        names.add(sub.attr)
            elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
                names.add(item.target.id)
        out[node.name] = names
    return out


def missing(names) -> list:
    text = reference_text()
    return sorted(n for n in set(names) if n not in text)


def test_python_names():
    names = list(module_all(PKG / "__init__.py")) + ["__version__"]
    for mod in ("cad", "materials", "torch", "jax"):
        names += module_all(PKG / f"{mod}.py")
    for members in class_members(PKG / "__init__.py").values():
        names += members
    assert not missing(names), "add to docs/handbuch/anh_befehle.tex"


def test_command_line_options():
    opts = []
    for app in ("fgsolve", "fgbench"):
        opts += re.findall(r'"(--[a-z-]+)"', (ROOT / "apps" / app / "main.cpp").read_text(encoding="utf-8"))
    assert {"--profile", "--sizes"} <= set(opts)
    assert not missing(opts)


def test_environment_variables():
    env = []
    for f in list((ROOT / "src").rglob("*.cpp")) + list((ROOT / "tests").rglob("*.cpp")):
        env += re.findall(r'getenv\("([A-Z0-9_]+)"\)', f.read_text(encoding="utf-8", errors="replace"))
    env += re.findall(r"FG_REGRESSION_[A-Z]+", (ROOT / "tests" / "regression" / "test_regression.cpp").read_text())
    for f in PKG.glob("*.py"):
        env += re.findall(r'environ\.get\("([A-Z0-9_]+)"', f.read_text(encoding="utf-8"))
    assert "FG_NUM_THREADS" in env and "FG_CORE" in env
    assert not missing(env)


def test_cmake_options_and_scripts():
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    names = re.findall(r"option\((FG_[A-Z0-9_]+)", cmake) + re.findall(r"set\((FG_[A-Z0-9_]+) .*CACHE", cmake)
    names += [p.name for p in (ROOT / "scripts").iterdir() if p.is_file()]
    assert "FG_ENABLE_NURBFORGE" in names
    assert not missing(names)
