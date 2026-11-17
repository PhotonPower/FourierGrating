"""Release metadata (static, no compiled core needed): one version for both distributions.

fouriergrating-cuda12 ships only the compiled core fouriergrating_cuda12._core and uses the Python part of
fouriergrating; it must require exactly the version it was built from (CMakeLists.txt).
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def project_version() -> str:
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    return re.search(r"project\(FourierGrating VERSION ([0-9]+\.[0-9]+\.[0-9]+)", cmake).group(1)


def test_cuda12_requires_the_same_version():
    toml = (ROOT / "python" / "cuda12" / "pyproject.toml").read_text(encoding="utf-8")
    pins = re.findall(r'"fouriergrating==([0-9.]+)"', toml)
    assert pins == [project_version()]


def test_changelog_has_a_section_for_the_version():
    log = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{project_version()}]" in log
