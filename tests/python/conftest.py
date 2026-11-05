import os
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
EXAMPLES = [s.strip() for s in (ROOT / "tests" / "regression" / "examples.txt").read_text().splitlines()
            if s.strip() and not s.startswith("#")]


@pytest.fixture(scope="session")
def root() -> pathlib.Path:
    return ROOT


@pytest.fixture(scope="session")
def fgsolve():
    """Path of an fgsolve executable (env FG_FGSOLVE) for direct comparisons, else skip."""
    exe = os.environ.get("FG_FGSOLVE")
    if not exe or not os.path.isfile(exe):
        pytest.skip("FG_FGSOLVE not set")
    return exe
