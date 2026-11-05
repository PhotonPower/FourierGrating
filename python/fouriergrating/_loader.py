"""Selection of the compiled core (Plan section "M12: Build").

The CUDA distribution ``fouriergrating-cuda12`` ships only ``fouriergrating_cuda12._core`` (same API,
built with the CUDA backend). It is used when it is installed and a GPU is usable; otherwise the CPU
core ``fouriergrating._core`` is loaded. CUDA DLLs are never bundled: on Windows the directories of the
pip packages ``nvidia-*-cu12`` (and, for builds from source, ``CUDA_PATH``) are registered with
``os.add_dll_directory`` before the import.
"""

from __future__ import annotations

import glob
import os
import sys
from types import ModuleType


def _add_cuda_dll_directories() -> None:
    if sys.platform != "win32":
        return
    dirs = []
    for base in sys.path:
        dirs += glob.glob(os.path.join(base, "nvidia", "*", "bin"))
    cuda = os.environ.get("CUDA_PATH")
    if cuda:
        dirs += [os.path.join(cuda, "bin"), os.path.join(cuda, "bin", "x64")]
    for d in dirs:
        if os.path.isdir(d):
            try:
                os.add_dll_directory(d)
            except OSError:
                pass


def load_core() -> ModuleType:
    if os.environ.get("FG_CORE", "auto") != "cpu":
        _add_cuda_dll_directories()
        try:
            from fouriergrating_cuda12 import _core as cuda_core  # type: ignore[import-not-found]

            if cuda_core.cuda_available() or os.environ.get("FG_CORE") == "cuda12":
                return cuda_core
        except ImportError:
            pass
    from . import _core

    return _core
