"""Each compiled TIGRE module can be imported before `tigre` itself.

The Cython modules used to import tigre.utilities.errors at load time, which
imports the whole `tigre` package, which in turn imports the compiled modules
back. Importing e.g. `_minTV` first then failed with "cannot import name
'minTV' from partially initialized module '_minTV' (most likely due to a
circular import)".

Each case runs in a fresh interpreter, because the import order only matters
the first time a module is loaded. No GPU is needed.
"""
import subprocess
import sys

import pytest

COMPILED_MODULES = [
    "_Ax",
    "_Atb",
    "_tv_proximal",
    "_minTV",
    "_AwminTV",
    "_minPICCS",
    "_gpuUtils",
    "_RandomNumberGenerator",
]


def _compiled_modules_available():
    try:
        import tigre  # noqa: F401
        import _Ax  # noqa: F401
    except Exception:
        return False
    return True


@pytest.mark.skipif(not _compiled_modules_available(), reason="TIGRE's compiled modules are not built")
@pytest.mark.parametrize("module", COMPILED_MODULES)
def test_compiled_module_can_be_imported_before_tigre(module):
    code = f"import {module}; import tigre; import tigre.algorithms"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
