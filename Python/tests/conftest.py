"""Shared pytest configuration for the TIGRE Python test suite.

Tests that need a CUDA device are marked with ``@pytest.mark.gpu`` (or a
module-level ``pytestmark = pytest.mark.gpu``). They are skipped automatically
when no GPU is visible, so the rest of the suite can run on CPU-only CI runners.
"""
import importlib.util
import os
import sys

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

def _gpu_count():
    try:
        from tigre.utilities.gpu import getGpuNames

        return len(getGpuNames())
    except Exception:  # extension not built, no driver, etc.
        return 0


HAS_GPU = _gpu_count() > 0


def _not_built():
    # The SYCL build does not have every extension yet (SYCL_PENDING in setup.py).
    try:
        from tigre.utilities.cuda_interface import _gpuUtils
    except ImportError:
        return set()
    if getattr(_gpuUtils, "getBackend", lambda: "cuda")() != "sycl":
        return set()
    names = ("_Atb", "_tv_proximal", "_minTV", "_AwminTV", "_minPICCS", "_RandomNumberGenerator")
    names = ["tigre.utilities.cuda_interface." + name for name in names]
    return {name for name in names if importlib.util.find_spec(name) is None}


NOT_BUILT = _not_built()
# Importing them then fails at once (and no test file puts a mock in their place).
sys.modules.update(dict.fromkeys(NOT_BUILT))


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    # On that build, a test that needs a missing extension is skipped, not failed.
    outcome = yield
    report = outcome.get_result()
    error = call.excinfo.value if call.excinfo is not None else None
    if report.failed and isinstance(error, ImportError) and error.name in NOT_BUILT:
        report.outcome = "skipped"
        reason = "Skipped: {} is not in the SYCL build yet".format(error.name)
        report.longrepr = (str(item.path), item.location[1] or 0, reason)


def pytest_collection_modifyitems(config, items):
    for item in items:
        # Tests parametrized with the name of an extension that is not built.
        callspec = getattr(item, "callspec", None)
        missing = NOT_BUILT.intersection(map(str, callspec.params.values())) if callspec else ()
        if missing:
            reason = "{} is not in the SYCL build yet".format(", ".join(sorted(missing)))
            item.add_marker(pytest.mark.skip(reason=reason))
    if HAS_GPU:
        return
    skip_gpu = pytest.mark.skip(reason="no CUDA GPU available")
    for item in items:
        if "gpu" in item.keywords:
            item.add_marker(skip_gpu)
