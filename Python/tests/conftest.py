"""Shared pytest configuration for the TIGRE Python test suite.

Tests that need a CUDA device are marked with ``@pytest.mark.gpu`` (or a
module-level ``pytestmark = pytest.mark.gpu``). They are skipped automatically
when no GPU is visible, so the rest of the suite can run on CPU-only CI runners.
"""
import os

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

# Helper modules and scripts in this folder that are not pytest test files:
# - algorithm_test.py / test_config_gen.py are driven by runscript.sh and need
#   generated configuration files (targetdir.npy, configuration*.npy).
# - stdout_test.py is a legacy Python 2 script.
collect_ignore = ["algorithm_test.py", "stdout_test.py", "test_config_gen.py"]


def _gpu_count():
    try:
        from tigre.utilities.gpu import getGpuNames

        return len(getGpuNames())
    except Exception:  # extension not built, no driver, etc.
        return 0


HAS_GPU = _gpu_count() > 0


def pytest_collection_modifyitems(config, items):
    if HAS_GPU:
        return
    skip_gpu = pytest.mark.skip(reason="no CUDA GPU available")
    for item in items:
        if "gpu" in item.keywords:
            item.add_marker(skip_gpu)

"""Matplotlib compatibility helpers shared by TIGRE's tests."""
import numpy as np


def poly_vertices(collection):
    """Return a Poly3DCollection's vertices as a (coordinates, N) array.

    Matplotlib 3.10 and earlier keep them in ``_vec``; 3.11 keeps the faces in ``_faces``.
    """
    if hasattr(collection, "_vec"):
        return collection._vec
    return np.asarray(collection._faces).reshape(-1, 3).T
