"""Execute the demo's control flow with deterministic CPU-only operators."""
import runpy
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pytest

@pytest.fixture
def demo(monkeypatch):
    head = np.zeros((2, 3, 4), dtype=np.float32)
    projections = np.ones((3, 4, 5), dtype=np.float32)
    noisy = projections + 5
    calls = []
    histories = {}
    tigre = ModuleType("tigre")
    tigre.__path__ = []
    tigre.geometry_default = lambda **kwargs: SimpleNamespace(nVoxel=np.array(head.shape))
    tigre.Ax = lambda *args: projections
    utilities = ModuleType("tigre.utilities")
    utilities.__path__ = []
    utilities.sample_loader = SimpleNamespace(load_head_phantom=lambda shape: head)
    utilities.CTnoise = SimpleNamespace(add=lambda *args, **kwargs: noisy)
    algorithms = ModuleType("tigre.algorithms")
    def reconstruct(name):
        def run(data, geo, angles, niter, **kwargs):
            calls.append((name, data))
            quality = np.arange(4*niter).reshape(4, niter) + 1000*len(calls)
            histories[name] = quality
            return head.copy(), quality
        return run
    for name in ("sirt", "sart", "ossart"):
        setattr(algorithms, name, reconstruct(name))
    metrics = ModuleType("tigre.utilities.Measure_Quality")
    metrics.Measure_Quality = lambda *args: [0.]
    for name, module in (("tigre", tigre), ("tigre.utilities", utilities),
                         ("tigre.algorithms", algorithms),
                         ("tigre.utilities.Measure_Quality", metrics)):
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(plt, "show", lambda: None)
    plt.close("all")
    path = Path(__file__).resolve().parents[1] / "demos/d10_MeasureQuality.py"
    runpy.run_path(str(path))
    yield calls, noisy, histories
    plt.close("all")

def test_reconstructions_use_the_generated_noisy_projections(demo):
    calls, noisy, histories = demo
    assert len(calls) == 3
    assert all(data is noisy for name, data in calls)

def test_each_metric_plot_uses_its_own_history_row(demo):
    calls, noisy, histories = demo
    assert len(plt.get_fignums()) == 4
    for row, number in enumerate(plt.get_fignums()):
        ax = plt.figure(number).axes[0]
        assert len(ax.lines) == 3
        for line, name in zip(ax.lines, ("sirt", "ossart", "sart")):
            np.testing.assert_array_equal(line.get_ydata(), histories[name][row])
