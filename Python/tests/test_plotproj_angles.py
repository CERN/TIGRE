"""Exercise the projection animation callback with public TIGRE angle shapes."""
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pytest
import importlib.util
from pathlib import Path
path = Path(__file__).resolve().parents[1] / "tigre/utilities/visualization/plotproj.py"
spec = importlib.util.spec_from_file_location("plotproj_under_test", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
plotProj = module.plotProj


@pytest.mark.parametrize("angles", [
    [0, np.pi/2, np.pi],
    np.array([0, np.pi/2, np.pi]),
    np.array([[0], [np.pi/2], [np.pi]]),
    np.array([[0, .1, .2], [np.pi/2, .3, .4], [np.pi, .5, .6]]),
])
@pytest.mark.parametrize("dim", [None, "T", "t"])
def test_projection_animation_draws_first_and_last_view(angles, dim, created_animations):
    plot = plotProj(np.ones((3, 4, 6)), angles, dim=dim, show_plot=False)
    animation = created_animations[-1]
    try:
        animation._fig.canvas.draw()  # Triggers the same _init_draw as the GUI.
        animation._func(2, *animation._args)
        animation._fig.canvas.draw()
        assert "T:2, alpha=+1.000 pi" == animation._fig.axes[0].get_title()
    finally:
        animation._draw_was_started = True
        plt.close(animation._fig)


@pytest.mark.parametrize("angles", [np.arange(3), np.zeros((3, 3))])
def test_static_view_uses_scalar_alpha(angles):
    try:
        plotProj(np.ones((3, 4, 6)), angles, slice=1, show_plot=False)
        assert "alpha=" in plt.gca().get_title()
    finally:
        plt.close("all")


@pytest.mark.parametrize("dim,last", [("U", 5), ("V", 3)])
def test_detector_slices_do_not_index_angles_as_detector_coordinates(dim, last, created_animations):
    plot = plotProj(np.ones((3, 4, 6)), np.arange(3), dim=dim, show_plot=False)
    animation = created_animations[-1]
    try:
        animation._func(last, *animation._args)
        assert "alpha=" not in animation._fig.axes[0].get_title()
    finally:
        animation._draw_was_started = True
        plt.close(animation._fig)


@pytest.mark.parametrize("angles", [np.zeros(2), np.zeros((3, 2)), np.zeros((1, 3)),
                                    np.zeros((3, 0)), [0, np.nan, 1]])
def test_invalid_angles_fail_before_animation(angles):
    with pytest.raises(ValueError, match="angles must be finite"):
        plotProj(np.ones((3, 4, 6)), angles, show_plot=False)


def test_empty_gif_filename_keeps_animation_without_saving(created_animations):
    plot = plotProj(np.ones((3, 4, 6)), savegif="", show_plot=False)
    try:
        assert plot.savegif is None
        created_animations[-1]._fig.canvas.draw()
    finally:
        created_animations[-1]._draw_was_started = True
        plt.close("all")


@pytest.fixture
def created_animations(monkeypatch):
    created = []
    original = module.animation.FuncAnimation
    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        created.append(result)
        return result
    monkeypatch.setattr(module.animation, "FuncAnimation", capture)
    yield created
    for animation in created:
        animation._draw_was_started = True
    plt.close("all")
