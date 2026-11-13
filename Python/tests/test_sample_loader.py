"""CPU-only regression checks for bundled head phantom loading."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

path = Path(__file__).resolve().parents[1] / "tigre/utilities/sample_loader.py"
spec = importlib.util.spec_from_file_location("sample_loader_under_test", path)
sample_loader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_loader)

def test_missing_head_phantom_has_an_actionable_error(monkeypatch):
    monkeypatch.setattr(sample_loader.os.path, "isfile", lambda path: False)
    with pytest.raises(FileNotFoundError, match="head.mat"):
        sample_loader.load_head_phantom([8, 10, 12])

def test_head_phantom_resampling_preserves_requested_shape_and_axis_order(monkeypatch):
    image = np.arange(5*6*7, dtype=float).reshape(5, 6, 7)
    monkeypatch.setattr(sample_loader.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(sample_loader.scipy.io, "loadmat", lambda path: {"img": image})
    head = sample_loader.load_head_phantom([7, 6, 5])
    np.testing.assert_allclose(head, sample_loader.scipy.ndimage.zoom(
        image.transpose(2, 1, 0), (1, 1, 1), order=3, prefilter=False))
    resized = sample_loader.load_head_phantom([8, 10, 12])
    assert resized.shape == (8, 10, 12)
    assert np.isfinite(resized).all()
