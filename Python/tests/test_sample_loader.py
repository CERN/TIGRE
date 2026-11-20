"""CPU-only regression checks for bundled head phantom loading."""
import importlib.util
import site
import sys
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


def _installed_sample_loader(site_packages):
    # Import a copy of sample_loader.py from where pip installs it: site_packages/tigre/utilities.
    target = site_packages / "tigre" / "utilities" / "sample_loader.py"
    target.parent.mkdir(parents=True)
    target.write_bytes(path.read_bytes())
    spec = importlib.util.spec_from_file_location("sample_loader_installed", target)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("site_packages", ["lib/python3.12/site-packages", "Lib/site-packages"])
def test_head_phantom_is_found_where_pip_installs_data_files(tmp_path, monkeypatch, site_packages):
    # setup.py installs head.mat as data_files=[("data", ...)], which pip puts in <prefix>/data
    # (Linux and Windows layouts of site-packages).
    prefix = tmp_path / "env"
    loader = _installed_sample_loader(prefix / site_packages)
    (prefix / "data").mkdir()
    loader.scipy.io.savemat(prefix / "data" / "head.mat", {"img": np.ones((4, 5, 6))})
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.setattr(site, "getuserbase", lambda: str(tmp_path / "user"))
    assert loader.load_head_phantom([6, 5, 4]).shape == (6, 5, 4)


def test_head_phantom_is_found_after_pip_install_user(tmp_path, monkeypatch):
    user_base = tmp_path / "user"
    loader = _installed_sample_loader(user_base / "lib/python3.12/site-packages")
    (user_base / "data").mkdir()
    loader.scipy.io.savemat(user_base / "data" / "head.mat", {"img": np.ones((4, 5, 6))})
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "env"))
    monkeypatch.setattr(site, "getuserbase", lambda: str(user_base))
    assert loader.load_head_phantom([6, 5, 4]).shape == (6, 5, 4)
