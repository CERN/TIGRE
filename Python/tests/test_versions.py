"""tigre.__version__ and tigre.show_versions() (CPU only)."""
import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "Python/tigre/utilities/versions.py"
spec = importlib.util.spec_from_file_location("tigre_versions_under_test", path)
versions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(versions)


def _pyproject_version():
    text = (ROOT / "pyproject.toml").read_text()
    return re.search(r'^version\s*=\s*"([^"]+)"', text, re.M).group(1)


def test_installed_version_matches_pyproject():
    try:
        from importlib.metadata import version

        installed = version("tigre")
    except Exception:
        pytest.skip("TIGRE is not installed")
    assert installed == _pyproject_version()
    assert versions._tigre_version() == installed


def test_version_is_unknown_when_not_installed(monkeypatch):
    import importlib.metadata as md

    def missing(name):
        raise md.PackageNotFoundError(name)

    monkeypatch.setattr(md, "version", missing)
    assert versions._tigre_version() == "unknown"


def test_show_versions_survives_a_broken_gpu_query(monkeypatch, capsys):
    def broken():
        raise RuntimeError("no CUDA driver")

    monkeypatch.setattr(versions, "_gpu_names", broken)
    versions.show_versions()
    out = capsys.readouterr().out
    assert out.splitlines()[0].startswith("TIGRE")
    assert "Python" in out and "NumPy" in out
    assert "could not query them (no CUDA driver)" in out


def test_show_versions_lists_gpus(monkeypatch, capsys):
    monkeypatch.setattr(versions, "_gpu_names", lambda: ["GPU A", "GPU B"])
    versions.show_versions()
    assert "GPUs     GPU A, GPU B" in capsys.readouterr().out
