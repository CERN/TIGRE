"""The supported interface for external geometry-design applications."""

import subprocess
import sys

import pytest

from tigre import design_api


def test_every_exported_name_exists_and_is_unique():
    assert len(design_api.__all__) == len(set(design_api.__all__))
    for name in design_api.__all__:
        assert hasattr(design_api, name), name


def test_api_version_and_requirement_check():
    major, minor = design_api.API_VERSION
    design_api.require_api(major, minor)
    design_api.require_api(major)
    with pytest.raises(RuntimeError, match="required"):
        design_api.require_api(major, minor + 1)
    with pytest.raises(RuntimeError, match="required"):
        design_api.require_api(major + 1)


def test_facade_does_not_import_a_desktop_toolkit():
    script = (
        "import sys; import tigre.design_api; "
        "bad = [m for m in sys.modules if m == 'tkinter' or m.startswith('tigre_geometry_design_tool')]; "
        "assert not bad, bad"
    )
    subprocess.run([sys.executable, "-c", script], check=True)
