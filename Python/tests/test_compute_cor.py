"""computeCOR: does it run, does it recover a known centre of rotation.

The reference is MATLAB/Utilities/computeCOR.m; see the module docstring of
tigre/utilities/computeCOR.py for the line-by-line comparison and the one
deliberate divergence (score normalised by the count of valid pixels, which
is what MATLAB's `length(find(test > 0))` computes).
"""
import numpy as np
import pytest

import tigre

TOL_MM = 0.2


def _problem(cor_mm, nang=120):
    geo = tigre.geometry(mode="cone", nVoxel=np.array([128, 128, 128]), default=True)
    angles = np.linspace(0, 2 * np.pi, nang, endpoint=False, dtype=np.float32)
    geo.check_geo(angles)
    geo.COR = np.ones(nang, dtype=np.float32) * np.float32(cor_mm)
    z, y, x = np.mgrid[:128, :128, :128].astype(np.float32)
    ph = np.zeros((128, 128, 128), dtype=np.float32)
    ph[((x - 70) ** 2 + (y - 60) ** 2 + (z - 64) ** 2) < 26 ** 2] = 1.0
    ph[((x - 50) ** 2 + (y - 74) ** 2 + (z - 64) ** 2) < 11 ** 2] = 2.0
    return geo, angles, tigre.Ax(ph, geo, angles, "Siddon")


@pytest.mark.parametrize("true_cor", [0.0, 0.5, -0.5, 1.5, -2.0])
def test_recovers_a_known_cor(true_cor):
    geo, angles, proj = _problem(true_cor)
    est = float(tigre.computeCOR(proj, geo, angles))
    assert abs(est - true_cor) < TOL_MM, f"true {true_cor}, estimated {est:.4f}"


def test_descending_angles_give_the_same_answer():
    geo, angles, proj = _problem(0.75)
    a = float(tigre.computeCOR(proj, geo, angles))
    b = float(tigre.computeCOR(proj[::-1], geo, angles[::-1]))
    assert a == pytest.approx(b, abs=1e-6)
