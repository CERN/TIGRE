"""The same input must give the same result, every call (#860).

`Ax(img, geo, angles, "interpolated")` used to return, about once in 1000-2000
calls, a result with a band of detector rows (almost) zero: the first kernel
read the image's texture before the copy into it had finished. Siddon on cone
beam waited for that copy; interpolated (cone and parallel) and parallel
Siddon did not. The reproduction below is the one from the issue: a uniform
cube, 98^3, a few angles, many calls.

The race was rare, so a passing run is evidence, not proof: with the defaults
(2000 calls per case) a 0.1 % failure rate is caught about 86 % of the time.
Set TIGRE_REPEAT_CALLS higher for a more thorough run.
"""
import os

import numpy as np
import pytest

import tigre

pytestmark = [pytest.mark.gpu, pytest.mark.extended]  # extended: not run by a plain `pytest`, see contributing.md

CALLS = int(os.environ.get("TIGRE_REPEAT_CALLS", "2000"))


def cube_problem(mode):
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([98, 98, 98]))
    img = np.zeros((98, 98, 98), np.float32)
    img[20:78, 20:78, 20:78] = 1.0
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    return geo, img, angles


def assert_repeatable(call, calls=CALLS):
    reference = call()
    bad = []
    for k in range(calls):
        result = call()
        if not np.array_equal(result, reference):
            rows = np.flatnonzero((result != reference).any(axis=(0, 2)))
            bad.append((k, rows.min(), rows.max()))
    assert not bad, "{} of {} calls differed from the first; first: call {}, detector rows {}-{}".format(
        len(bad), calls, bad[0][0], bad[0][1], bad[0][2])


@pytest.mark.parametrize("mode", ["cone", "parallel"])
@pytest.mark.parametrize("projection_type", ["interpolated", "Siddon"])
def test_ax_gives_the_same_result_every_call(mode, projection_type):
    geo, img, angles = cube_problem(mode)
    assert_repeatable(lambda: tigre.Ax(img, geo, angles, projection_type))


@pytest.mark.parametrize("mode", ["cone", "parallel"])
def test_ax_with_one_angle_gives_the_same_result_every_call(mode):
    """One angle per call, as in the issue: the kernel starts right after the
    upload, so the race window is widest."""
    geo, img, _ = cube_problem(mode)
    angle = np.array([0.3])
    assert_repeatable(lambda: tigre.Ax(img, geo, angle, "interpolated"))


@pytest.mark.parametrize("backprojection_type", ["FDK", "matched"])
def test_atb_gives_the_same_result_every_call(backprojection_type):
    geo, img, angles = cube_problem("cone")
    proj = tigre.Ax(img, geo, angles)
    assert_repeatable(lambda: tigre.Atb(proj, geo, angles, backprojection_type), calls=CALLS // 4)
