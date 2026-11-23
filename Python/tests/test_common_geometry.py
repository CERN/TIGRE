"""The geometry helpers in tigre.utilities.common_geometry build the scans they describe.

Each test reads the source and detector of every view back with from_TIGRE_to_vector (checked against the projector in
test_vector_geometry.py) and compares them with what the helper's documentation promises. Before these helpers were
rebuilt on from_vector_to_TIGRE, the "fixed" detector of staticDetLinearSourceGeo and ArbitrarySourceDetMoveGeo moved by
several mm, ArbitrarySourceDetMoveGeo was hundreds of mm off for sources around the object, staticDetLinearSourceGeo
failed on the default geometry, and ArbitrarySourceDetectorFixedObject scaled detector offsets by 2*DSO/DSD, misplaced
the source of offset detectors and failed for a source on the -x axis.
"""
import numpy as np
import pytest

import tigre
from tigre.utilities.common_geometry import (
    ArbitrarySourceDetectorFixedObject,
    ArbitrarySourceDetMoveGeo,
    rotation_from_vecs,
    staticDetectorGeo,
    staticDetLinearSourceGeo,
)

ATOL = 1e-9


def _geo(**fields):
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([32, 32, 32]))
    geo.DSO, geo.DSD = 300.0, 500.0
    for name, value in fields.items():
        setattr(geo, name, value)
    return geo


def _vectors(geo):
    return tigre.from_TIGRE_to_vector(geo, geo.angles)


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _central(geo):
    """Vectors of the central position (angle 0) of ``geo``."""
    return tigre.from_TIGRE_to_vector(geo, np.zeros((1, 3)))[0]


# --------------------------------------------------------------------------- staticDetectorGeo (unchanged)


def test_static_detector_arc():
    a = np.linspace(-0.4, 0.4, 9)
    vec = _vectors(staticDetectorGeo(_geo(), a))
    np.testing.assert_allclose(vec[:, 0:3], np.stack([300 * np.cos(a), 300 * np.sin(a), 0 * a], 1), atol=ATOL)
    np.testing.assert_allclose(vec[:, 3:12], np.tile(_central(_geo()), (9, 1))[:, 3:12], atol=ATOL)


# --------------------------------------------------------------------------- staticDetLinearSourceGeo


@pytest.mark.parametrize("s_rot", [0.0, 20.0])
def test_linear_source_fixed_detector_on_default_geometry(s_rot):
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([32, 32, 32]))  # integer rotDetector
    s = np.linspace(-120, 120, 9)
    vec = _vectors(staticDetLinearSourceGeo(geo, s, s_rot=s_rot))
    sr = np.radians(s_rot)
    np.testing.assert_allclose(
        vec[:, 0:3], np.stack([geo.DSO + s * np.sin(sr), s * np.cos(sr), 0 * s], 1), atol=ATOL
    )
    np.testing.assert_allclose(vec[:, 3:12], np.tile(_central(geo)[3:12], (9, 1)), atol=ATOL)


def test_linear_source_keeps_detector_offset_and_rotation_and_rotates():
    geo = _geo(offDetector=np.array([3.0, -5.0]), rotDetector=np.array([0.05, -0.1, 0.2]))
    s = np.linspace(-80, 80, 7)
    rot = np.linspace(0, 90, 7)
    vec = _vectors(staticDetLinearSourceGeo(geo, s, rot=rot))
    central = _central(geo)
    for i in range(7):
        r = _rz(np.radians(rot[i]))
        np.testing.assert_allclose(vec[i, 0:3], r @ (central[0:3] + [0, s[i], 0]), atol=ATOL)
        for k in (3, 6, 9):
            np.testing.assert_allclose(vec[i, k : k + 3], r @ central[k : k + 3], atol=ATOL)


def test_linear_source_refuses_mismatched_lengths():
    with pytest.raises(ValueError):
        staticDetLinearSourceGeo(_geo(), np.zeros(5), rot=np.zeros(4))


# --------------------------------------------------------------------------- ArbitrarySourceDetMoveGeo


def test_arbitrary_sources_default_detector_stays_put():
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([32, 32, 32]))
    src = np.stack([np.full(9, 300.0), np.linspace(-120, 120, 9), np.linspace(-20, 20, 9)], 1)
    vec = _vectors(ArbitrarySourceDetMoveGeo(geo, src))
    np.testing.assert_allclose(vec[:, 0:3], src, atol=ATOL)
    np.testing.assert_allclose(vec[:, 3:12], np.tile(_central(geo)[3:12], (9, 1)), atol=ATOL)


def test_arbitrary_sources_around_object_detector_faces_origin():
    a = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    src = np.stack([300 * np.cos(a), 300 * np.sin(a), np.linspace(-30, 30, 8)], 1)
    det = -src * (200 / 300)
    geo = _geo(offOrigin=np.array([2.0, -3.0, 4.0]))  # (z, y, x)
    vec = _vectors(ArbitrarySourceDetMoveGeo(geo, src, d_pos=det))
    image_centre = geo.offOrigin[::-1]
    np.testing.assert_allclose(vec[:, 0:3], src - image_centre, atol=ATOL)
    np.testing.assert_allclose(vec[:, 3:6], det - image_centre, atol=ATOL)
    normal = np.cross(vec[:, 6:9], vec[:, 9:12])
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    np.testing.assert_allclose(normal, -det / np.linalg.norm(det, axis=1, keepdims=True), atol=ATOL)
    assert np.all(np.abs(vec[:, 8]) < ATOL)  # rows stay horizontal


def test_arbitrary_detector_rotation_matches_rotDetector():
    d_rot = np.array([[4.0, -6.0, 10.0]])
    vec = _vectors(ArbitrarySourceDetMoveGeo(_geo(), np.array([[300.0, 0, 0]]), d_rot=d_rot))
    expected = _central(_geo(rotDetector=np.radians(d_rot[0])))
    np.testing.assert_allclose(vec[0], expected, atol=ATOL)


@pytest.mark.parametrize(
    "kwargs",
    [dict(s_pos=np.zeros((3, 2))), dict(s_pos=np.ones((3, 3)), d_pos=np.ones((2, 3))), dict(s_pos=np.ones((1, 3)), d_pos=np.zeros((1, 3)))],
    ids=["s_pos shape", "d_pos shape", "detector at origin"],
)
def test_arbitrary_refuses_bad_input(kwargs):
    with pytest.raises(ValueError):
        ArbitrarySourceDetMoveGeo(_geo(), **kwargs)


# --------------------------------------------------------------------------- ArbitrarySourceDetectorFixedObject


@pytest.mark.parametrize("dsd", [400.0, 500.0, 600.0, 900.0])
def test_fixed_object_reproduces_vectors(dsd):
    rng = np.random.default_rng(int(dsd))
    angles = np.concatenate([[np.pi, 0.0], rng.uniform(0, 2 * np.pi, 6)])  # pi: source on the -x axis
    geo = _geo(DSD=dsd, offDetector=np.array([4.0, -6.0]), rotDetector=np.array([0.05, -0.1, 0.2]))
    vec = tigre.from_TIGRE_to_vector(geo, angles)
    du, dv = geo.dDetector[1], geo.dDetector[0]
    target = _geo(DSD=dsd)
    out = ArbitrarySourceDetectorFixedObject(
        target, vec[:, 0:3], vec[:, 3:6], vec[:, 6:9] / du, vec[:, 9:12] / dv, use_center_correction=False
    )
    assert out is target
    np.testing.assert_allclose(_vectors(out), vec, atol=1e-8)


def test_fixed_object_centre_correction_and_origin():
    angles = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    vec = tigre.from_TIGRE_to_vector(_geo(), angles)
    shift = np.array([25.0, -10.0, 7.0])
    du, dv = _geo().dDetector[1], _geo().dDetector[0]
    args = (vec[:, 0:3] + shift, vec[:, 3:6] + shift, vec[:, 6:9] / du, vec[:, 9:12] / dv)
    np.testing.assert_allclose(_vectors(ArbitrarySourceDetectorFixedObject(_geo(), *args)), vec, atol=1e-6)
    out = ArbitrarySourceDetectorFixedObject(_geo(), *args, origin_mm=shift, use_center_correction=False)
    np.testing.assert_allclose(_vectors(out), vec, atol=1e-8)


# --------------------------------------------------------------------------- rotation_from_vecs


@pytest.mark.parametrize("v", [[1.0, 0, 0], [0, 2.0, 0], [0.3, -0.4, 0.5]])
def test_rotation_from_opposite_vectors_is_a_rotation(v):
    v = np.array(v)
    r = rotation_from_vecs(v, -v)
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(r), 1.0)
    np.testing.assert_allclose(r @ (v / np.linalg.norm(v)), -v / np.linalg.norm(v), atol=1e-12)
