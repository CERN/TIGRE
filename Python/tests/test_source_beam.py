"""Flux-contour semantics, elliptical collimation and matched native operators."""
from types import SimpleNamespace

import numpy as np
import pytest

from tigre.utilities.source_beam import SourceBeam
from tigre.utilities.detector_geometry import DetectorSurfaceGeometry
from tigre.utilities.geometry_io import save_geometry_yaml, load_geometry_file


def beam(aperture=None, fraction=.7):
    return SourceBeam([[1, 0, 0]], [0, 1, 0], [60, 30], fraction, aperture)


def test_reference_ellipse_is_not_a_cutoff():
    b = beam()
    wide, narrow = np.tan(np.deg2rad([30, 15]))
    directions = [[1, 0, 0], [1, wide, 0], [1, 0, narrow],
                  [1, wide / np.sqrt(2), narrow / np.sqrt(2)], [1, 2*wide, 0], [-1, 0, 0]]
    np.testing.assert_allclose(b.relative_flux(directions), [1, .7, .7, .7, .7**4, 0])
    pixels = np.asarray(directions)[None, None].astype(np.float32)
    surface = DetectorSurfaceGeometry([[0, 0, 0]], pixels, source_beam=b)
    assert surface.valid_mask.all()  # Reference ellipse and rear rays do not create an aperture.


def test_aperture_is_elliptical_and_preserves_gaps():
    b = SourceBeam([[1, 0, 0]], [0, 1, 0], [60, 30], .5, [90, 60])
    t = np.tan(np.deg2rad(30))
    pixels = np.array([[[1, .8, .8*t], [1, .5, .5*t], [1, 0, 0], [-1, 0, 0]]])
    s = DetectorSurfaceGeometry([[0, 0, 0]], pixels, [[True, True, False, True]], source_beam=b)
    np.testing.assert_array_equal(s.valid_mask[0, 0], [False, True, False, False])
    assert s.detector_valid_mask[0, 0, 0]  # Outside ellipse, inside rectangular bounds.
    assert not s.subset([0]).valid_mask[0, 0, 0]


@pytest.mark.parametrize('changes', [
    {'axis_directions': [[0, 0, 0]]}, {'axis_directions': [[np.nan, 0, 0]]},
    {'wide_axis_directions': [1, 1, 0]}, {'opening_angles_deg': [20, 40]},
    {'opening_angles_deg': [180, 20]}, {'opening_angles_deg': [60, 0]},
    {'reference_fraction': 0}, {'reference_fraction': 1}, {'reference_fraction': np.nan},
    {'aperture_angles_deg': [40, -10]}, {'opening_angles_deg': [60, 30, 10]}])
def test_invalid_beams(changes):
    values = dict(axis_directions=[[1, 0, 0]], wide_axis_directions=[0, 1, 0], opening_angles_deg=[60, 30])
    with pytest.raises(ValueError):
        SourceBeam(**(values | changes))


def surface():
    yy, zz = np.meshgrid(np.linspace(-80, 80, 25), np.linspace(-80, 80, 21))
    pixels = np.stack([np.full_like(yy, -100), yy, zz], axis=-1)
    b = SourceBeam([[-1, 0, 0]] * 3, [0, 1, 0], [60, 30], .7, [40, 20])
    mask = np.ones((21, 25), bool)
    mask[:, 12] = False
    return DetectorSurfaceGeometry([[100, 0, 0]] * 3, pixels, mask, source_beam=b)


@pytest.mark.parametrize('exact', [False, True])
@pytest.mark.parametrize('aperture', [None, [40, 20]])
def test_beam_archives_and_mat_bridge(tmp_path, exact, aperture):
    from tigre.utilities.geometry_matlab_bridge import to_mat, from_mat
    target = tmp_path / 'beam.yaml'
    original = surface()
    if aperture is None:
        original = DetectorSurfaceGeometry(original.source_positions_mm, original.pixel_centers_mm,
            original.detector_valid_mask, source_beam=SourceBeam([[-1, 0, 0]] * 3, [0, 1, 0], [60, 30], .7))
    if exact:
        save_geometry_yaml(target, surface=original)
    else:
        from tigre.utilities.detector_geometry_builder import build_modular_surface
        settings = dict(layout="Piecewise linear", n_views=3, source_radius_mm=300,
            source_z_mm=0, module_shape=[5, 7], pixel_spacing_mm=[5, 5],
            detector_radius_mm=200, column_gap=0, n_modules=1, angular_span_deg=20,
            transverse_centers_mm=[0], normal_offsets_mm=[0], panel_angles_deg=[0], z_offsets_mm=[0], source_beam=dict(
            opening_angles_deg=[60, 30], reference_fraction=.7, aperture_angles_deg=aperture))
        original = build_modular_surface(**settings)
        save_geometry_yaml(target, settings=settings)
    restored = load_geometry_file(target).surface
    np.testing.assert_array_equal(restored.valid_mask, original.valid_mask)
    np.testing.assert_array_equal(restored.detector_valid_mask, original.detector_valid_mask)
    np.testing.assert_allclose(restored.source_beam.opening_angles_deg, original.source_beam.opening_angles_deg)
    mat = tmp_path / 'beam.mat'
    to_mat(target, mat)
    from_mat(mat, tmp_path / 'roundtrip.yaml')
    restored = load_geometry_file(tmp_path / 'roundtrip.yaml').surface
    np.testing.assert_array_equal(restored.valid_mask, original.valid_mask)
    np.testing.assert_allclose(restored.source_beam.reference_fraction, .7)


def test_object_rotation_preserves_beam_and_aperture():
    from tigre.utilities.object_motion import fixed_object_geometry
    s = surface()
    g = SimpleNamespace(sVoxel=np.array([32.] * 3), offOrigin=np.zeros(3), object_angles_deg=[0, 30, 90])
    transformed, _ = fixed_object_geometry(s, g)
    np.testing.assert_array_equal(transformed.valid_mask, s.valid_mask)
    np.testing.assert_allclose(transformed.source_beam.axis_directions[-1], [0, 1, 0], atol=1e-7)


def test_native_collimation_matches_masked_operator_and_adjoint():
    from tigre.utilities.Ax_surface import Ax_surface, Atb_surface
    s = surface()
    plain = DetectorSurfaceGeometry(s.source_positions_mm, s.pixel_centers_mm, s.detector_valid_mask)
    g = SimpleNamespace(nVoxel=np.array([16] * 3), dVoxel=np.array([2.] * 3), offOrigin=np.zeros(3))
    rng = np.random.default_rng(13)
    x = rng.random((16, 16, 16), dtype=np.float32)
    y = rng.random(s.valid_mask.shape, dtype=np.float32)
    ax = Ax_surface(x, g, s)
    np.testing.assert_allclose(ax, Ax_surface(x, g, plain) * s.valid_mask, atol=1e-5)
    assert np.all(ax[~s.valid_mask] == 0)
    aty = Atb_surface(y, g, s)
    np.testing.assert_allclose(aty, Atb_surface(y * s.valid_mask, g, plain), rtol=2e-5, atol=2e-5)
    np.testing.assert_allclose(np.vdot(ax, y), np.vdot(x, aty), rtol=2e-5)
