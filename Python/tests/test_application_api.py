"""Independent public-client checks; no UI, implicit device or solver state."""
import importlib
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from tigre import design_api as api
from test_projection_operator import MatrixOperator


def test_preserves_design_api_10_and_reports_features_without_probing_cuda():
    api.require_api(1, 0)
    api.require_api(1, 1)
    assert api.API_VERSION == (1, 1)
    assert len(api.__all__) == len(set(api.__all__))
    with patch('tigre.utilities.gpu.getGpuNames', side_effect=AssertionError('implicit device query')):
        capabilities = api.backend_capabilities()
    assert len(capabilities['algorithms']) == 26
    assert not capabilities['surface_fdk']
    assert not capabilities['iteration_progress']
    assert not capabilities['in_flight_cancellation']
    assert not capabilities['cuda_runtime_checked']


def test_owned_result_diagnostics_and_masked_data_match_dense_least_squares():
    op = MatrixOperator()
    data = op.forward(np.array([1, 2, 3, 0], np.float32).reshape(op.volume_shape))
    data[~op.valid_mask] = np.nan
    initial = np.full(op.volume_shape, .25, np.float32)
    options = {'init': initial, 'return_residuals': True}
    events = []
    result = api.reconstruct(data, operator=op, niter=3, options=options, progress=events.append)
    expected = np.linalg.lstsq(op.matrix, np.nan_to_num(data).ravel(), rcond=None)[0]
    expected[-1] = .25
    np.testing.assert_allclose(result.volume.ravel(), expected, atol=2e-5)
    assert [event.stage for event in events] == ['validating', 'running', 'completed']
    assert all(event.algorithm == 'cgls' and event.requested_iterations == 3 for event in events)
    assert result.diagnostics[-1] < result.diagnostics[0]
    assert result.elapsed_seconds >= 0
    assert not np.shares_memory(result.volume, initial)
    np.testing.assert_array_equal(initial, .25)
    assert np.isnan(data[~op.valid_mask]).all()
    assert options['init'] is initial


def test_cancellation_before_native_and_after_native_discards_result():
    op = MatrixOperator()
    data = np.zeros(op.projection_shape, np.float32)
    with patch('tigre.algorithms.cgls_operator', side_effect=AssertionError('must not execute')):
        with pytest.raises(api.ReconstructionCancelled):
            api.reconstruct(data, operator=op, cancelled=lambda: True)
    stop = False
    def backend(*args, **kwargs):
        nonlocal stop
        stop = True
        return np.ones(op.volume_shape, np.float32)
    events = []
    with patch('tigre.algorithms.cgls_operator', side_effect=backend):
        with pytest.raises(api.ReconstructionCancelled, match='discarded'):
            api.reconstruct(data, operator=op, cancelled=lambda: stop, progress=events.append)
    assert [event.stage for event in events] == ['validating', 'running']


def test_progress_callback_can_cancel_before_execution_and_errors_propagate():
    op = MatrixOperator()
    stop = False
    def progress(event):
        nonlocal stop
        if event.stage == 'running': stop = True
    with patch('tigre.algorithms.cgls_operator', side_effect=AssertionError('must not execute')):
        with pytest.raises(api.ReconstructionCancelled):
            api.reconstruct(np.zeros(op.projection_shape), operator=op, progress=progress, cancelled=lambda: stop)
    class CallbackError(Exception): pass
    def broken(event): raise CallbackError('caller callback')
    with pytest.raises(CallbackError):
        api.reconstruct(np.zeros(op.projection_shape), operator=op, progress=broken)


def test_backend_failure_retains_cause_and_does_not_report_completion():
    op = MatrixOperator()
    cause = RuntimeError('device allocation failed')
    events = []
    with patch('tigre.algorithms.cgls_operator', side_effect=cause):
        with pytest.raises(api.BackendExecutionError, match='cgls') as caught:
            api.reconstruct(np.zeros(op.projection_shape), operator=op, progress=events.append)
    assert caught.value.__cause__ is cause
    assert [event.stage for event in events] == ['validating', 'running']


@pytest.mark.parametrize('options', [dict(algorithm='fdk'), dict(niter=True), dict(niter=-1),
    dict(progress=1), dict(cancelled=1), dict(options=[]), dict(options={'gpuids': 0}),
    dict(algorithm='pcsd'), dict(algorithm='pcsd', options={'maxl2err': np.nan})])
def test_invalid_requests_do_not_execute_backend(options):
    op = MatrixOperator()
    with patch('tigre.algorithms.cgls_operator', side_effect=AssertionError('must not execute')):
        with pytest.raises((ValueError, TypeError)):
            api.reconstruct(np.zeros(op.projection_shape), operator=op, **options)


def test_explicit_device_discovery_and_binding_ownership():
    with patch('tigre.utilities.gpu.getGpuNames', return_value=['one', 'two']):
        assert api.available_devices() == (api.DeviceInfo(0, 'one'), api.DeviceInfo(1, 'two'))
    with patch('tigre.utilities.gpu.getGpuNames', side_effect=RuntimeError('driver missing')):
        with pytest.raises(api.BackendExecutionError, match='discovery') as error:
            api.available_devices()
        assert str(error.value.__cause__) == 'driver missing'
    from test_surface_projectors import fixture
    geo, surface = fixture()
    op = api.create_operator(geo, surface=surface, gpuids=[1, 0], n_slabs=3)
    assert op.gpuids == [1, 0] and op.n_slabs == 3
    geo.offOrigin[:] = 999
    assert not np.all(op.geo.offOrigin == 999)
    with pytest.raises(ValueError, match='exactly one'):
        api.create_operator(geo, surface=surface, angles=[0])
    with pytest.raises(ValueError, match='omit geometry'):
        api.forward_project(np.ones(op.volume_shape), geo, operator=op)


def test_forward_service_preserves_inputs_and_pocs_diagnostics():
    op = MatrixOperator()
    image = np.ones(op.volume_shape, np.float32)
    data = api.forward_project(image, operator=op)
    np.testing.assert_array_equal(image, 1)
    result = api.reconstruct(data, operator=op, algorithm='pcsd', niter=2,
                             options={'maxl2err': 0., 'tviter': 0, 'return_history': True})
    assert isinstance(result.diagnostics, dict) and len(result.diagnostics['residual']) == 2
    alias = api.reconstruct(data, operator=op, algorithm='ossart', niter=1)
    assert alias.algorithm == 'os_sart'


def test_frozen_design_api_10_names_are_retained():
    expected={'load_head_phantom', 'volume_from_mapping', 'rotation_zyz', 'volume_box', 'MAX_MODULE_DIMENSION', 'build_modular_surface', 'remap_projections_for_fdk', 'MAX_TOTAL_PIXELS', 'require_api', 'surface_from_geometry', 'view_angles', 'MAX_MODULES', 'panel_tilts', 'panel_rotation', 'API_VERSION', 'SurfaceRenderer', 'save_geometry_yaml', 'validate_surface_size', 'check_tigre_compatibility', 'checked_module_shape', 'cargo_module_layout', 'with_source_beam', 'panel_count_for_layout', 'source_segments', 'MAX_RECIPE_VIEWS', 'fixed_object_geometry', 'load_geometry_archive', 'validate_transmission_geometry', 'draw_world_axes', 'surface_archive_arrays', 'acquisition_timing', 'piecewise_source_positions', 'DetectorSurfaceGeometry', 'load_geometry_file', 'linear_source_surface', 'object_angles'}
    assert expected.issubset(set(api.__all__))
    assert all(hasattr(api,name) for name in expected)
