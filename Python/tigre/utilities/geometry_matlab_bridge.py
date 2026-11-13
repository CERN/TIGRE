"""MAT/YAML/NPZ interchange for MATLAB geometry workflows; run with the recon interpreter."""
import argparse
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat, savemat

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tigre.utilities.geometry_io import load_geometry_file, save_geometry_yaml
from tigre.utilities.detector_geometry import DetectorSurfaceGeometry


def to_mat(source, target):
    document = load_geometry_file(source)
    surface = document.surface
    arrays = dict(source_positions_mm=surface.source_positions_mm.T,
                  pixel_centers_mm=surface.pixel_centers_mm.transpose(1, 2, 3, 0),
                  valid_mask=surface.valid_mask.transpose(1, 2, 0))
    if surface.source_beam is not None:
        arrays['detector_valid_mask'] = surface.detector_valid_mask.transpose(1, 2, 0)
        arrays['source_beam'] = {name: getattr(surface.source_beam, name) for name in
            ('axis_directions', 'wide_axis_directions', 'opening_angles_deg', 'reference_fraction', 'aperture_angles_deg')
            if getattr(surface.source_beam, name) is not None}
    if surface._module_slices:
        arrays['module_slices'] = np.array([[r.start, r.stop, c.start, c.stop]
                                          for r, c in surface._module_slices])
    result = {'surface': arrays, 'geo': np.empty((0, 0))}
    if document.settings is not None:
        def mat_value(value):
            if value is None:
                return np.empty((0, 0))
            if isinstance(value, dict):
                return {key: mat_value(item) for key, item in value.items()}
            return value
        result['settings'] = mat_value(document.settings)
    if document.geometry is not None:
        result['geo'] = dict(sVoxel=np.asarray(document.geometry.sVoxel)[::-1],
                             offOrigin=np.atleast_2d(document.geometry.offOrigin)[:, ::-1].T)
        if hasattr(document.geometry, 'nVoxel'):
            result['geo']['nVoxel'] = document.geometry.nVoxel[::-1]
            result['geo']['dVoxel'] = document.geometry.dVoxel[::-1]
        if hasattr(document.geometry, 'object_angles_deg'):
            from tigre.utilities.object_motion import object_angles
            result['geo']['object_angles_deg'] = object_angles(document.geometry, len(surface.source_positions_mm)).copy()
        for name in ('acquisition', 'motion_mode', 'view_times_s', 'simulation_mode'):
            if hasattr(document.geometry, name):
                result['geo'][name] = getattr(document.geometry, name)
    savemat(target, result)


def from_mat(source, target):
    document = loadmat(source, simplify_cells=True)
    geo = document.get('geo')
    volume = None
    if isinstance(geo, dict):
        from types import SimpleNamespace
        offsets = np.asarray(geo['offOrigin'])
        offsets = offsets[::-1] if offsets.ndim == 1 else offsets.T[:, ::-1]
        volume = SimpleNamespace(sVoxel=np.asarray(geo['sVoxel'])[::-1],
                                 offOrigin=offsets)
        if 'nVoxel' in geo:
            volume.nVoxel = np.asarray(geo['nVoxel'])[::-1]
        if 'object_angles_deg' in geo:
            volume.object_angles_deg = np.asarray(geo['object_angles_deg'], dtype=float).reshape(-1, 3)
        for name in ('acquisition', 'motion_mode', 'view_times_s', 'simulation_mode'):
            if name in geo:
                setattr(volume, name, geo[name])
    if isinstance(document.get('settings'), dict):
        settings = document['settings']
        for key in ('transverse_centers_mm', 'normal_offsets_mm', 'z_offsets_mm',
                    'panel_angles_deg', 'panel_x_angles_deg', 'panel_y_angles_deg'):
            settings[key] = np.atleast_1d(settings[key])
        if np.asarray(settings.get('hardware_angles_deg', [])).size == 0:
            settings['hardware_angles_deg'] = None
        else:
            settings['hardware_angles_deg'] = np.asarray(settings['hardware_angles_deg'], dtype=float).reshape(-1, 3)
        if not isinstance(settings.get('acquisition'), dict):
            settings['acquisition'] = None
        if not isinstance(settings.get('source_beam'), dict):
            settings['source_beam'] = None
        elif np.asarray(settings['source_beam'].get('aperture_angles_deg', [])).size == 0:
            settings['source_beam']['aperture_angles_deg'] = None
        for key in ('n_views', 'n_modules', 'column_gap'):
            value = float(settings[key])
            if not value.is_integer():
                raise ValueError(f'{key} must be an integer')
            settings[key] = int(value)
        shape = np.asarray(settings['module_shape'])
        if not np.equal(shape, np.floor(shape)).all():
            raise ValueError('module_shape must contain integers')
        settings['module_shape'] = shape.astype(int)
        if 'top_arm_rows' in settings:
            value = float(settings['top_arm_rows'])
            if not value.is_integer():
                raise ValueError('top_arm_rows must be an integer')
            settings['top_arm_rows'] = int(value)
        if np.asarray(settings.get('source_start_mm', [])).size == 0:
            settings['source_start_mm'] = None
        for key in ('source_waypoints_mm', 'source_segment_views'):
            if np.asarray(settings.get(key, [])).size == 0:
                settings[key] = None
        if settings['source_segment_views'] is not None:
            counts = np.atleast_1d(np.asarray(settings['source_segment_views'], dtype=float))
            if (not np.isfinite(counts).all() or np.any(counts != np.floor(counts))
                    or np.any(counts < 2) or np.any(counts > 6000)):
                raise ValueError('source_segment_views must contain integers from 2 to 6000')
            settings['source_segment_views'] = counts.astype(np.int64)
        save_geometry_yaml(target, settings=settings, geometry=volume)
    else:
        # Read unsqueezed arrays: singleton rows/views carry physical meaning.
        raw = loadmat(source, struct_as_record=False, squeeze_me=False)['surface'][0, 0]
        pixels = raw.pixel_centers_mm
        if pixels.ndim == 3:
            pixels = pixels[..., None]
        mask = raw.valid_mask
        if mask.ndim == 2:
            mask = mask[..., None]
        modules = getattr(raw, 'module_slices', None)
        if modules is not None:
            modules = [(slice(int(a), int(b)), slice(int(c), int(d))) for a, b, c, d in modules]
        beam = None
        if hasattr(raw, 'source_beam'):
            from tigre.utilities.source_beam import SourceBeam
            saved = raw.source_beam[0, 0]
            values = {name: getattr(saved, name) for name in
                ('axis_directions', 'wide_axis_directions', 'opening_angles_deg', 'reference_fraction', 'aperture_angles_deg')
                if hasattr(saved, name)}
            values['reference_fraction'] = values['reference_fraction'].ravel()
            beam = SourceBeam(**values)
            mask = raw.detector_valid_mask
            if mask.ndim == 2:
                mask = mask[..., None]
        surface = DetectorSurfaceGeometry(raw.source_positions_mm.T,
            pixels.transpose(3, 0, 1, 2), mask.transpose(2, 0, 1), module_slices=modules, source_beam=beam)
        save_geometry_yaml(target, surface=surface, geometry=volume)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['to-mat', 'from-mat'])
    parser.add_argument('source')
    parser.add_argument('target')
    args = parser.parse_args()
    (to_mat if args.operation == 'to-mat' else from_mat)(args.source, args.target)
