"""Small explicit-detector forward/adjoint example; no reconstruction or gap filling."""

from types import SimpleNamespace

import numpy as np

from tigre.utilities.Ax_surface import Atb_surface, Ax_surface
from tigre.utilities.detector_geometry import DetectorModule, DetectorSurfaceGeometry


def main():
    geo = SimpleNamespace(nVoxel=np.array([8, 16, 16]), dVoxel=np.ones(3), offOrigin=np.zeros(3))
    volume = np.ones(tuple(geo.nVoxel), dtype=np.float32)
    source = np.array([[40.0, 0.0, 0.0]])

    # Curved detector: explicit world coordinates for each sample.
    theta = np.linspace(-0.3, 0.3, 24)
    curved = np.zeros((8, 24, 3))
    curved[..., 0] = 40.0 - 80.0 * np.cos(theta)
    curved[..., 1] = 80.0 * np.sin(theta)
    curved[..., 2] = (np.arange(8) - 3.5)[:, None]

    # Three perpendicular planar panels arranged as a U, with missing layout cells.
    local = np.zeros((8, 12, 3))
    local[..., 1] = np.linspace(-6.0, 6.0, 12)
    local[..., 2] = (np.arange(8) - 3.5)[:, None]
    modules = []
    for index, (angle, centre) in enumerate(
        [
            (-np.pi / 2, [-34.0, -6.5, 0.0]),
            (0.0, [-40.0, 0.0, 0.0]),
            (np.pi / 2, [-34.0, 6.5, 0.0]),
        ]
    ):
        pose = np.eye(4)
        c, s = np.cos(angle), np.sin(angle)
        pose[:2, :2] = [[c, -s], [s, c]]
        pose[:3, 3] = centre
        modules.append(DetectorModule(local, pose, column_offset=13 * index))

    surfaces = {
        "curved": DetectorSurfaceGeometry(source, curved),
        "three-panel U": DetectorSurfaceGeometry.from_modules(source, modules),
    }
    for name, surface in surfaces.items():
        projections = Ax_surface(volume, geo, surface, gpuids=0)
        back = Atb_surface(projections, geo, surface, gpuids=0)
        assert projections.shape == surface.valid_mask.shape
        assert np.isfinite(back).all()
        assert np.all(projections[~surface.valid_mask] == 0)
        print(
            f"{name}: {projections.shape}, {surface.valid_mask.sum()} valid rays, "
            f"projection norm {np.linalg.norm(projections):.6g}"
        )


if __name__ == "__main__":
    main()
