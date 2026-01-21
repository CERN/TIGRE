from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from Python.tigre.utilities.geometry import Geometry
from Python.tigre.utilities.io.varian.utils import PathLike
from Python.tigre.utilities.io.varian.varian_io import (
    ProjData,
    ScanParams,
    ReconParams,
    read_varian_geometry,
    load_projections,
    load_blank_projections,
)
from Python.tigre.utilities.io.varian.scatter import (
    DetScattParams,
    ScattParams,
    correct_detector_scatter,
    correct_scatter,
)
from scipy.ndimage import median_filter
from tqdm import tqdm


def VarianDataLoader(filepath: PathLike, **kwargs) -> tuple[NDArray, Geometry, NDArray]:

    acdc, dps, sc = parse_inputs(**kwargs)

    scan_params = ScanParams(filepath)
    recon_params = ReconParams(filepath)
    if recon_params.root is None:
        recon_params = None
    geometry = read_varian_geometry(scan_params, recon_params)

    blank_proj_data = load_blank_projections(filepath, scan_params)
    if acdc:
        angular_threshold = scan_params.calculate_angular_threshold()
        proj_data = load_projections(filepath, angular_threshold)
    else:
        proj_data = load_projections(filepath)

    if dps:
        dps_calib = DetScattParams(filepath)
        blank_proj_data.projs = correct_detector_scatter(blank_proj_data.projs, geometry, dps_calib)
        proj_data.projs = correct_detector_scatter(proj_data.projs, geometry, dps_calib)

    if sc:
        sc_calib = ScattParams(filepath)
        proj_data.projs = correct_scatter(proj_data, blank_proj_data, geometry, sc_calib)

    log_projs = log_normalize(proj_data, blank_proj_data)
    log_projs = correct_ring_artifacts(log_projs)

    return log_projs, geometry, np.deg2rad(proj_data.angles)


def correct_ring_artifacts(log_projs: NDArray, kernel_size: tuple[int, int] = (1, 9)) -> NDArray:
    print("Performing ring artifact correction:")
    log_projs = np.array([median_filter(p, size=kernel_size) for p in tqdm(log_projs)])
    return log_projs


def log_normalize(proj_data: ProjData, blank_proj_data: ProjData) -> NDArray:
    proj_data.projs = proj_data.projs.astype("float32")
    blank_proj_data.projs = blank_proj_data.projs.astype("float32")
    log_projs = np.zeros_like(proj_data.projs)
    eps = np.finfo(proj_data.projs.dtype).eps

    print("Performing log normalization:")
    for i in tqdm(range(len(proj_data.angles))):
        blank_interp = blank_proj_data.interp_proj(proj_data.angles[i])
        ratio = (blank_interp / (proj_data.projs[i] + eps)) + eps
        ratio[ratio < 1] = 1
        log_projs[i] = np.log(ratio)
    return log_projs


def parse_inputs(**kwargs) -> tuple[bool, bool, bool]:
    """
    Returns tags.
    ACDC: acceleration-deceleration correction (default: True)
    DPS: detector point scatter correction (default: True)
    SC: kernel-based scatter correction (default: True)
    """

    acdc = kwargs["acdc"] if "acdc" in kwargs else True
    dps = kwargs["dps"] if "dps" in kwargs else True
    sc = kwargs["sc"] if "sc" in kwargs else True
    if sc and not dps:
        RuntimeWarning("dps should be enabled when sc=true.")
    return acdc, dps, sc
