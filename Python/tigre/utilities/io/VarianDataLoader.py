from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from tigre.utilities.geometry import Geometry
from tigre.utilities.io.varian.utils import PathLike
from tigre.utilities.io.varian.varian_io import (
    ProjData,
    ScanParams,
    ReconParams,
    read_varian_geometry,
    load_projections,
    load_blank_projections,
)
from tigre.utilities.io.varian.scatter import (
    DetScattParams,
    ScattParams,
    correct_detector_scatter,
    correct_scatter,
)
from scipy.ndimage import median_filter
from tqdm import tqdm


def VarianDataLoader(filepath: PathLike, **kwargs) -> tuple[NDArray, Geometry, NDArray]:
    """Loads raw projection data (xim format) for CBCT scans acquired with Varian OBI
    (Truebeam 2.0 or 2.7). Option to perform detector scatter correction (dps) and FASKS scatter
    correction (sc) based on algorithm described in Sun & Star-Lack 2010
    (doi: 10.1088/0031-9155/55/22/007). Ring artifact correction is also applied.
    NOTE: This function has only been tested on clinical data from Varian Truebeam (Ver 2.7).

    Args:
        filepath (PathLike): folder containing projection data and calibration files.
        kwargs:
            acdc (bool): acceleration-deceleration correction (default: True)
            dps (bool): detector point scatter correction (default: True)
            sc (bool): kernel-based scatter correction (default: True)

    Returns:
        tuple[NDArray, Geometry, NDArray]: log-normalized projections, geometry, projection angles (in radians)
    """
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
    """Applies median filter along column dimension to reduce ring artifacts."""
    print("Performing ring artifact correction:")
    log_projs = np.array([median_filter(p, size=kernel_size) for p in tqdm(log_projs)])
    return log_projs


def log_normalize(proj_data: ProjData, blank_proj_data: ProjData) -> NDArray:
    """Applies log normalization: p = -log(I/I0).

    Returns:
        NDArray: log normalized projections
    """
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
