import numpy as np
import glob
import os
import xml.etree.ElementTree as ET
from scipy.ndimage import gaussian_filter, median_filter
from scipy.fft import fft2, ifft2
from scipy.signal import decimate, convolve2d
from scipy.interpolate import interpn
from Python.tigre.utilities.io.varian.utils import get_xmlns, interpolate_blank_scan
from tqdm import tqdm


def read_scatter_calib(filepath):

    file = glob.glob(os.path.join(filepath, "Calibrations", "SC-*", "Factory", "Calibration.xml"))
    if len(file) == 0:
        raise RuntimeError("Scatter calibration file not found.")
    elif len(file) > 1:
        raise RuntimeError("Multiple calibration files found")

    tree = ET.parse(file[0])
    sc_calib = tree.getroot()

    return sc_calib


def calculate_grid_response(u, v, grid_efficiency):
    # TODO: read params from Calibration.xml
    K = -0.15
    B = 1
    t_ratio = K * abs(v.T) + B
    kernel = np.tile(t_ratio[:, np.newaxis], [1, len(u)])
    kernel[kernel < grid_efficiency] = grid_efficiency
    return kernel


def _log_normalize(blank, proj):
    eps = np.finfo(proj.dtype).eps
    ratio = blank / (proj + eps)
    ratio[ratio < 1] = 1  # TODO: replace with interpolation of these values?
    log_norm = np.log(ratio)
    return log_norm


def mm2cm(x):
    return 0.1 * x


def estimate_water_equiv_thickness(blank, proj, mu_water, smooth=False, sigma=None):
    log_norm = _log_normalize(blank, proj)
    thickness = log_norm / mu_water
    if smooth:
        if sigma is None:
            raise ValueError("sigma cannot be 'None' if smooth=True.")
        thickness = gaussian_filter(thickness, sigma)
    return thickness


def get_thickness_masks(thickness_map, lower_bounds):
    mask = np.zeros_like(thickness_map)
    masks = np.tile(mask[np.newaxis, :, :], [len(lower_bounds), 1, 1])
    for i in range(len(lower_bounds) - 1):
        masks[i] = (thickness_map > lower_bounds[i]) * (thickness_map < lower_bounds[i + 1])
    masks[-1] = thickness_map > lower_bounds[-1]
    return masks


def _get_form_func_params(sc_calib):
    ns = get_xmlns(sc_calib)
    obj_scatter_models = sc_calib.find("CalibrationResults/ObjectScatterModels", ns)
    obj_scatter_fits = obj_scatter_models.findall("ObjectScatterModel/ObjectScatterFit", ns)
    sigma1 = [float(elem.find("sigma1", ns).text) for elem in obj_scatter_fits]  # (cm^-1)
    sigma2 = [float(elem.find("sigma2", ns).text) for elem in obj_scatter_fits]  # (cm^-1)
    B = [float(elem.find("B", ns).text) for elem in obj_scatter_fits]
    return sigma1, sigma2, B


def calculate_form_functions(thickness_map, sc_calib, grid_coords):
    sigma1, sigma2, B = _get_form_func_params(sc_calib)
    num_groups = len(sigma1)
    U, V = grid_coords
    grid = U**2 + V**2
    gform = np.zeros_like(thickness_map)
    gforms = np.tile(gform[np.newaxis, :, :], [num_groups, 1, 1])

    for i in range(num_groups):
        gforms[i] = np.exp(-0.5 * grid / (sigma1[i] ** 2)) + B[i] * np.exp(
            -0.5 * grid / (sigma2[i] ** 2)
        )
    return gforms


def _rescale(x, min_val=0.0, max_val=1.0):
    x_rescale = (max_val - min_val) * (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)) + min_val
    return x_rescale


def calculate_edge_response(thickness_map, threshold=50, filter_size=25, num_iter=5):
    # TODO: fetch params from calibration.xml file
    edge_weight = thickness_map > threshold
    edge_weight = np.array(edge_weight, dtype="float")
    tmp_mask = edge_weight
    avg_kernel = np.ones([filter_size, filter_size]) / filter_size

    for i in range(num_iter):
        edge_weight = convolve2d(edge_weight, avg_kernel, mode="same")

    tmp = tmp_mask * edge_weight
    tmp[tmp == 0] = np.nan
    tmp = _rescale(tmp, min_val=0.6)
    tmp[np.isnan(tmp)] = 0
    edge_weight = tmp
    return edge_weight


def _get_ampl_params(sc_calib):
    ns = get_xmlns(sc_calib)
    obj_scatter_models = sc_calib.find("CalibrationResults/ObjectScatterModels", ns)
    obj_scatter_fits = obj_scatter_models.findall("ObjectScatterModel/ObjectScatterFit", ns)
    A = [float(elem.find("A", ns).text) for elem in obj_scatter_fits]
    alpha = [float(elem.find("alpha", ns).text) for elem in obj_scatter_fits]
    beta = [float(elem.find("beta", ns).text) for elem in obj_scatter_fits]
    return A, alpha, beta


def calculate_amplitudes(blank, proj, edge_weight, sc_calib):
    A, alpha, beta = _get_ampl_params(sc_calib)
    num_groups = len(A)
    log_norm = _log_normalize(blank, proj)
    norm = proj / blank
    # TODO: improve thresholding with interpolation?
    norm[norm > 1] = 1
    norm[norm < 0] = 0
    amplitudes = []
    for i in range(num_groups):
        amplitude = mm2cm(A[i]) * edge_weight * (norm ** alpha[i]) * (log_norm ** beta[i])
        amplitudes.append(amplitude)
    return np.array(amplitudes)


def update_scatter_estimate(
    proj, thickness_map, thickness_masks, amplitudes, gforms, grid_kernel, gamma
):
    num_groups = len(thickness_masks)
    rep_proj = np.tile(proj[np.newaxis, :, :], [num_groups, 1, 1])
    term1 = rep_proj * thickness_masks * amplitudes
    rep_grid_kernel = np.tile(grid_kernel[np.newaxis, :, :], [num_groups, 1, 1])
    term2 = fft2(gforms * rep_grid_kernel)
    rep_thickness = np.tile(thickness_map[np.newaxis, :, :], [num_groups, 1, 1])
    term3 = fft2(rep_thickness * term1)
    tmp1 = np.sum(fft2(term1) * term2, axis=0)
    tmp2 = np.sum(term2 * term3, axis=0)
    itmp1 = ifft2(tmp1)
    itmp2 = ifft2(tmp2)
    scatter = (1 - gamma * thickness_map) * np.real(itmp1) + gamma * np.real(itmp2)
    return scatter


def update_primary_estimate(primary, scatter_old, scatter, lam=0.6):
    primary += lam * (scatter_old - scatter)
    eps = np.finfo(primary.dtype).eps
    primary[primary < eps] = eps  # TODO improve thresholding
    return primary


def calculate_primary(proj, scatter, max_scatter_frac=0.95):
    scatter_frac = scatter / (proj + np.finfo(proj.dtype).eps)
    # TODO: add nan check?
    scatter_frac = median_filter(scatter_frac, size=3)
    scatter_frac = np.minimum(scatter_frac, max_scatter_frac)
    return proj * (1 - scatter_frac)


def correct_scatter(
    projs,
    angles,
    airnorms,
    blank_projs,
    blank_angles,
    blank_airnorms,
    geometry,
    sc_calib,
    downsample=12,
    num_iter=8,
    lam=0.005,
):
    # Detector coords,centered (cm)
    u, v = _get_detector_coords(geometry)
    U, V = np.meshgrid(u, v)

    # Downsampled coords
    du, dv = _get_detector_coords(geometry, downsample=downsample)
    DU, DV = np.meshgrid(du, dv)

    ns = get_xmlns(sc_calib)
    mu_water = float(sc_calib.find("CalibrationResults/Globals/muH2O", ns).text)  # mm^-1
    # for thickness smoothing
    sigma_u = float(sc_calib.find("CalibrationResults/Globals/AsymPertSigmaMMu", ns).text)
    sigma_v = float(sc_calib.find("CalibrationResults/Globals/AsymPertSigmaMMv", ns).text)
    step_du = np.mean(np.diff(du))
    step_dv = np.mean(np.diff(dv))
    sigma = (mm2cm(sigma_v) / step_dv, mm2cm(sigma_u) / step_du)  # dimensionless

    obj_scatter_models = sc_calib.find("CalibrationResults/ObjectScatterModels", ns)
    gamma = float(obj_scatter_models.find("ObjectScatterModel/ObjectScatterFit/gamma", ns).text)
    thickness_bounds_mm = [
        float(thickness.text)
        for thickness in obj_scatter_models.findall("ObjectScatterModel/Thickness", ns)
    ]
    grid_efficiency = float(
        obj_scatter_models.find("ObjectScatterModel/GridEfficiency/LamellaTransmission", ns).text
    )
    grid_kernel = calculate_grid_response(du, dv, grid_efficiency)

    primaries = np.zeros_like(projs)
    print("Performing ASKS scatter correction: ")
    for i, proj in tqdm(enumerate(projs)):
        blank_interp, airnorm_interp = interpolate_blank_scan(
            angles[i], blank_projs, blank_angles, blank_airnorms
        )

        cf_air = airnorms[i] / airnorm_interp
        blank = interpn((v, u), blank_interp * cf_air, (DV, DU))
        primary = interpn((v, u), proj, (DV, DU))
        scatter = np.zeros_like(primary)  # initialize scatter

        for n in range(num_iter):
            scatter_old = scatter
            thickness_map = estimate_water_equiv_thickness(
                blank, primary, mu_water, smooth=True, sigma=sigma
            )
            edge_weight = calculate_edge_response(thickness_map)
            thickness_masks = get_thickness_masks(thickness_map, thickness_bounds_mm)
            gforms = calculate_form_functions(thickness_map, sc_calib, (DU, DV))
            amplitudes = calculate_amplitudes(blank, primary, edge_weight, sc_calib)

            scatter = update_scatter_estimate(
                primary,
                mm2cm(thickness_map),
                thickness_masks,
                amplitudes,
                gforms,
                grid_kernel,
                gamma,
            )
            primary = update_primary_estimate(primary, scatter_old, scatter, lam=lam)

        # Upsample
        scatter_est = interpn(
            (dv, du), scatter, (V, U), method="cubic", bounds_error=False, fill_value=None
        )
        scatter_est[scatter_est < 0] = 0
        # TODO: check for nans and interpolate
        primaries[i] = calculate_primary(proj, scatter_est)
    return primaries


def _get_dps_params(sc_calib):
    ns = get_xmlns(sc_calib)
    det_scatter_model = sc_calib.find("CalibrationResults/Globals/DetectorScatterModel", ns)
    return [float(elem.text) for elem in det_scatter_model]


def calculate_dps_kernel(sc_calib, U, V):
    grid = np.sqrt(U**2 + V**2)
    a = _get_dps_params(sc_calib)
    det_kernel = a[0] * np.exp(-a[1] * grid) + a[2] * (np.exp(-a[3] * (grid - a[4]) ** 3))
    det_kernel = a[5] * det_kernel / np.sum(det_kernel)
    return det_kernel


def _get_detector_coords(geometry, downsample=0):
    # Detector coords,centered (cm)
    u = (
        0.5
        * (geometry.sDetector[1] - geometry.dDetector[1])
        * np.linspace(-1, 1, int(geometry.nDetector[1]))
    )
    v = (
        0.5
        * (geometry.sDetector[0] - geometry.dDetector[0])
        * np.linspace(-1, 1, int(geometry.nDetector[0]))
    )

    if downsample:
        u = decimate(u, downsample)
        v = decimate(v, downsample)

    return mm2cm(u), mm2cm(v)


def correct_detector_scatter(projs, geometry, sc_calib, downsample=8):

    u, v = _get_detector_coords(geometry)
    U, V = np.meshgrid(u, v)
    du, dv = _get_detector_coords(geometry, downsample=downsample)
    DU, DV = np.meshgrid(du, dv)

    det_kernel = calculate_dps_kernel(sc_calib, DU, DV)

    corrected_projs = np.zeros_like(projs)
    print("Performing detector point scatter correction: ")
    for i, proj in tqdm(enumerate(projs)):
        proj_down = interpn((v, u), projs[0], (DV, DU))
        scatter_down = convolve2d(proj_down, det_kernel, mode="same")
        scatter = interpn(
            (dv, du), scatter_down, (V, U), method="cubic", bounds_error=False, fill_value=None
        )
        corrected_projs[i] = proj - scatter

    return corrected_projs
