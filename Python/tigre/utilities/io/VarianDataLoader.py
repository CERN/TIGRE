import numpy as np
import os
import glob
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET
from Python.tigre.utilities.io.xim_io import XIM
from scipy.ndimage import gaussian_filter
from scipy.fft import fft2, ifft2
from scipy.signal import decimate, convolve2d
from scipy.interpolate import interpn
from scipy.ndimage import median_filter
from tqdm import tqdm


def VarianDataLoader(filepath, **kwargs):

    acdc, dps, sc = parse_inputs(**kwargs)

    geometry, scan_params = read_varian_geometry(filepath)
    if dps or sc:
        sc_calib = read_scatter_calib(filepath)

    angular_threshold = 0
    if acdc:
        ns = get_xmlns(scan_params)
        rot_velocity = float(scan_params.find("Acquisitions/Velocity", ns).text)
        frame_rate = float(scan_params.find("Acquisitions/FrameRate", ns).text)
        angular_threshold = calculate_angular_threshold(rot_velocity, frame_rate)

    projs, angles, airnorms = load_projections(filepath, angular_threshold)

    blank_projs, blank_angles, blank_airnorms = load_blank_projections(filepath, scan_params)

    if dps:
        blank_projs = correct_detector_scatter(blank_projs, geometry, sc_calib)
        projs = correct_detector_scatter(projs, geometry, sc_calib)

    if sc:
        projs = correct_scatter(
            projs, angles, airnorms, blank_projs, blank_angles, blank_airnorms, geometry, sc_calib
        )

    log_projs = log_normalize(projs, angles, airnorms, blank_projs, blank_angles, blank_airnorms)
    log_projs = enforce_positive(log_projs)
    log_projs = correct_ring_artifacts(log_projs)

    return log_projs, geometry, np.deg2rad(angles)


def get_xmlns(xml_root):
    return {"": xml_root.tag.split("}")[0].strip("{")}


def read_varian_geometry(filepath):

    # read Scan.xml
    file = os.path.join(filepath, "Scan.xml")
    tree = ET.parse(file)
    scan_params = tree.getroot()
    ns = get_xmlns(scan_params)
    acquisition_params = scan_params.find("Acquisitions", ns)
    if acquisition_params is None:
        raise RuntimeError("Cone beam scan parameters not found")

    geometry = Geometry()
    geometry.mode = "cone"
    geometry.DSD = float(acquisition_params.find("SID", ns).text)
    geometry.DSO = float(acquisition_params.find("SAD", ns).text)
    geometry.nDetector = np.array(
        [
            int(acquisition_params.find("ImagerSizeY", ns).text),
            int(acquisition_params.find("ImagerSizeX", ns).text),
        ]
    )
    geometry.dDetector = np.array(
        [
            float(acquisition_params.find("ImagerResY", ns).text),
            float(acquisition_params.find("ImagerResX", ns).text),
        ]
    )
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    offset = float(acquisition_params.find("ImagerLat", ns).text)
    geometry.offDetector = np.array([0, -offset])
    geometry.offOrigin = np.array([0, 0, 0])

    # Auxiliary variable to define accuracy of 'interpolated' projection
    # Recommended <=0.5 (vx/sample)
    geometry.accuracy = 0.5

    # read Reconstruction.xml
    file = glob.glob(os.path.join(filepath, "**", "Reconstruction.xml"), recursive=True)
    if len(file) > 1:
        raise RuntimeError("Multiple Reconstruction.xml files found.")
    tree = ET.parse(file[0])
    recon_params = tree.getroot()
    ns = get_xmlns(recon_params)

    if recon_params is None:
        print("Estimating acceptable image size...")
        geometry.dVoxel = np.array(
            [geometry.dDetector[1], geometry.dDetector[0], geometry.dDetector[0]]
        ) * (geometry.DSO / geometry.DSD)
        geometry.nVoxel = np.ceil(
            np.array(
                [
                    geometry.nDetector[1],
                    geometry.nDetector[0] + abs(geometry.offDetector[0]) / geometry.dDetector[0],
                    geometry.nDetector[0] + abs(geometry.offDetector[0]) / geometry.dDetector[0],
                ]
            )
        ).astype("int")
        geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    else:
        # Retrieve existing reconstruction parameters
        geometry.sVoxel = np.array(
            [
                float(recon_params.find("VOISizeZ", ns).text),
                float(recon_params.find("VOISizeY", ns).text),
                float(recon_params.find("VOISizeX", ns).text),
            ]
        )
        slice_num = np.ceil(
            geometry.sVoxel[0] / float(recon_params.find("SliceThickness", ns).text)
        )
        geometry.nVoxel = np.array(
            [
                slice_num,
                int(recon_params.find("MatrixSize", ns).text),
                int(recon_params.find("MatrixSize", ns).text),
            ]
        )
        geometry.dVoxel = geometry.sVoxel / geometry.nVoxel

    return geometry, scan_params


def calculate_angular_threshold(rot_velocity, frame_rate):
    return 0.9 * (rot_velocity / frame_rate)


def correct_ring_artifacts(log_projs, kernel_size=(1, 9)):
    log_projs = np.array([median_filter(p, size=kernel_size) for p in tqdm(log_projs)])
    return log_projs


def enforce_positive(x):
    return np.clip(x, a_min=0, a_max=None)


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
    K = -0.15  # TODO: determine where this number is in the cited article
    B = 1
    t_ratio = K * abs(v.T) + B
    kernel = np.tile(t_ratio[:, np.newaxis], [1, len(u)])
    kernel[kernel < grid_efficiency] = grid_efficiency
    return kernel


def log_norm_func(blank, proj):
    # TODO CHANGE FUNCTION NAME
    eps = np.finfo(proj.dtype).eps
    ratio = blank / (proj + eps)
    ratio[ratio < 1] = 1  # TODO: replace with interpolation of these values?
    log_norm = np.log(ratio)
    return log_norm


def estimate_water_equiv_thickness(blank, proj, mu_water, smooth=False, sigma=None):
    log_norm = log_norm_func(blank, proj)
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


def get_form_func_params(sc_calib):
    ns = get_xmlns(sc_calib)
    obj_scatter_models = sc_calib.find("CalibrationResults/ObjectScatterModels", ns)
    obj_scatter_fits = obj_scatter_models.findall("ObjectScatterModel/ObjectScatterFit", ns)
    # for i in range(obj_scatter_fits):
    sigma1 = [float(elem.find("sigma1", ns).text) for elem in obj_scatter_fits]  # (cm^-1)
    sigma2 = [float(elem.find("sigma2", ns).text) for elem in obj_scatter_fits]  # (cm^-1)
    B = [float(elem.find("B", ns).text) for elem in obj_scatter_fits]
    return sigma1, sigma2, B


def calculate_form_functions(thickness_map, sc_calib, grid_coords):
    sigma1, sigma2, B = get_form_func_params(sc_calib)
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


def calculate_edge_response(thickness_map, threshold=50, filter_size=25, num_iter=5):
    # TODO: fetch params from calibration.xml file
    edge_weight = thickness_map > threshold
    edge_weight = np.array(edge_weight, dtype="float")
    tmp_mask = edge_weight
    avg_kernel = np.ones([filter_size, filter_size]) / filter_size

    for i in range(num_iter):
        edge_weight = convolve2d(edge_weight, avg_kernel, mode="same")

    tmp = tmp_mask * edge_weight
    min_val = 0.6
    edge_weight = (1 - min_val) * (tmp - tmp.min()) / (tmp.max() - tmp.min()) + min_val
    return edge_weight


def get_ampl_params(sc_calib):
    ns = get_xmlns(sc_calib)
    obj_scatter_models = sc_calib.find("CalibrationResults/ObjectScatterModels", ns)
    obj_scatter_fits = obj_scatter_models.findall("ObjectScatterModel/ObjectScatterFit", ns)
    A = [float(elem.find("A", ns).text) for elem in obj_scatter_fits]
    alpha = [float(elem.find("alpha", ns).text) for elem in obj_scatter_fits]
    beta = [float(elem.find("beta", ns).text) for elem in obj_scatter_fits]
    return A, alpha, beta


def calculate_amplitudes(blank, proj, edge_weight, sc_calib):
    A, alpha, beta = get_ampl_params(sc_calib)
    num_groups = len(A)
    log_norm = log_norm_func(blank, proj)
    norm = proj / blank
    # TODO: improve thresholding with interpolation?
    norm[norm > 1] = 1
    norm[norm < 0] = 0
    amplitudes = []
    for i in range(num_groups):
        amplitude = A[i] * edge_weight * (norm ** alpha[i]) * (log_norm ** beta[i])
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
    u, v = get_detector_coords_cm(geometry)
    U, V = np.meshgrid(u, v)

    # Downsampled coords
    du, dv = get_detector_coords_cm(geometry, downsample_rate=downsample)
    DU, DV = np.meshgrid(du, dv)

    ns = get_xmlns(sc_calib)
    mu_water = float(sc_calib.find("CalibrationResults/Globals/muH2O", ns).text)  # mm^-1
    # for thickness smoothing
    sigma_u = float(sc_calib.find("CalibrationResults/Globals/AsymPertSigmaMMu", ns).text)
    sigma_v = float(sc_calib.find("CalibrationResults/Globals/AsymPertSigmaMMv", ns).text)
    step_du = np.mean(np.diff(du))
    step_dv = np.mean(np.diff(dv))
    sigma = (mm2cm(sigma_v) / step_dv, mm2cm(sigma_u) / step_du)

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
                thickness_map,
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
        # TODO: add nan check
        primaries[i] = calculate_primary(proj, scatter_est)
    return primaries


def read_dps_params(sc_calib):
    ns = get_xmlns(sc_calib)
    det_scatter_model = sc_calib.find("CalibrationResults/Globals/DetectorScatterModel", ns)
    return [float(elem.text) for elem in det_scatter_model]


def calculate_dps_kernel(sc_calib, U, V):
    grid = np.sqrt(U**2 + V**2)
    a = read_dps_params(sc_calib)
    det_kernel = a[0] * np.exp(-a[1] * grid) + a[2] * (np.exp(-a[3] * (grid - a[4]) ** 3))
    det_kernel = a[5] * det_kernel / np.sum(det_kernel)
    return det_kernel


def mm2cm(x):
    return 0.1 * x


def get_detector_coords_cm(geometry, downsample_rate=0):
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

    u = mm2cm(u)
    v = mm2cm(v)
    if downsample_rate:
        u = decimate(u, downsample_rate)
        v = decimate(v, downsample_rate)

    return u, v


def correct_detector_scatter(projs, geometry, sc_calib, downsample=8):

    u, v = get_detector_coords_cm(geometry)
    U, V = np.meshgrid(u, v)
    du, dv = get_detector_coords_cm(geometry, downsample_rate=downsample)
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


def interp_weight(x, xp, N=360):
    # linear interpolation of xp at x, mod(N).
    # xp must be monotonically increasing.
    i_min = np.argmin(abs(x - xp))
    if x - xp[i_min] >= 0:
        i_lower = i_min
    else:
        i_lower = (i_min - 1) % len(xp)

    i_upper = (i_lower + 1) % len(xp)

    xp_diff = (xp[i_upper] - xp[i_lower]) % N
    weights = (np.array([(xp[i_upper] - x), (x - xp[i_lower])]) % N) / float(xp_diff)
    return i_lower, i_upper, weights


def interpolate_blank_scan(angle, blank_projs, blank_angles, blank_airnorms):
    i_lower, i_upper, w = interp_weight(angle, blank_angles)
    blank_interp = w[0] * blank_projs[i_lower] + w[1] * blank_projs[i_upper]
    airnorm_interp = w[0] * blank_airnorms[i_lower] + w[1] * blank_airnorms[i_upper]
    return blank_interp, airnorm_interp


def log_normalize(projs, angles, airnorms, blank_projs, blank_angles, blank_airnorms):
    log_projs = np.zeros_like(projs)
    eps = np.finfo(projs.dtype).eps
    if blank_angles.size == 0:
        for i in range(len(angles)):
            cf_air = airnorms[i] / blank_airnorms
            ratio = cf_air * blank_projs / (projs[i] + eps)
            ratio[ratio < 1] = 1
            log_projs[i] = np.log(ratio)
    else:
        for i in range(len(angles)):
            blank_interp, airnorm_interp = interpolate_blank_scan(
                angles[i], blank_projs, blank_angles, blank_airnorms
            )
            cf_air = airnorms[i] / airnorm_interp
            ratio = cf_air * blank_interp / (projs[i] + eps)
            ratio[ratio < 1] = 1
            log_projs[i] = np.log(ratio)
    return log_projs


def load_blank_projections(filepath, scan_params):
    ns = get_xmlns(scan_params)
    start_angle = float(scan_params.find("Acquisitions/StartAngle", ns).text)
    stop_angle = float(scan_params.find("Acquisitions/StopAngle", ns).text)

    # Determine rotation direction
    if stop_angle - start_angle:
        rot_dir = "CC"
    else:
        rot_dir = "CW"

    # Get Truebeam version
    ALLOWED_VERSIONS = [2.0, 2.7]
    version_str = scan_params.attrib["Version"].rstrip(".0.0")
    try:
        version = float(version_str)
    except:
        print(f"Error retrieving version number, value found was '{version_str}'")
    else:

        if version not in ALLOWED_VERSIONS:
            raise ValueError(f"Version {version} not in allowed versions: {ALLOWED_VERSIONS}.")

    bowtie = scan_params.find("Acquisitions/Bowtie", ns)

    if bowtie is None:
        blank_fname = "Filter.xim"
    else:
        rot_str = f"_{rot_dir}*" if version == 2.7 else ""
        blank_fname = f"FilterBowtie{rot_str}.xim"

    blank_filepath = glob.glob(
        os.path.join(filepath, "Calibrations", "Air-*", "**", blank_fname), recursive=True
    )

    GANTRY_KVSOURCE_ANGLE_SEP = 90
    if version == 2.0:
        xim_img = XIM(blank_filepath[0])
        try:
            blank_projection = xim_img.array
        except AttributeError:
            print("Blank image not found.")
        else:
            blank_airnorm = xim_img.properties["KVNormChamber"]
            blank_projection = np.fliplr(np.array(blank_projection, dtype="float"))
            return blank_projection, np.array([]), float(blank_airnorm)

    elif version == 2.7:
        blank_projections = []
        blank_angles = []
        blank_airnorms = []
        for xim_filepath in blank_filepath:
            xim_img = XIM(xim_filepath)
            try:
                blank = xim_img.array
            except AttributeError:
                pass
            else:
                blank = np.fliplr(np.array(blank, dtype="float"))
                blank_projections.append(blank)
                blank_angles.append(xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP)
                blank_airnorms.append(xim_img.properties["KVNormChamber"])

        blank_projections = np.array(blank_projections)
        blank_angles = np.array(blank_angles)
        blank_airnorms = np.array(blank_airnorms, dtype="float")

        # Sort for easy interpolation later
        N = 360  # deg
        i_sort = np.argsort(blank_angles % N)
        blank_angles = blank_angles[i_sort] % N
        blank_projections = blank_projections[i_sort]
        blank_airnorms = blank_airnorms[i_sort]

        return blank_projections, blank_angles, blank_airnorms


def load_projections(filepath, threshold=0):

    # separation between the gantry head and the KV x-ray tube source (in deg)
    GANTRY_KVSOURCE_ANGLE_SEP = 90

    folder = os.path.join(filepath, "Acquisitions")

    ximfilelist = []
    for path in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, path)):
            subdir = os.path.join(folder, path)
            for file in os.listdir(subdir):
                if file.startswith("Proj_") and file.endswith(".xim"):
                    ximfilelist.append(os.path.join(subdir, file))

    projections = []
    angles = []
    airnorms = []

    print("Loading Varian CBCT dataset: " + folder)
    for xim_filepath in tqdm(ximfilelist):
        xim_img = XIM(xim_filepath)
        try:
            proj = xim_img.array
        except AttributeError:
            pass
        else:
            angle = xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP
            if not angles or abs(angle - angles[-1]) > threshold:
                proj = np.fliplr(np.array(proj, dtype="float"))
                projections.append(proj)
                angles.append(angle)
                airnorms.append(xim_img.properties["KVNormChamber"])

    projections = np.array(projections, dtype="float")
    angles = np.array(angles)
    airnorms = np.array(airnorms, dtype="float")
    return projections, angles, airnorms


def parse_inputs(**kwargs):
    """
    Returns tags.
    ACDC: acceleration-deceleration correction (default: True)
    DPS: detector point scatter correction (default: True)
    SC: kernel-based scatter correction (default: True)
    """

    acdc = kwargs["acdc"] if "acdc" in kwargs else True
    dps = kwargs["dps"] if "dps" in kwargs else True
    sc = kwargs["sc"] if "sc" in kwargs else True

    return acdc, dps, sc
