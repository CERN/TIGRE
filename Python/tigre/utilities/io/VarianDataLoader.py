import numpy as np
import os
import glob
from Python.tigre.utilities.geometry import Geometry
from Python.tigre.utilities.io.varian.utils import get_xmlns, interpolate_blank_scan
from Python.tigre.utilities.io.varian.scatter_correct import (
    read_scatter_calib,
    correct_detector_scatter,
    correct_scatter,
)
from Python.tigre.utilities.io.varian.xim_io import XIM
import xml.etree.ElementTree as ET
from scipy.ndimage import median_filter
from tqdm import tqdm


def VarianDataLoader(filepath, **kwargs):
    # TODO: add docstring
    acdc, dps, sc = parse_inputs(**kwargs)

    geometry, scan_params = read_varian_geometry(filepath)
    if dps or sc:
        sc_calib = read_scatter_calib(
            filepath
        )  # TODO: move to scatter correct module and pass filepath to scatter correction functions instead

    angular_threshold = 0
    if acdc:
        angular_threshold = calculate_angular_threshold(scan_params)
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
    log_projs = correct_ring_artifacts(log_projs)

    return log_projs, geometry, np.deg2rad(angles)


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


def _get_angular_threshold_params(scan_params):
    ns = get_xmlns(scan_params)
    rot_velocity = float(scan_params.find("Acquisitions/Velocity", ns).text)
    frame_rate = float(scan_params.find("Acquisitions/FrameRate", ns).text)
    return rot_velocity, frame_rate


def calculate_angular_threshold(scan_params):
    rot_velocity, frame_rate = _get_angular_threshold_params(scan_params)
    return 0.9 * (rot_velocity / frame_rate)


def correct_ring_artifacts(log_projs, kernel_size=(1, 9)):
    log_projs = np.array([median_filter(p, size=kernel_size) for p in tqdm(log_projs)])
    return log_projs


def log_normalize(projs, angles, airnorms, blank_projs, blank_angles, blank_airnorms):
    log_projs = np.zeros_like(projs)
    eps = np.finfo(projs.dtype).eps
    if blank_angles.size == 0:
        for i in range(len(angles)):
            cf_air = airnorms[i] / blank_airnorms
            ratio = cf_air * blank_projs / (projs[i] + eps)
            ratio[ratio < 1] = 1  # TODO: replace with
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
