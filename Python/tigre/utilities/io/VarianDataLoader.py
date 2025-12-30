import numpy as np
import os
import glob
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET
from Python.tigre.utilities.io.xim_io import XIM
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
        frame_rate = rot_velocity = float(scan_params.find("Acquisitions/FrameRate", ns).text)
        angular_threshold = calculate_angular_threshold(rot_velocity, frame_rate)

    projs, angles, airnorms = load_projections(filepath, angular_threshold)

    blank_projs, blank_angles, blank_airnorms = load_blank_projections(filepath, scan_params)

    if dps:
        projs = correct_detector_point_scatter(projs, geometry, sc_calib)
        blank_projs = correct_detector_point_scatter(blank_projs, geometry, sc_calib)

    if sc:
        projs = correct_scatter(
            sc_calib, blank_projs, blank_angles, blank_airnorms, projs, airnorms, geometry
        )

    log_projs = log_normalize(projs, angles, airnorms, blank_projs, blank_angles, blank_airnorms)

    return log_projs, geometry, angles


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
            float(acquisition_params.find("ImagerSizeX", ns).text),
            float(acquisition_params.find("ImagerSizeY", ns).text),
        ]
    )
    geometry.dDetector = np.array(
        [
            float(acquisition_params.find("ImagerResX", ns).text),
            float(acquisition_params.find("ImagerResY", ns).text),
        ]
    )
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    offset = float(acquisition_params.find("ImagerLat", ns).text)
    geometry.offDetector = np.array(
        [-1 * offset, 0]
    )  # offset orientation in VarianCBCT is opposite to TIGRE
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
            [geometry.dDetector[0], geometry.dDetector[0], geometry.dDetector[1]]
        ) * (geometry.DSO / geometry.DSD)
        geometry.nVoxel = np.ceil(
            np.array(
                [
                    geometry.nDetector[0] + abs(geometry.offDetector[0]) / geometry.dDetector[0],
                    geometry.nDetector[0] + abs(geometry.offDetector[0]) / geometry.dDetector[0],
                    geometry.nDetector[1],
                ]
            )
        )
        geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    else:
        # Retrieve existing reconstruction parameters
        geometry.sVoxel = np.array(
            [
                float(recon_params.find("VOISizeX", ns).text),
                float(recon_params.find("VOISizeY", ns).text),
                float(recon_params.find("VOISizeZ", ns).text),
            ]
        )
        slice_num = round(geometry.sVoxel[2] / float(recon_params.find("SliceThickness", ns).text))
        geometry.nVoxel = np.array(
            [
                float(recon_params.find("MatrixSize", ns).text),
                float(recon_params.find("MatrixSize", ns).text),
                slice_num,
            ]
        )
        geometry.dVoxel = geometry.sVoxel / geometry.nVoxel

    return geometry, scan_params


def calculate_angular_threshold(rot_velocity, frame_rate):
    angular_interval = rot_velocity / frame_rate
    return 0.9 * angular_interval


def read_scatter_calib(filepath):

    file = os.path.join(filepath, "Calibrations", "SC-*", "Factory", "Calibration.xml")
    tree = ET.parse(file)
    sc_calib = tree.getroot().find("Calibration")

    return sc_calib


def correct_scatter(sc_calib, blank, sec, blk_airnorm, projs, airnorm, geo):
    pass
    # TODO


def correct_detector_point_scatter(proj, geo, sc_calib):
    pass
    # TODO


def log_normalize(projs, angles, airnorm, blk, sec, blk_airnorm):
    pass
    # TODO


def remove_ring_artifacts(log_proj):
    pass
    # TODO


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
            return np.array(blank_projection, dtype="float32"), np.array([]), float(blank_airnorm)

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
                blank_projections.append(blank)
                blank_angles.append(xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP)
                blank_airnorms.append(xim_img.properties["KVNormChamber"])
        blank_projections = np.array(blank_projections, dtype="float32")
        blank_angles = np.array(blank_angles)
        blank_airnorms = np.array(blank_airnorms, dtype="float32")

        # TODO: sort to be monotonically increasing
        i_sort = np.argsort(blank_angles % 360)
        blank_angles = blank_angles[i_sort]
        blank_projections = blank_projections[i_sort,]
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
    for xim_filepath in tqdm(ximfilelist):
        xim_img = XIM(xim_filepath)
        try:
            proj = xim_img.array
        except AttributeError:
            pass
        else:
            angle = xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP
            if not angles or abs(angle - angles[-1]) > threshold:
                angles.append(angle)
                projections.append(proj)
                airnorms.append(xim_img.properties["KVNormChamber"])

    projections = np.array(projections, dtype="float32")
    angles = np.array(angles)
    airnorms = np.array(airnorms, dtype="float32")
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
