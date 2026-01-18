from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
import os
import glob
from Python.tigre.utilities.geometry import Geometry
from Python.tigre.utilities.io.varian.utils import interp_weight, XML
from Python.tigre.utilities.io.varian.scatter_correct import (
    DetScattParams,
    ScattParams,
    read_scatt_calib_xml,
    correct_detector_scatter,
    correct_scatter,
)
from Python.tigre.utilities.io.varian.xim_io import XIM
import xml.etree.ElementTree as ET
from scipy.ndimage import median_filter
from tqdm import tqdm
from dataclasses import dataclass

# separation between the gantry head and the KV x-ray tube source (in deg)
GANTRY_KVSOURCE_ANGLE_SEP = 90
# Varian TrueBeam versions
ALLOWED_VERSIONS = [2.0, 2.7]


def VarianDataLoader(filepath, **kwargs):

    acdc, dps, sc = parse_inputs(**kwargs)

    scan_params = ScanParams(filepath, xml_reader=read_scan_xml)
    recon_params = ReconParams(filepath, xml_reader=read_recon_xml)
    geometry = read_varian_geometry(scan_params, recon_params)

    angular_threshold = 0.0
    if acdc:
        angular_threshold = scan_params.calculate_angular_threshold()

    proj_data = load_projections(filepath, angular_threshold)
    blank_proj_data = load_blank_projections(filepath, scan_params)

    if dps:
        dps_calib = DetScattParams(filepath, xml_reader=read_scatt_calib_xml)
        blank_proj_data.projs = correct_detector_scatter(blank_proj_data.projs, geometry, dps_calib)
        proj_data.projs = correct_detector_scatter(proj_data.projs, geometry, dps_calib)

    if sc:
        sc_calib = ScattParams(filepath, xml_reader=read_scatt_calib_xml)
        proj_data.projs = correct_scatter(proj_data, blank_proj_data, geometry, sc_calib)

    log_projs = log_normalize(proj_data, blank_proj_data)
    log_projs = correct_ring_artifacts(log_projs)

    return log_projs, geometry, np.deg2rad(proj_data.angles)


@dataclass
class ProjData:
    projs: NDArray
    angles: NDArray
    airnorms: NDArray

    def num_projs(self):
        return np.shape(self.projs)[0]

    def interp_proj(self, angle):
        """Interpolate between projections and airnorm values at a given angle.
        If self.num_projs=1, the projection and airnorm are returned unchanged."""
        if self.num_projs == 0:
            raise RuntimeError("Interpolation failed. No projections found.")
        if self.num_projs == 1:
            return self.projs, self.airnorms
        else:
            i_lower, i_upper, w = interp_weight(angle, self.angles)
            blank_interp = w[0] * self.projs[i_lower] + w[1] * self.projs[i_upper]
            airnorm_interp = w[0] * self.airnorms[i_lower] + w[1] * self.airnorms[i_upper]
            return blank_interp, airnorm_interp


class ScanParams(XML):
    def __init__(self, filepath, xml_reader):
        super().__init__(filepath, xml_reader)
        self.acquisition_params = self._get_acq_params()
        self.SID = self._get_SID()
        self.SAD = self._get_SAD()
        self.imager_size = self._get_imager_size()
        self.imager_res = self._get_imager_res()
        self.imager_lat = self._get_imager_lat()
        self.start_angle = self._get_start_angle()
        self.stop_angle = self._get_stop_angle()
        self.bowtie_filter = self._get_bowtie_filter()
        self.rot_velocity = self._get_rot_velocity()
        self.frame_rate = self._get_frame_rate()
        self.version = self._get_varian_version()

    def _get_acq_params(self):
        acquisition_params = self.root.find("Acquisitions", self.ns)
        return acquisition_params

    def _get_SID(self):
        return float(self.acquisition_params.find("SID", self.ns).text)

    def _get_SAD(self):
        return float(self.acquisition_params.find("SAD", self.ns).text)

    def _get_imager_size(self):
        imager_size = [
            int(self.acquisition_params.find("ImagerSizeX", self.ns).text),
            int(self.acquisition_params.find("ImagerSizeY", self.ns).text),
        ]
        return imager_size

    def _get_imager_res(self):
        imager_res = [
            float(self.acquisition_params.find("ImagerResX", self.ns).text),
            float(self.acquisition_params.find("ImagerResY", self.ns).text),
        ]
        return imager_res

    def _get_imager_lat(self):
        return float(self.acquisition_params.find("ImagerLat", self.ns).text)

    def _get_start_angle(self):
        return float(self.root.find("Acquisitions/StartAngle", self.ns).text)

    def _get_stop_angle(self):
        return float(self.root.find("Acquisitions/StopAngle", self.ns).text)

    def _get_bowtie_filter(self):
        return self.root.find("Acquisitions/Bowtie", self.ns)

    def _get_rot_velocity(self):
        return float(self.root.find("Acquisitions/Velocity", self.ns).text)

    def _get_frame_rate(self):
        return float(self.root.find("Acquisitions/FrameRate", self.ns).text)

    def _get_varian_version(self):
        version_str = self.root.attrib["Version"].rstrip(".0.0")
        try:
            version = float(version_str)
        except Exception:
            raise RuntimeError(f"Error retrieving version number. Value found was '{version_str}'")
        else:
            if version not in ALLOWED_VERSIONS:
                raise ValueError(f"Version {version} not in allowed versions: {ALLOWED_VERSIONS}.")
            else:
                return version

    def calculate_angular_threshold(self, frac=0.9):
        return frac * (self.rot_velocity / self.frame_rate)

    def rot_direction(self):
        if self.stop_angle - self.start_angle:
            return "CC"
        else:
            return "CW"


class ReconParams(XML):
    def __init__(self, filepath, xml_reader):
        super().__init__(filepath, xml_reader)
        self.VOI_size = self._get_VOI_size()
        self.matrix_size = self._get_matrix_size()
        self.slice_thickness = self._get_slice_thickness()

    def _get_VOI_size(self):
        VOI_size = [
            float(self.root.find("VOISizeX", self.ns).text),
            float(self.root.find("VOISizeY", self.ns).text),
            float(self.root.find("VOISizeZ", self.ns).text),
        ]
        return VOI_size

    def _get_matrix_size(self):
        return (int(self.root.find("MatrixSize", self.ns).text),)

    def _get_slice_thickness(self):
        return float(self.root.find("SliceThickness", self.ns).text)


def read_scan_xml(filepath):
    file = os.path.join(filepath, "Scan.xml")
    tree = ET.parse(file)
    root = tree.getroot()
    if root is None:
        raise RuntimeError("Cone beam scan parameters not found.")
    else:
        return root


def read_recon_xml(filepath):
    file = glob.glob(os.path.join(filepath, "**", "Reconstruction.xml"), recursive=True)
    if len(file) > 1:
        raise RuntimeError("Multiple Reconstruction.xml files found.")
    tree = ET.parse(file[0])
    return tree.getroot()


def read_varian_geometry(scan_params, recon_params):

    geometry = Geometry()
    geometry.mode = "cone"
    geometry.DSD = scan_params.SID
    geometry.DSO = scan_params.SAD
    geometry.nDetector = np.array(scan_params.imager_size[::-1])
    geometry.dDetector = np.array(scan_params.imager_res[::-1])
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    offset = scan_params.imager_lat
    geometry.offDetector = np.array([0, -offset])
    geometry.offOrigin = np.array([0, 0, 0])

    # Auxiliary variable to define accuracy of 'interpolated' projection
    # Recommended <=0.5 (vx/sample)
    geometry.accuracy = 0.5

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
        geometry.sVoxel = np.array(recon_params.VOI_size[::-1])
        num_slices = np.ceil(geometry.sVoxel[0] / recon_params.slice_thickness)
        geometry.nVoxel = np.array(
            [
                num_slices,
                recon_params.matrix_size,
                recon_params.matrix_size,
            ]
        )
        geometry.dVoxel = geometry.sVoxel / geometry.nVoxel

    return geometry


def correct_ring_artifacts(log_projs, kernel_size=(1, 9)):
    print("Performing ring artifact correction:")
    log_projs = np.array([median_filter(p, size=kernel_size) for p in tqdm(log_projs)])
    return log_projs


def log_normalize(proj_data, blank_proj_data):
    log_projs = np.zeros_like(proj_data.projs)
    eps = np.finfo(proj_data.projs.dtype).eps

    print("Performing log normalization:")
    for i in tqdm(range(len(proj_data.angles))):
        blank_interp, airnorm_interp = blank_proj_data.interp_proj(proj_data.angles[i])
        blank = blank_interp * proj_data.airnorms[i] / airnorm_interp

        ratio = (blank + eps) / (proj_data.projs[i] + eps)
        ratio[ratio < 1] = 1
        log_projs[i] = np.log(ratio)
    return log_projs


def load_blank_projections(filepath, scan_params):

    if scan_params.bowtie_filter is None:
        blank_fname = "Filter.xim"
    else:
        rot_str = (
            f"_{scan_params.rot_direction()}*" if scan_params.version == ALLOWED_VERSIONS[1] else ""
        )
        blank_fname = f"FilterBowtie{rot_str}.xim"

    blank_filepath = glob.glob(
        os.path.join(filepath, "Calibrations", "Air-*", "**", blank_fname), recursive=True
    )

    if scan_params.version == ALLOWED_VERSIONS[0]:
        xim_img = XIM(blank_filepath[0])
        try:
            blank_projection = xim_img.array
        except AttributeError:
            raise RuntimeError(f"blank_proj_data image not found in {blank_filepath}.")
        else:
            blank_airnorm = xim_img.properties["KVNormChamber"]
            blank_projection = np.fliplr(np.array(blank_projection, dtype="float"))
        return ProjData(
            projs=blank_projection, angles=np.array([]), airnorms=np.array(blank_airnorm)
        )

    elif scan_params.version == ALLOWED_VERSIONS[1]:
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

        blank_projections /= np.expand_dims(blank_airnorms, axis=(1, 2))
        blank_airnorms = np.ones_like(blank_airnorms)

        return ProjData(projs=blank_projections, angles=blank_angles, airnorms=blank_airnorms)


def load_projections(filepath, threshold=0.0):

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

    projections /= np.expand_dims(airnorms, axis=(1, 2))
    airnorms = np.ones_like(airnorms)

    return ProjData(projs=projections, angles=angles, airnorms=airnorms)


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
    if sc:
        RuntimeWarning("dps must be enabled when sc=True. Setting dps=True.")
        dps = True
    return acdc, dps, sc
