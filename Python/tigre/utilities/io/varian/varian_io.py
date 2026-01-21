from __future__ import annotations
import numpy as np
import os
from numpy.typing import NDArray
from Python.tigre.utilities.io.varian.utils import (
    XML,
    PathLike,
    XMLReader,
    interp_weight,
    sort_mod_N,
)
import xml.etree.ElementTree as ET
from Python.tigre.utilities.geometry import Geometry
from Python.tigre.utilities.io.varian.xim import XIM
from dataclasses import dataclass
import glob
from tqdm import tqdm

# separation between the gantry head and the KV x-ray tube source (in deg)
GANTRY_KVSOURCE_ANGLE_SEP = 90
# Varian TrueBeam versions
ALLOWED_VERSIONS = [2.0, 2.7]
# CC=Counterclockwise, CW=Clockwise
ROTATION_DIR = ["CC", "CW"]
SCAN_FILE = "Scan.xml"
RECON_FILE = "Reconstruction.xml"


@dataclass
class ProjData:
    projs: NDArray
    angles: NDArray

    def num_projs(self) -> int:
        return np.shape(self.projs)[0]

    def interp_proj(self, angle) -> NDArray:
        """
        Interpolate between projections at a given angle.
        If self.num_projs=1, the projection is returned unchanged.
        """
        if self.num_projs == 0:
            raise RuntimeError("Interpolation failed. No projections found.")
        if self.num_projs == 1:
            return self.projs
        else:
            i_lower, i_upper, w = interp_weight(angle, self.angles)
            blank_interp = w[0] * self.projs[i_lower] + w[1] * self.projs[i_upper]
            return blank_interp


def _read_scan_xml(filepath: PathLike) -> ET.Element:
    file = os.path.join(filepath, SCAN_FILE)
    tree = ET.parse(file)
    root = tree.getroot()
    if root is None:
        raise RuntimeError("Cone beam scan parameters not found.")
    else:
        return root


def _read_recon_xml(filepath: PathLike) -> ET.Element:
    file = glob.glob(os.path.join(filepath, "**", RECON_FILE), recursive=True)
    if len(file) > 1:
        raise RuntimeError(f"Multiple {RECON_FILE} files found.")
    tree = ET.parse(file[0])
    return tree.getroot()


class ScanParams(XML):
    def __init__(self, filepath: PathLike, xml_reader: XMLReader = _read_scan_xml) -> None:
        super().__init__(filepath, xml_reader)
        self.version: float = self._get_varian_version()
        self.root: ET.Element = self._get_acq_params()  # override root
        self.SID: float = self._get_SID()
        self.SAD: float = self._get_SAD()
        self.imager_size: list[int] = self._get_imager_size()
        self.imager_res: list[float] = self._get_imager_res()
        self.imager_lat: float = self._get_imager_lat()
        self.start_angle: float = self._get_start_angle()
        self.stop_angle: float = self._get_stop_angle()
        self.bowtie_filter: ET.Element | None = self._get_bowtie_filter()
        self.rot_velocity: float = self._get_rot_velocity()
        self.frame_rate: float = self._get_frame_rate()

    def _get_acq_params(self) -> ET.Element:
        acquisition_params = self.root.find("Acquisitions", self.ns)
        if acquisition_params is None:
            raise RuntimeError("Acquistion parameters not found in Scan.xml")
        else:
            return acquisition_params

    def _get_SID(self) -> float:
        field = "SID"
        elem = self._get_field(field)
        try:
            SID = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return SID

    def _get_SAD(self) -> float:
        field = "SAD"
        elem = self._get_field(field)
        try:
            SAD = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return SAD

    def _get_imager_size(self) -> list[int]:
        fields = ["ImagerSizeX", "ImagerSizeY"]
        elems = [self._get_field(field) for field in fields]
        try:
            imager_size = [int(elem.text) for elem in elems]
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return imager_size

    def _get_imager_res(self) -> list[float]:
        fields = ["ImagerResX", "ImagerResY"]
        res_elem = [self._get_field(field) for field in fields]
        try:
            imager_res = [float(elem.text) for elem in res_elem]
        except Exception:
            raise ValueError(f"Error converting {fields} to float")
        else:
            return imager_res

    def _get_imager_lat(self) -> float:
        field = "ImagerLat"
        elem = self._get_field(field)
        try:
            lat = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return lat

    def _get_start_angle(self) -> float:
        field = "StartAngle"
        elem = self._get_field(field)
        try:
            start_angle = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return start_angle

    def _get_stop_angle(self) -> float:
        field = "StopAngle"
        elem = self._get_field(field)
        try:
            stop_angle = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return stop_angle

    def _get_bowtie_filter(self) -> ET.Element | None:
        return self._get_field("Bowtie")

    def _get_rot_velocity(self) -> float:
        field = "Velocity"
        elem = self._get_field(field)
        try:
            velocity = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return velocity

    def _get_frame_rate(self) -> float:
        field = "FrameRate"
        elem = self._get_field(field)
        try:
            frame_rate = float(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return frame_rate

    def _get_varian_version(self) -> float:
        field = "Version"
        version_str = self.root.attrib[field].rstrip(".0.0")
        try:
            version = float(version_str)
        except Exception:
            raise ValueError(f"Error retrieving {field}. Found '{version_str}'")
        else:
            if version not in ALLOWED_VERSIONS:
                raise ValueError(f"Invalid {field}: {version}. Must be one of {ALLOWED_VERSIONS}.")
            else:
                return version

    def calculate_angular_threshold(self, frac: float = 0.9) -> float:
        return frac * (self.rot_velocity / self.frame_rate)

    def rot_direction(self) -> str:
        if self.stop_angle - self.start_angle:
            return ROTATION_DIR[0]
        else:
            return ROTATION_DIR[1]


class ReconParams(XML):
    def __init__(self, filepath: PathLike, xml_reader: XMLReader = _read_recon_xml) -> None:
        super().__init__(filepath, xml_reader)
        if self.root is not None:
            self.VOI_size = self._get_VOI_size()
            self.matrix_size = self._get_matrix_size()
            self.slice_thickness = self._get_slice_thickness()

    def _get_VOI_size(self) -> list[float]:
        fields = ["VOISizeX", "VOISizeY", "VOISizeZ"]
        elems = [self._get_field(field) for field in fields]
        try:
            VOI_size = [float(elem.text) for elem in elems]
        except Exception:
            raise ValueError(f"Invalid {fields}.")
        else:
            return VOI_size

    def _get_matrix_size(self) -> int:
        field = "MatrixSize"
        elem = self._get_field(field)
        try:
            matrix_size = int(elem.text)
        except Exception:
            raise ValueError(f"Invalid {field}.")
        else:
            return matrix_size

    def _get_slice_thickness(self) -> float:
        field = "SliceThickness"
        elem = self._get_field(field)
        try:
            slice_thickness = float(elem.text)
        except Exception:
            raise RuntimeError(f"Invalid {field}.")
        else:
            return slice_thickness


def read_varian_geometry(scan_params: ScanParams, recon_params: ReconParams | None) -> Geometry:

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


def air_normalize(proj: NDArray, airnorm: float) -> NDArray:
    return proj / airnorm


def load_blank_projections(filepath: PathLike, scan_params: ScanParams) -> ProjData:

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

    if scan_params.version not in ALLOWED_VERSIONS:
        raise ValueError(
            f"Version {scan_params.version} not in allowed versions. Must be one of {ALLOWED_VERSIONS}."
        )
    elif scan_params.version == ALLOWED_VERSIONS[0]:
        xim_img = XIM(blank_filepath[0])
        try:
            blank_proj = xim_img.array
        except AttributeError:
            raise RuntimeError(f"Blank scan image not found in {blank_filepath}.")
        else:
            blank_airnorm = float(xim_img.properties["KVNormChamber"])
            blank_proj = np.fliplr(np.array(blank_proj, dtype="float"))
            blank_proj = air_normalize(blank_proj, blank_airnorm)
            return ProjData(projs=blank_proj, angles=np.array([]))

    elif scan_params.version == ALLOWED_VERSIONS[1]:
        blank_projs = []
        blank_angles = []
        blank_airnorms = []
        for xim_filepath in blank_filepath:
            xim_img = XIM(xim_filepath)
            try:
                blank_proj = xim_img.array
            except AttributeError:
                pass
            else:
                blank_proj = np.fliplr(np.array(blank_proj, dtype="float"))
                blank_projs.append(blank_proj)
                blank_angles.append(
                    float(xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP)
                )
                blank_airnorms.append(float(xim_img.properties["KVNormChamber"]))

        blank_angles = np.array(blank_angles, dtype="float")
        blank_projs = np.array(
            [air_normalize(blank_projs[i], air) for i, air in enumerate(blank_airnorms)]
        )

        blank_angles, i_sort = sort_mod_N(blank_angles, N=360)
        blank_projs = blank_projs[i_sort]

        return ProjData(projs=blank_projs, angles=blank_angles)


def load_projections(filepath: PathLike, threshold: float = 0.0) -> ProjData:

    folder = os.path.join(filepath, "Acquisitions")

    ximfilelist = []
    for path in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, path)):
            subdir = os.path.join(folder, path)
            for file in os.listdir(subdir):
                if file.startswith("Proj_") and file.endswith(".xim"):
                    ximfilelist.append(os.path.join(subdir, file))

    projs = []
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
            angle = float(xim_img.properties["GantryRtn"] + GANTRY_KVSOURCE_ANGLE_SEP)
            if not angles or abs(angle - angles[-1]) > threshold:
                proj = np.fliplr(np.array(proj, dtype="float"))
                projs.append(proj)
                angles.append(angle)
                airnorms.append(float(xim_img.properties["KVNormChamber"]))

    projs = np.array([air_normalize(projs[i], air) for i, air in enumerate(airnorms)])

    return ProjData(projs=projs, angles=np.array(angles))
