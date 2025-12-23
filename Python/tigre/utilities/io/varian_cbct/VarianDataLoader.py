import numpy as np
import os
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET
from xim_io import XIM
from tqdm import tqdm


def VarianDataLoader(filepath, **kwargs):

    acdc, dps, sc = parse_inputs(**kwargs)  #

    folder, geometry = read_varian_geometry(filepath)  # read from xml. equiv: GeometryFromXML
    if dps or sc:
        sc_calib = read_scatter_calib(
            filepath
        )  # read scatter calib from xml. equiv: ScCalibFromXML

    angular_threshold = 0
    if acdc:

        angular_threshold = angular_interval * 0.9

    projs, airnorm = load_varian_projections(folder, geometry, **kwargs)

    return projections, geometry, angles


def read_varian_geometry(filepath):
    # read from xml
    geometry = Geometry()

    # read Scan.xml
    file = os.path.join(filepath, "Scan.xml")
    tree = ET.parse(file)
    scan_params = tree.getroot().find("Scan")
    if scan_params is None:
        raise RuntimeError("Cone beam scan parameters not found")

    # CBCT geo params
    geometry.mode = "cone"
    geometry.DSD = float(scan_params.Acquisitions.SID.text)
    geometry.DSO = float(scan_params.Acquisitions.SAD.text)
    geometry.nDetector = np.array(
        [
            float(scan_params.Acquisitions.ImagerSizeX.text),
            float(scan_params.Acquisitions.ImagerSizeY.text),
        ]
    )
    geometry.dDetector = np.array(
        [
            float(scan_params.Acquisitions.ImagerResX.text),
            float(scan_params.Acquisitions.ImagerResY.text),
        ]
    )
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    offset = float(scan_params.Acquisitions.ImagerLat.text)
    geometry.offDetector = np.array(
        [-1 * offset, 0]
    )  # offset orientation in VarianCBCT is opposite to TIGRE
    geometry.offOrigin = np.array([0, 0, 0])

    # Auxiliary
    # Variable to define accuracy of 'interpolated' projection
    # It defines the amoutn of samples per voxel.
    # Recommended <=0.5          (vx/sample)
    geometry.accuracy = 0.5

    # read Reconstruction.xml
    file = os.path.join(filepath, "**", "Reconstruction.xml")
    tree = ET.parse(file)
    recon_params = tree.getroot().find("Reconstruction")
    if recon_params is None:
        print("Estimating acceptable image size...")
        geometry.dVoxel = np.array(
            [geometry.dDetector(0), geometry.dDetector(0), geometry.dDetector(1)]
        ) * (geometry.DSO / geometry.DSD)
        geometry.nVoxel = np.ceil(
            np.array(
                [
                    geometry.nDetector(0) + abs(geometry.offDetector(0)) / geometry.dDetector(0),
                    geometry.nDetector(0) + abs(geometry.offDetector(0)) / geometry.dDetector(0),
                    geometry.nDetector(1),
                ]
            )
        )
        geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    else:
        pass
        # TODO

    return folder, geometry, scan_xml


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


def load_blank_projections(filepath):
    pass
    # TODO


def load_projections(filepath, geometry, angles, **kwargs):

    xim_files = os.dir(os.path.join(filepath, "**", "Proj_*.xim"))

    # load images
    xim_files = sorted(
        [file for file in os.listdir(folder) if file.startswith("Proj_") and file.endswith(".xim")]
    )

    return projections, angles, airnorm


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
