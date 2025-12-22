import numpy as np
import os
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET


def VarianDataLoader(filepath, **kwargs):
    # add kwargs for acdc, dps, and sc. Disregard BH correction.
    # parse_inputs
    folder, geometry, angles = readVarianGeometry(filepath)  # read from xml. equiv: GeometryFromXML
    sc_calib = readScatterCalib(filepath)  # read scatter calib from xml. equiv: ScCalibFromXML
    return loadVarianProjections(folder, geometry, angles, **kwargs)


def readVarianGeometry(filepath):
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
    file = os.path.join(filepath, "Reconstruction.xml")
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

    return folder, geometry, angles, scan_xml


def readScatterCalib(filepath):

    return sc_calib


def loadVarianProjections(folder, geometry, angles, **kwargs):
    # load projections
    # if tag_acdc: remove over-sampled projections
    # load blank scan
    # if tag_dps: perform detector point scatter correction
    # if tag_sc: perform asks scatter correction
    # perform log norm with airnorm
    # enforce positive values
    # ring removal
    # convert angles from deg to rad
    return np.asarray(projections), geometry, angles
