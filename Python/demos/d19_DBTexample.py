#%% Demo 19: Reconstruction for DBT
#
# This script demonstrates simple DBT reconstruction using OS-SART and FDK
# Adapted from MATLAB/Demos/d19_DBTexample.m
#
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
#
# Copyright (c) 2024, University of Bath and
#                     CERN-European Organization for Nuclear Research
#                     All rights reserved.
#
# License:            Open Source under BSD.
#                     See the full license at
#                     https://github.com/CERN/TIGRE/blob/master/LICENSE
#
# Contact:            tigre.toolbox@gmail.com
# Codes:              https://github.com/CERN/TIGRE/
# Coded by:           Tim Cogan
# --------------------------------------------------------------------------
#%%Initialize
import tigre
import tigre.algorithms as algs
from tigre.utilities import sample_loader
from tigre.utilities.common_geometry import staticDetectorGeo
from tigre.utilities.geometry import Geometry
from tigre.utilities.geometry_default import ConeGeometryDefault

import numpy as np
from numpy import ndarray



def get_angles(tube_angle: int = 25, num_projections: int = 9) -> ndarray:
    angles = np.linspace(-tube_angle / 2, tube_angle / 2, num_projections)
    return np.radians(angles)


def get_geometry() -> ConeGeometryDefault:
    geo = tigre.geometry_default(high_resolution=False)
    pixel_mm, slice_mm = 0.1, 0.5
    row_pixels, col_pixels = 3064, 2396
    source_to_detector_mm = 660
    source_to_origin_mm = 620

    geo.DSD = source_to_detector_mm
    geo.DSO = source_to_origin_mm

    geo.nDetector = np.array(([row_pixels, col_pixels])) // 2
    geo.dDetector = np.array([pixel_mm, pixel_mm])
    geo.sDetector = geo.nDetector * geo.dDetector

    num_voxels_z, num_voxels_y, num_voxels_x = 128, 2054, 784
    geo.nVoxel = np.array([num_voxels_z, num_voxels_y, num_voxels_x]) // 2
    geo.sVoxel = np.array([num_voxels_z * slice_mm, num_voxels_y * pixel_mm, num_voxels_x * pixel_mm]) // 2
    geo.dVoxel = geo.sVoxel / geo.nVoxel

    airgap_mm = 22
    # These two can be also defined per angle auxiliary
    geo.offOrigin = np.array([geo.sVoxel[0] / 2 - (geo.DSD - geo.DSO) + airgap_mm, 0, geo.sVoxel[2]/2])
    geo.offDetector = np.array([0, geo.sDetector[1] / 2])

    # Variable to define accuracy of 'interpolated' projection.
    # It defines the amount of samples per voxel.
    # Recommended <=0.5 (vx/sample)
    geo.accuracy=0.5

    # Can also be defined per angle
    center_of_rotation_y_displacement_mm = 0.
    geo.COR = center_of_rotation_y_displacement_mm  # type: ignore

    # Rotation of the detector, by X,Y and Z axis respectively.
    # This can also be defined per angle.
    geo.rotDetector = np.array([0, 0, 0])
    geo.mode= "cone"  # Can be "parallel" as well
    return geo


def get_projections(geo: Geometry, angles: ndarray) -> ndarray:
    # Example of data (CT Head). This is not a true DBT, but it works as an
    # example in TIGRE. Make sure to use a true tomosynthesis data.
    head = sample_loader.load_head_phantom(geo.nVoxel)
    projections = tigre.Ax(head, geo, angles, 'interpolated')
    # NOTE:
    # If you use a true DBT projection, use the following lines to adapt your data to TIGRE CT:
    # max_pixel_val = 2**14 - 1
    # projections = -np.log(projections / max_pixel_val)
    return projections


def compare_reconstruction_methods(projections: ndarray, geo: Geometry, angles: ndarray) -> None:
    imgFDK = algs.fdk(projections, geo, angles)

    num_iterations = 2
    imgOSSART = algs.ossart(projections, geo, angles, num_iterations)

    tigre.plotimg(np.concatenate([imgFDK, imgOSSART], axis=1), dim="z")


def main() -> None:
    angles = get_angles()
    geo = get_geometry()
    geo = staticDetectorGeo(geo, angles)  # Adapt CT geo to DBT
    projections = get_projections(geo, angles)
    compare_reconstruction_methods(projections, geo, angles)


if __name__ == "__main__":
    main()
