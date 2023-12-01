from __future__ import print_function
from __future__ import with_statement

import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy

from tigre.utilities.geometry import Geometry


def DiondoDataLoader(filepath, **kwargs):
    #DiondoDataLoader(filepath) Loads Diondo uCT datasets into TIGRE standard
    #
    # DiondoDataLoader(filepath, OPT=VAL, ...) uses options and values.
    #    These are options in case you don't want to load the entire
    #    dataset, but only particular sets of projections.
    #    The possible arguments are:
    #         'sampling': type of sampling. default 'equidistant' Can be:
    #                'equidistant': equidistantly sample the entire set
    #                     of angles. 'num_angles' should be set for partial
    #                     loading of the data.
    #                'step': sample the entire set of projections every
    #                       'sampling_step' angles.
    #                'continous': Load the first 'num_angles' amount of
    #                           angles only.
    #
    #         'num_angles': Number of total angles to load. Default all of
    #                  them. Useful for 'equidistant' and 'continous' loading
    #
    #         'sampling_step': step to load when loading projections.
    #                 Default=1. Useful for 'step' loading.

    folder, geometry, angles, height = readDiondoGeometry(filepath)
    return loadDiondoProjections(folder, geometry, angles, height, **kwargs)


def readDiondoGeometry(filepath):
    # P. Basford
    # Based on Nikon data loader developed by A. Biguri and W. Sun
    # W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)
    filepath = Path(filepath)
    if filepath.suffix == ".xml":
        folder = filepath.parent
        config = filepath
    else:
        folder = filepath
        files = list(folder.glob("*.xml"))
        if not files:
            raise ValueError(f"No .xml file found in folder: {folder}")
        config = files[0]

    xml_tree = ET.parse(config)
    xml_root = xml_tree.getroot()

    xml_scan_params = xml_root.find("ScanParameter")
    if xml_scan_params.find("HelixMode").text == "true":
        raise NotImplementedError("Not yet written Helical")
    if xml_scan_params.find("Laminography").text == "true":
        raise NotImplementedError("Not yet written Laminography")

    geometry = Geometry()
    geometry.accuracy = 0.5

    ## Detector information
    xml_recon = xml_root.find("Recon")
    pixels_x = int(xml_recon.find("ProjectionDimX").text)
    pixels_y = int(xml_recon.find("ProjectionDimY").text)
    pixel_size_x = float(xml_recon.find("ProjectionPixelSizeX").text)
    pixel_size_y = float(xml_recon.find("ProjectionPixelSizeY").text)
    # Number of pixel in the detector
    geometry.nDetector = numpy.array((pixels_y, pixels_x))
    # Size of pixels in the detector
    geometry.dDetector = numpy.array(
        (pixel_size_y, pixel_size_x)
    )
    # Total size of the detector
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    ## Offset of the detector:
    geometry.offDetector = numpy.array(
        (0,0)
    )

    ## Image information
    voxels_x = int(xml_recon.find("VolumeDimX").text)
    voxels_y = int(xml_recon.find("VolumeDimY").text)
    voxels_z = int(xml_recon.find("VolumeDimZ").text)
    voxel_size_x = float(xml_recon.find("VolumeVoxelSizeX").text)
    voxel_size_y = float(xml_recon.find("VolumeVoxelSizeY").text)
    voxel_size_z = float(xml_recon.find("VolumeVoxelSizeZ").text)
    # Number of voxels for the volume
    # the algos require these to be integers.
    geometry.nVoxel = numpy.array(
        (voxels_z, voxels_y, voxels_x)
    )
    # Size of each voxel
    geometry.dVoxel = numpy.array(
        (voxel_size_z, voxel_size_y, voxel_size_x)
    )
    # Size of the image in mm
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    geometry.offOrigin = numpy.array((0, 0, 0))

    #%% Global geometry
    xml_geometry = xml_root.find("Geometrie")
    src_det = float(xml_geometry.find("SourceDetectorDist").text)
    src_obj = float(xml_geometry.find("SourceObjectDist").text)
    height = float(xml_geometry.find("Slices/string").text.split(",")[0].strip("mm"))
    projections = int(xml_recon.find("ProjectionCount").text)
    cor = float(xml_recon.find("ProjectionCenterOffsetX").text)
    geometry.DSO = src_obj
    geometry.DSD = src_det
    geometry.COR = -cor

    if geometry.COR == 0:
        print(
            "Centre of Rotation seems to be zero. Make sure that it is true and that the machine did not omit that information."
        )
    else:
        print(
            "TIGRE doesn't know if the sign of COR is the right one. Consider trying both and reporting to tigre.toolbox@gmail.com."
        )

    ## whitelevel
    geometry.whitelevel = float(2**16-1)

    ## angles
    angular_step = -math.radians(360.0/projections) # -ve due to the way it rotates
    angles = numpy.empty(projections)
    for i in range(0, projections):
        angles[i] = i * angular_step

    return folder, geometry, angles, height


def loadDiondoProjections(folder, geometry, angles, height, **kwargs):
    # [projections, geometry, angles] = loadDiondoProjections(filepath, geometry, angles, height, **kwargs)
    #    loads Diondo uCT machine projections
    #
    #       See DiondoDataLoader() for additional options.
    #
    #  Developed by P. Basford  07/2022
    #  Based on work by A. Biguri and W. Sun 06.07.2020

    # parse inputs
    angles, indices = parse_inputs(geometry, angles, **kwargs)

    # load images
    projection_folder = Path(folder, "Projection")
    if not projection_folder.exists():
        raise ValueError("Not importing the raw data from a diondo scan")
    filename_prefix=f"_{height:07.2f}".replace(".",",")
    (pixels_y, pixels_x) = geometry.nDetector
    projection_data = numpy.zeros([len(indices), pixels_y, pixels_x], dtype="<u2")
    whitelevel = 
    for i in tqdm(indices):
        filename = f"{filename_prefix}_{i:04d}.raw"
        file_path = Path(projection_folder, filename)
        image = numpy.fromfile(file_path, dtype="uint16") #read data in as 1D array
        image = image.reshape((pixels_y, pixels_x)) # reshape to be image size
        projection_data[i, :, : ] = -numpy.log(image / float(geometry.whitelevel))
    del geometry.whitelevel

    return projection_data, geometry, angles


def parse_inputs(geometry, angles, **kwargs):

    # TODO: warn user about invalid options or values
    sampling = kwargs["sampling"] if "sampling" in kwargs else "equidistant"
    nangles = int(kwargs["num_angles"]) if "num_angles" in kwargs else len(angles)
    step = int(kwargs["sampling_step"]) if "sampling_step" in kwargs else 1

    indices = numpy.arange(0, len(angles))

    if sampling == "equidistant":
        step = int(round(len(angles) / nangles))
        indices = indices[::step]
        angles = angles[::step]
    elif sampling == "continuous":
        indices = indices[:nangles]
        angles = angles[:nangles]
    elif sampling == "step":
        indices = indices[::step]
        angles = angles[::step]
    else:
        raise ValueError("Unknown sampling type: " + str(sampling))
    return angles, indices

def main():
    projections, geometry, angles = DiondoDataLoader(
        "DiondoDataLoaderTest/", sampling="continuous", num_angles=100
    )
    print("", projections.shape, "", geometry, "", angles, "", sep="\n")

if __name__ == "__main__":
    main()
