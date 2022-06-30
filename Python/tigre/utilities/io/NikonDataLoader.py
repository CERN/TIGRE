from __future__ import print_function
from __future__ import with_statement

import os
import math
import numpy
from tqdm import tqdm

from PIL import Image
from configparser import ConfigParser
from tigre.utilities.geometry import Geometry


def NikonDataLoader(filepath, **kwargs):
    # NikonDataLoader(filepath) Loads Nikon uCT datasets into TIGRE standard
    #
    # NikonDataLoader(filepath, OPT=VAL, ...) uses options and values.
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

    folder, geometry, angles = readXtekctGeometry(filepath)
    return loadNikonProjections(folder, geometry, angles, **kwargs)


def readXtekctGeometry(filepath):

    # Developed by A. Biguri and W. Sun
    # W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)

    if filepath.endswith(".xtekct"):
        folder, ini = os.path.split(filepath)
    else:
        folder = filepath
        files = [file for file in os.listdir(folder) if file.endswith(".xtekct")]
        if not files:
            raise ValueError("No .xtekct file found in folder: " + folder)
        ini = files[0]

    cfg = ConfigParser()
    cfg.read(os.path.join(folder, ini))
    cfg = cfg["XTekCT"]

    geometry = Geometry()
    geometry.accuracy = 0.5

    ## Detector information
    # Number of pixel in the detector
    geometry.nDetector = numpy.array((int(cfg["DetectorPixelsX"]), int(cfg["DetectorPixelsY"])))
    # Size of pixels in the detector
    geometry.dDetector = numpy.array(
        (float(cfg["DetectorPixelSizeX"]), float(cfg["DetectorPixelSizeY"]))
    )
    # Total size of the detector
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    ## Offset of the detector:
    geometry.offDetector = numpy.array(
        (float(cfg["DetectorOffsetX"]), float(cfg["DetectorOffsetY"]))
    )

    ## Image information
    # Number of voxels for the volume
    # the algos require these to be integers.
    geometry.nVoxel = numpy.array((int(cfg["VoxelsX"]), int(cfg["VoxelsY"]), int(cfg["VoxelsZ"])))
    # Size of each voxel
    geometry.dVoxel = numpy.array(
        (float(cfg["VoxelSizeX"]), float(cfg["VoxelSizeY"]), float(cfg["VoxelSizeZ"]))
    )
    # Size of the image in mm
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    geometry.offOrigin = numpy.array((0, 0, 0))

    #%% Global geometry
    geometry.DSO = float(cfg["SrcToObject"])
    geometry.DSD = float(cfg["SrcToDetector"])
    geometry.COR = -float(cfg["CentreOfRotationTop"])

    if geometry.COR == 0:
        print(
            "Centre of Rotation seems to be zero. Make sure that it is true and that the machine did not omit that information."
        )
    else:
        print(
            "TIGRE doesn't know if the sign of COR is the right one. Consider trying both and reporting to tigre.toolbox@gmail.com."
        )

    ## whitelevel
    geometry.whitelevel = float(cfg["WhiteLevel"])

    ## angles
    angles = []

    # It can be either an .ang  or .txt file
    # .ang
    files = [file for file in os.listdir(folder) if file.endswith(".ang")]
    if files:
        with open(os.path.join(folder, files[0]), "r") as file:
            # TODO: detect header lines automatically
            file.readline()
            for line in file:
                angles.append(math.radians(line.split("\t")[1]))
        return folder, geometry, numpy.array(angles)

    # .txt
    files = [file for file in os.listdir(folder) if file.endswith("_ctdata.txt")]
    if files:
        with open(os.path.join(folder, files[0]), "r") as file:
            # TODO: detect header lines automatically
            for i in range(3):
                file.readline()
            for line in file:
                angles.append(math.radians(float(line.split("\t")[1])))
        return folder, geometry, numpy.array(angles)

    print("File with definition of angles not found, estimating them from geometry info.")
    n_angles = int(cfg["Projections"])
    angle_step = float(cfg["AngularStep"])
    initial_angle = float(cfg["InitialAngle"])
    angles = numpy.arange(n_angles) * math.radians(angle_step) + math.radians(initial_angle)

    return folder, geometry, angles


def loadNikonProjections(folder, geometry, angles, **kwargs):
    # [projections, geometry, angles] = loadNikonProjections(filepath, geometry, angles, **kwargs)
    #    loads Nikon uCT machine projections
    #
    #    loadNikonData(filepath, geometry, angles) Loads a dataset given its FILEPATH,
    #       GEOMETRY and ANGLES (loaded from readXtekctGeometry())
    #       See NikonDataLoader() for additional options.
    #
    #% developed by A. Biguri and W. Sun 06.07.2020

    # parse inputs
    angles, indices = parse_inputs(geometry, angles, **kwargs)

    # load images
    files = sorted([file for file in os.listdir(folder) if file.lower().endswith(".tif")])

    image = Image.open(os.path.join(folder, files[indices[0]]))
    image = numpy.asarray(image).astype(numpy.float32)
    projections = numpy.zeros([len(indices), image.shape[0], image.shape[1]], dtype=numpy.single)
    projections[0, :, :] = -numpy.log(image / float(geometry.whitelevel))
    # Python uses zero-based indexing
    # if we start with index = 1 we get an error
    # use enumerate, it's cleaner
    print("Loading Nikon dataset: " + folder)
    for index, i in enumerate(tqdm(indices[1:]), 1):
        image = Image.open(os.path.join(folder, files[i]))
        image = numpy.asarray(image).astype(numpy.float32)
        projections[index, :, :] = -numpy.log(image / float(geometry.whitelevel))

    del geometry.whitelevel

    return numpy.asarray(projections), geometry, angles


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


if __name__ == "__main__":
    projections, geometry, angles = NikonDataLoader(
        "NikonDataLoaderTest/", sampling="continuous", num_angles=100
    )
    print("", projections.shape, "", geometry, "", angles, "", sep="\n")
