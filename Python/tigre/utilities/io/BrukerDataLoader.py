from __future__ import print_function
from __future__ import with_statement

import os
import math
import numpy
from tqdm import tqdm

from PIL import Image
from configparser import ConfigParser, Error
from tigre.utilities.geometry import Geometry


def BrukerDataLoader(filepath, **kwargs):
    # BrukerDataLoader(filepath) Loads Bruker Skyscan datasets into TIGRE standard
    #
    # BrukerDataLoader(filepath, OPT=VAL, ...) uses options and values.
    #    These are options in case you don't want to load the entire
    #    dataset, but only particular sets of projections.
    #    The possible arguments are:
    #         'dataset_number': Some fodlers will have several scans.
    #                 Set to 'all' to load all of them, or give a number (starting from 0)
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

    dataset_number = kwargs["dataset_number"] if "dataset_number" in kwargs else None
    if dataset_number == "all":
        print("Loading all scans in folder, assuming the same geometry for all \n\n")
        del kwargs["dataset_number"]
        num_scans = find_number_of_datasets(filepath)
        angles = [None] * num_scans
        projections = [None] * num_scans
        for scan in range(num_scans):
            print("Loading scan number " + str(scan) + "...\n")
            folder, geometry, angles[scan] = read_Bruker_geometry(
                filepath, dataset_number=scan, **kwargs
            )
            projections[scan], geometry, angles[scan] = load_Bruker_projections(
                folder, geometry, angles[scan], dataset_number=scan, **kwargs
            )
        projections = numpy.concatenate(projections)
        angles = numpy.concatenate(angles)

        return projections, geometry, angles

    else:

        folder, geometry, angles = read_Bruker_geometry(filepath, **kwargs)
        return load_Bruker_projections(folder, geometry, angles, **kwargs)


def read_Bruker_geometry(filepath, **kwargs):

    # check if input was log file itself, or just the folder
    if filepath.endswith(".log"):
        folder, ini = os.path.split(filepath)
    else:
        folder = filepath
        files = [file for file in os.listdir(folder) if file.endswith(".log")]
        if not files:
            raise ValueError("No .log file found in folder: " + folder)

        num_scans = find_number_of_datasets(filepath)

        if num_scans is not None:
            dataset_number = kwargs["dataset_number"] if "dataset_number" in kwargs else None
            if dataset_number is None:
                raise ValueError(
                    "This folder contains many datasets, please select which one to load with BrukerDataLoader(..., dataset_number=a_number)"
                )
            if dataset_number >= num_scans:
                raise ValueError("Dataset number given larger than total number of datasets")

            matching = [s for s in files if "{0:0=2d}".format(dataset_number) + ".log" in s]
            if len(matching) > 1:
                raise AssertionError(
                    "More than 1 file for the same dataset found, confused what to do, so I error"
                )
            ini = matching[0]
        else:
            ini = files[0]
    # create configureation parser
    cfg = ConfigParser()
    cfg.read(os.path.join(folder, ini))
    cfg_system = cfg["System"]

    # start empty geometry
    geometry = Geometry()
    geometry.accuracy = 0.5

    ## Detector information
    # Size of pixels in the detector

    try:
        ratio = float(cfg_system["CameraXYRatio"])
    except KeyError:
        ratio = float(cfg_system["Camera X/Y Ratio"])
    except:
        ratio = 1
    # need to verify if this is correct.
    ratio = 1

    cfg_aq = cfg["Acquisition"]
    try:
        binu = float(cfg_aq["Camera binning"][0])
        binv = float(cfg_aq["Camera binning"][2])
    except:
        binu = 1
        binv = 1

    geometry.dDetector = numpy.array(
        (
            float(cfg_system["Camera Pixel Size (um)"]) / 1000.0 * binv * ratio,
            float(cfg_system["Camera Pixel Size (um)"]) / 1000.0 * binu,
        )
    )
    # Number of pixel in the detector
    geometry.nDetector = numpy.array(
        (float(cfg_aq["Number of Rows"]), float(cfg_aq["Number of Columns"]))
    )

    # Total size of the detector
    geometry.sDetector = geometry.nDetector * geometry.dDetector
    geometry.nDetector = geometry.nDetector.astype(int)
    ## Offset of the detector:
    try:
        offs = (
            -(geometry.nDetector[0] / 2 - float(cfg_aq["Optical Axis (line)"]))
            * geometry.dDetector[0]
        )
    except:
        offs = 0.0

    geometry.offDetector = numpy.array((offs, 0.0))

    # Size of each voxel
    geometry.dVoxel = numpy.array(
        (
            float(cfg_aq["Image Pixel Size (um)"]) / 1000,
            float(cfg_aq["Image Pixel Size (um)"]) / 1000,
            float(cfg_aq["Image Pixel Size (um)"]) / 1000,
        )
    )
    geometry.nVoxel = numpy.array(
        (geometry.nDetector[0], geometry.nDetector[1], geometry.nDetector[1])
    )
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel
    geometry.nVoxel = geometry.nVoxel.astype(int)

    #% Global geometry
    geometry.DSO = float(cfg_aq["Object to Source (mm)"])
    geometry.DSD = float(cfg_aq["Camera to Source (mm)"])

    geometry.whitelevel = 2 ** int(cfg_aq["Depth (bits)"])

    angles = numpy.arange(
        0.0,
        float(cfg_aq["Number of Files"]) * float(cfg_aq["Rotation Step (deg)"]),
        float(cfg_aq["Rotation Step (deg)"]),
    )
    angles = angles[:-1] * numpy.pi / 180

    files = [file for file in os.listdir(folder) if file.endswith(".csv")]
    if files:
        offset = numpy.genfromtxt(
            os.path.join(folder, files[0]), delimiter=",", skip_header=5, dtype=float
        )
        offset = numpy.delete(offset, 0, 1)
        geometry.offDetector = geometry.offDetector + numpy.fliplr(offset) * (
            geometry.dDetector / [binv, binu]
        )

    return filepath, geometry, angles


def load_Bruker_projections(folder, geometry, angles, **kwargs):

    angles, indices, dataset_number = parse_inputs(geometry, angles, **kwargs)

    # load images
    files = sorted([file for file in os.listdir(folder) if file.lower().endswith(".tif")])
    if dataset_number is not None:
        files = [file for file in files if file[-10:-8] == "{0:0=2d}".format(dataset_number)]

    image = Image.open(os.path.join(folder, files[indices[0]]))
    image = numpy.asarray(image).astype(numpy.float32)
    projections = numpy.zeros([len(indices), image.shape[0], image.shape[1]], dtype=numpy.single)
    projections[0, :, :] = -numpy.log(image / float(geometry.whitelevel))
    index = 1

    print("Loading Bruker Skyscan dataset: " + folder)
    for i in tqdm(indices[1:]):
        image = Image.open(os.path.join(folder, files[i]))
        image = numpy.asarray(image).astype(numpy.float32)
        projections[index, :, :] = -numpy.log(image / float(geometry.whitelevel))
        index = index + 1
    del geometry.whitelevel

    return numpy.asarray(projections), geometry, angles


## This should be on a separate "io_common.py" file.
def parse_inputs(geometry, angles, **kwargs):

    # TODO: warn user about invalid options or values
    sampling = kwargs["sampling"] if "sampling" in kwargs else "equidistant"
    nangles = int(kwargs["num_angles"]) if "num_angles" in kwargs else len(angles)
    step = int(kwargs["sampling_step"]) if "sampling_step" in kwargs else 1
    dataset_number = kwargs["dataset_number"] if "dataset_number" in kwargs else None

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

    return angles, indices, dataset_number


def find_number_of_datasets(filepath):
    # check if input was log file itself, or just the folder
    if filepath.endswith(".log"):
        folder, ini = os.path.split(filepath)
    else:
        folder = filepath
        files = [file for file in os.listdir(folder) if file.endswith(".log")]
        if not files:
            raise ValueError("No .log file found in folder: " + folder)

        ini = min(files, key=len)  # shortest one is the main one?
        cfg = ConfigParser()
        cfg.read(os.path.join(folder, ini))
        cfg_aq = cfg["Acquisition"]

        if cfg.has_option("Acquisition", "Number of connected scans"):
            return int(cfg_aq["Number of connected scans"])
        else:
            return None
