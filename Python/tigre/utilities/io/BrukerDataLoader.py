from __future__ import print_function
from __future__ import with_statement

import os
import math
import numpy
from tqdm import tqdm

from PIL import Image
from configparser import ConfigParser
from tigre.utilities.geometry import Geometry

def BrukerDataLoader(filepath, **kwargs):
    # BrukerDataLoader(filepath) Loads Bruker Skyscan datasets into TIGRE standard
    #
    # BrukerDataLoader(filepath, OPT=VAL, ...) uses options and values.
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

    folder, geometry, angles = read_Bruker_geometry(filepath)
    return load_Bruker_projections(folder, geometry, angles, **kwargs)

def read_Bruker_geometry(filepath):
    
    # check if input was log file itself, or just the folder
    if filepath.endswith(".log"):
        folder, ini = os.path.split(filepath)
    else:
        folder = filepath
        files = [file for file in os.listdir(folder) if file.endswith(".log")]
        if not files:
            raise ValueError("No .log file found in folder: " + folder)
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

    geometry.dDetector = numpy.array(
        (float(cfg_system["Camera Pixel Size (um)"])/1000.0, float(cfg_system["Camera Pixel Size (um)"])/1000.0*float(cfg_system["CameraXYRatio"]))
    )

    cfg_aq = cfg["Acquisition"]
    # Number of pixel in the detector
    geometry.nDetector = numpy.array((float(cfg_aq["Number of Rows"]), float(cfg_aq["Number of Columns"])))
    
    # Total size of the detector
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    ## Offset of the detector:
    geometry.offDetector = numpy.array(
        (0.0, 0.0)
    )

    # Size of each voxel
    geometry.dVoxel = numpy.array(
        (float(cfg_aq["Image Pixel Size (um)"])/1000, float(cfg_aq["Image Pixel Size (um)"])/1000, float(cfg_aq["Image Pixel Size (um)"])/1000)
    )
    geometry.nVoxel=numpy.array((geometry.nDetector[0], geometry.nDetector[1], geometry.nDetector[1]))
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel

    #% Global geometry
    geometry.DSO = float(cfg_aq["Object to Source (mm)"])
    geometry.DSD = float(cfg_aq["Camera to Source (mm)"])


    # I dont like the image geometry bruker gives:
    mag=geometry.DSD/geometry.DSO
    geometry.dVoxel=numpy.array((geometry.dVoxel[0],geometry.dVoxel[0],geometry.dVoxel[1]))/mag
    geometry.nVoxel=numpy.array((geometry.nDetector[0], geometry.nDetector[1], geometry.nDetector[1]))
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel


    geometry.whitelevel=2**int(cfg_aq["Depth (bits)"])

    angles=numpy.arange(0.0, float(cfg_aq["Number of Files"])*float(cfg_aq["Rotation Step (deg)"]), float(cfg_aq["Rotation Step (deg)"]))
    angles=angles[:-1]*numpy.pi/180
   
    return filepath, geometry, angles

def load_Bruker_projections(folder, geometry, angles, **kwargs):

    angles, indices = parse_inputs(geometry, angles, **kwargs)

    # load images
    files = sorted([file for file in os.listdir(folder) if file.lower().endswith(".tif")])

    image = Image.open(os.path.join(folder, files[indices[0]]))
    image = numpy.asarray(image).astype(numpy.float32)
    projections = numpy.zeros([len(indices),image.shape[0],image.shape[1]],dtype=numpy.single)
    projections[0,:,:] = -numpy.log(image / float(geometry.whitelevel))
    index=1
 
    print("Loading Bruker Skyscan dataset: " + folder)
    for i in tqdm(indices[1:]):
        image = Image.open(os.path.join(folder, files[i]))
        image = numpy.asarray(image).astype(numpy.float32)
        projections[index,:,:]=(-numpy.log(image / float(geometry.whitelevel)))
        index=index+1
    del geometry.whitelevel

    return numpy.asarray(projections), geometry, angles


## This should be on a separate "io_common.py" file.
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