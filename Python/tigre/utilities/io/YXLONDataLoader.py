from __future__ import print_function
from __future__ import with_statement
from __future__ import division

import os
import math
from re import L
import numpy
import scipy.signal
from tqdm import tqdm

from PIL import Image
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET
import csv


def YXLONDataLoader(filepath, **kwargs):
    # YXLONDataLoader(filepath) Loads Nikon uCT datasets into YXLON standard
    #
    # YXLONDataLoader(filepath, OPT=VAL, ...) uses options and values.
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

    folder, geometry, angles = readYXLONGeometry(filepath)
    return loadYXLONProjections(folder, geometry, angles, **kwargs)


def readYXLONGeometry(filepath):

    # Developed by A. Biguri and W. Sun
    # W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)
    # Modified for YXLON readers by A. Biguri

    # There are 2 types of YXLON datasets that I've seen. 
    # the ones that have ScannerParameters.xml and the ones that have a txt with the filename as scanner parameters.

    if os.path.isfile(os.path.join(filepath, "ScanParameter.xml")):
        folder, geometry, angles = __YXLON_xml_geometry(filepath)
    else:
        files = [file for file in os.listdir(filepath) if file.endswith(".txt")]
        if not files:
            raise ValueError("No Geometry file found in folder: " + folder)

        folder, geometry, angles = __YXLON_txt_geometry(filepath)
    return  folder, geometry, angles

def __YXLON_xml_geometry(filepath):

    geometry = Geometry()

    # cone beam scan parameters
    file = os.path.join(filepath, "ScanParameter.xml")
    ns = dict([node for _, node in ET.iterparse(file, events=["start-ns"])])
    tree = ET.parse(file)

    params = tree.getroot().find("ConeBeamScanParameter")
    if params is None:
        raise RuntimeError("Cone beam scan parameters not found")

    # projection information

    rotDetector = params.find("DetectorOrientation")
    if rotDetector is not None:
        # TODO: handle DetectorOrientation
        if rotDetector.text.lower() != "original":
            print("Detector orientation found but ignored!")

    flipDetector = params.find("DetectorFlipDirections")
    if flipDetector is not None:
        # TODO: handle DetectorFlipDirections
        if flipDetector.text.lower() != "none":
            print("Detector flip direction found but ignored!")

    pitch = params.find("PixelSize")
    geometry.dDetector = numpy.array(
        [float(pitch.find("Height").text), float(pitch.find("Width").text)]
    )

    size = params.find("ProjectionSize")
    geometry.nDetector = numpy.array([int(size.find("Height").text), int(size.find("Width").text)])
    geometry.sDetector = geometry.nDetector * geometry.dDetector

    # trajectory information

    trajectory = params.find("ScanTrajectory")
    if trajectory is None:
        raise RuntimeError("Unable to identify scan trajectory")

    scan_type = trajectory.attrib["{{{}}}type".format(ns["xsi"])]
    if scan_type != "CircleTrajectory":
        raise ValueError("Unsupported trajectory type: " + scan_type)

    geometry.mode = "cone"
    geometry.filter = None
    geometry.accuracy = 0.5

    scale = float(trajectory.find("PixelPitch").text)
    geometry.DSD = scale * float(trajectory.find("Fdd").text)
    geometry.DSO = scale * float(trajectory.find("Fod").text)

    # TODO: check meaning of offset
    geometry.COR = -scale * float(trajectory.find("YoCenter").text)

    geometry.offOrigin = numpy.array((0, 0, 0))

    count_proj = int(trajectory.find("NumberOfProjections").text)
    start_angle = float(trajectory.find("StartAngle").text)
    range_angle = float(trajectory.find("Rotation").text)

    angles = math.radians(start_angle) - numpy.arange(count_proj) * math.radians(
        range_angle / count_proj
    )

    dejustment = params.find("DejustmentParameter")
    if dejustment is not None:

        offsX = float(dejustment.find("HorizontalDetectorOffset").text)
        offsY = float(dejustment.find("VerticalDetectorOffset").text)
        geometry.offDetector = -scale * numpy.array((offsY, offsX))

        # TODO: check meaning of parameters
        tiltA = math.radians(float(dejustment.find("DetectorTiltA").text))
        tiltB = math.radians(float(dejustment.find("DetectorTiltB").text))
        tiltC = math.radians(float(dejustment.find("DetectorTiltC").text))
        if tiltA != 0 or tiltB != 0 or tiltC != 0:
            print("Detector tilt found but ignored! Parameter matching unknown, if you have more info email tigre.toolbox@gmail.com")

        geometry.rotDetector = numpy.array((0, 0, 0))

    else:
        geometry.offDetector = numpy.array((0, 0))
        geometry.rotDetector = numpy.array((0, 0, 0))

    geometry.nVoxel = numpy.array(
        (geometry.nDetector[0], geometry.nDetector[1], geometry.nDetector[1])
    )
    geometry.sVoxel = (
        numpy.array((geometry.sDetector[0], geometry.sDetector[1], geometry.sDetector[1]))
        * geometry.DSO
        / geometry.DSD
    )
    geometry.dVoxel = geometry.sVoxel / geometry.nVoxel
    geometry.whitelevel = float(2**params.find("ProjectionBitdepth")-1)

    return filepath, geometry, angles


def __YXLON_txt_geometry(filepath):

    files = [file for file in os.listdir(filepath) if file.endswith(".txt")]
    ini = files[0]

    geometry = Geometry()
    geometry.accuracy = 0.5
    binning=1
    with open(os.path.join(filepath, ini)) as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"): #You can also use delimiter="\t" rather than giving a dialect.
            if not line:
                continue
            if line[-1]=='Number of Detector Elements ':
                geometry.nDetector=numpy.array((int(line[0]),int(line[0])))
            if line[-1]=='2D-Pixel Size [mm]':
                geometry.dDetector=numpy.array((float(line[0]),float(line[0])))
            if line[-1]=='Pixelbinning':
                binning=float(line[0])
            if line[-1]=='3D-XY-Pixel Size [mm]':
                geometry.dVoxel[1]=float(line[0])
                geometry.dVoxel[2]=float(line[0])
            if line[-1]=='3D-Z-Pixel Size [mm]':
                geometry.dVoxel[0]=float(line[0])
            if line[-1]=='Image Dimension':
                geometry.nVoxel[1]=int(line[0])
                geometry.nVoxel[2]=int(line[0])
            if line[-1]=='Number of Z-slices':
                geometry.nVoxel[0]=float(line[0])
            if line[-1]=='FOD [mm]':
                geometry.DSO=float(line[0])
            if line[-1]=='FDD [mm]':
                geometry.DSD=float(line[0])
            if line[-1]=='Center Offset [mm]':
                geometry.COR=float(line[0])
            if line[-1]=='CT: Number of Projections':
                n_angles = int(line[0])

    geometry.dDetector=geometry.dDetector * binning
    geometry.sDetector = geometry.nDetector * geometry.dDetector
    geometry.sVoxel = geometry.nVoxel * geometry.dVoxel


    ## whitelevel
    geometry.whitelevel = float(2**16-1)

    ## angles
    angles = -numpy.linspace(0,2*numpy.pi,n_angles)
    return filepath, geometry, angles


def loadYXLONProjections(folder, geometry, angles, **kwargs):
    # projections, geometry, angles = loadYXLONProjections(filepath, geometry, angles, **kwargs)
    #    loads YXLON uCT machine projections
    #
    #    loadYXLONProjections(filepath, geometry, angles) Loads a dataset given its FILEPATH,
    #       GEOMETRY and ANGLES (loaded from readYXLONGeometry())
    #       See YXLONDataLoader() for additional options.
    #
    #% developed by A. Biguri and W. Sun 06.07.2020
    # Modified by A. biguri for YXLON datasets

    # parse inputs
    angles, indices = parse_inputs(geometry, angles, **kwargs)

    # load images
    if os.path.isdir(os.path.join(folder,'Projections')):
        folder= os.path.join(folder,'Projections')
    elif os.path.isdir(os.path.join(folder,'Projektionen')):
        folder= os.path.join(folder,'Projektionen')
    else:
        raise ValueError("Folder with projections not found")

    files = sorted([file for file in os.listdir(folder) if file.lower().endswith(".raw")])
    istiff=False
    if not files:
        files = sorted([file for file in os.listdir(folder) if file.lower().endswith(".tif") or file.lower().endswith(".tiff")])
        if not files:
            raise ValueError("No projection founds with extension .raw or .tif")
        istiff=True
    
    if not istiff:
        image =  numpy.fromfile(os.path.join(folder,  files[indices[0]]), dtype=numpy.uint16)
        image = numpy.asarray(image).astype(numpy.float32)
        image = numpy.flipud(numpy.reshape(image,(geometry.nDetector)))
    else:
        image = Image.open(os.path.join(folder, files[indices[0]]))
        image = numpy.asarray(image).astype(numpy.float32)

    projections = numpy.zeros([len(indices), image.shape[0], image.shape[1]], dtype=numpy.single)
    projections[0, :, :] = -numpy.log(image / float(geometry.whitelevel))
    # Python uses zero-based indexing
    # if we start with index = 1 we get an error
    # use enumerate, it's cleaner
    print("Loading YXLON dataset: " + folder)
    for index, i in enumerate(tqdm(indices[1:]), 1):

        if not istiff:
            image =  numpy.fromfile(os.path.join(folder,  files[i]), dtype=numpy.uint16)
            image = numpy.asarray(image).astype(numpy.float32)
            image = numpy.flipud(numpy.reshape(image,(geometry.nDetector)))
        else:
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


