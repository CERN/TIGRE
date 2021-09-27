from __future__ import print_function
from __future__ import with_statement
from __future__ import division

import os
import math
import numpy
import scipy.signal

from PIL import Image
from tigre.utilities.geometry import Geometry
import xml.etree.ElementTree as ET

def LoadScan(directory, binning=1, step=1, count=None, average=True, histogram=True, verbose=False):

	geometry, angles = ParseGeometry(directory, verbose)

	if step > 1: angles = angles[::step]
	if count is not None: angles = angles[:count]

	projections = LoadProjections(os.path.join(directory, 'Projections'), binning=binning, count=count, step=step, average=average, histogram=histogram, verbose=verbose)

	if binning > 1:
		geometry.dDetector *= binning
		geometry.nDetector = numpy.array((projections.shape[1], projections.shape[2]))

	geometry.sDetector = geometry.nDetector * geometry.dDetector
	geometry.nVoxel = numpy.array((geometry.nDetector[0], geometry.nDetector[1], geometry.nDetector[1]))
	geometry.sVoxel = numpy.array((geometry.sDetector[0], geometry.sDetector[1], geometry.sDetector[1])) * geometry.DSO / geometry.DSD
	geometry.dVoxel = geometry.sVoxel / geometry.nVoxel

	if verbose: print(geometry, angles, projections.shape, sep='\n\n', end='\n\n')

	return geometry, angles, projections

def ParseGeometry(directory, verbose=False):

	geometry = Geometry()

	# cone beam scan parameters
	file = os.path.join(directory, 'ScanParameter.xml')
	ns = dict([node for _, node in ET.iterparse(file, events=['start-ns'])])
	tree = ET.parse(file)

	params = tree.getroot().find('ConeBeamScanParameter')
	if params is None:
		raise RuntimeError('Cone beam scan parameters not found')

	# projection information

	rotDetector = params.find('DetectorOrientation')
	if rotDetector is not None:
		# TODO: handle DetectorOrientation
		if rotDetector.text.lower() != 'original':
			print('Detector orientation found but ignored!')

	flipDetector = params.find('DetectorFlipDirections')
	if flipDetector is not None:
		# TODO: handle DetectorFlipDirections
		if flipDetector.text.lower() != 'none':
			print('Detector flip direction found but ignored!')

	pitch = params.find('PixelSize')
	geometry.dDetector = numpy.array([float(pitch.find('Height').text), float(pitch.find('Width').text)])

	size = params.find('ProjectionSize')
	geometry.nDetector = numpy.array([int(size.find('Height').text), int(size.find('Width').text)])

	# trajectory information

	trajectory = params.find('ScanTrajectory')
	if trajectory is None:
		raise RuntimeError('Unable to identify scan trajectory')

	scan_type = trajectory.attrib['{{{}}}type'.format(ns['xsi'])]
	if scan_type != 'CircleTrajectory':
		raise ValueError('Unsupported trajectory type: ' + scan_type)

	geometry.mode = 'cone'
	geometry.filter = None
	geometry.accuracy = 0.5

	scale = float(trajectory.find('PixelPitch').text)
	geometry.DSD = scale * float(trajectory.find('Fdd').text)
	geometry.DSO = scale * float(trajectory.find('Fod').text)

	# TODO: check meaning of offset
	offsY = -scale * float(trajectory.find('YoCenter').text)
	if offsY != 0:
		print("Object Y-axis center found but ignored!")

	geometry.offOrigin = numpy.array((0, 0, 0))

	count_proj = int(trajectory.find('NumberOfProjections').text)
	start_angle = float(trajectory.find('StartAngle').text)
	range_angle = float(trajectory.find('Rotation').text)

	angles = math.radians(start_angle) - numpy.arange(count_proj) * math.radians(range_angle / count_proj)

	dejustment = params.find('DejustmentParameter')
	if dejustment is not None:

		offsX = float(dejustment.find('HorizontalDetectorOffset').text)
		offsY = float(dejustment.find('VerticalDetectorOffset').text)
		geometry.offDetector = -scale * numpy.array((offsY, offsX))

		# TODO: check meaning of parameters
		tiltA = math.radians(float(dejustment.find('DetectorTiltA').text))
		tiltB = math.radians(float(dejustment.find('DetectorTiltB').text))
		tiltC = math.radians(float(dejustment.find('DetectorTiltC').text))
		if tiltA != 0 or tiltB != 0 or tiltC != 0:
			print("Detector tilt found but ignored!")

		geometry.rotDetector = numpy.array((0, 0, 0))

	else:
		geometry.offDetector = numpy.array((0, 0))
		geometry.rotDetector = numpy.array((0, 0, 0))

	return geometry, angles

def BinProjection(array, factor, average=True):
	if average:
		if array.shape[0] % factor or array.shape[1] % factor:
			array = array[:(array.shape[0] // factor) * factor, :(array.shape[1] // factor) * factor]
		shape = (array.shape[0] // factor, factor, array.shape[1] // factor, factor)
		return array.reshape(shape).mean(axis=(-1, 1), dtype=numpy.float32)
	else:
		return array[::factor, ::factor]

def LoadProjections(directory, binning=1, step=1, count=None, average=True, histogram=True, verbose=False):

	projections = []

	files = sorted([file for file in os.listdir(directory) if file.lower().startswith('projection_') and file.lower().endswith(('.tif', '.tiff'))])

	if step > 1: files = files[::step]
	if count is not None: files = files[:count]

	for file in files:

		image = Image.open(os.path.join(directory, file))
		if verbose: print(file, image.format, image.size, image.mode)

		array = numpy.asarray(image).astype(numpy.float32)
		if binning > 1: array = BinProjection(array, binning, average)

		# normalization & logarithmization
		if histogram:
			hist = numpy.histogram(array, bins='scott')
			peaks = scipy.signal.find_peaks_cwt(hist[0], numpy.arange(1, len(hist[0]) // 10))
			peak = peaks[-1]
			I0 = hist[1][peak]
		else:
			I0 = numpy.max(array)

		array = -numpy.log(array / I0)

		projections.append(array)

	return numpy.asarray(projections)
