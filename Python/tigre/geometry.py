from __future__ import division
import numpy as np


class geometry:

    def __init__(self):

        self.mode = None

    def check_geo(self):

        if self.mode =='cone':
            manditory_attribs = ['nVoxel','sVoxel','dVoxel'
                                'nDetector', 'sDetector',
                                'DSO','DSD']




            optional_fields = ['offOrigin','offDetector','rotDetector','COR',
                               'mode','accuracy']
        if self.mode == 'paralell':
            print("doing something")
        elif self.mode == None:
            raise AttributeError("Mode needs to be specified in order for tigre to assess whether the geometry is "
                                 "correctly implemented")




    def __str__(self):
        parameters = []
        parameters.append("TIGRE parameters")
        parameters.append("-----")
        parameters.append("Distance from source to detector = " + str(self.DSD) + "mm")
        parameters.append("Distance from source to origin = " + str(self.DSO) + "mm")

        parameters.append("-----")
        parameters.append("Detector parameters")
        parameters.append("Number of pixels = " + str(self.nDetector))
        parameters.append("Size of each pixel = " + str(self.dDetector) + "mm")
        parameters.append("Total size of the detector = " + str(self.sDetector) + "mm")

        parameters.append("-----")
        parameters.append("Image parameters")
        parameters.append("Number of voxels = " + str(self.nVoxel))
        parameters.append("Total size of the image = " + str(self.sVoxel) + "mm")
        parameters.append("Size of each voxel = " + str(self.dVoxel) + "mm")

        parameters.append("-----")
        parameters.append("Offset correction parameters")
        parameters.append("Offset of image from origin = " + str(self.offOrigin) + "mm")
        parameters.append("Offset of detector = " + str(self.offDetector) + "mm")

        parameters.append("-----")
        parameters.append("Auxillary parameters")
        parameters.append("Accuracy of forward projection = " + str(self.accuracy))

        return '\n'.join(parameters)
