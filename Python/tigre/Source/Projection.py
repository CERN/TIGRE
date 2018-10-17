import numpy as np
import setuptools

class Projection:

    def __init__(self,data, geo, angles, **kwargs):

        if self.checkgeo(geo):
            self.geo = geo

        if self.checkinputs(data,angles):
            self.angles = angles
            self.data = data

    def checkgeo(self,geo):
        if type(geo)==dict:





