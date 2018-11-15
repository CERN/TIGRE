import numpy as np
class Recon(object):

    def __init__(self,proj,geo,angles,niter,**kwargs):
        self.proj = proj
        self.geo = geo
        self.angles = angles
        self.niter = niter
        self.__dict__.update(**kwargs)

    def backprojection(self):
        res = np.zeros(self.proj.shape)
        return res

