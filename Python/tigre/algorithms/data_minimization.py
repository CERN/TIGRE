from __future__ import division
import numpy as np
from tigre.Atb import Atb
from tigre.Ax import Ax
import copy
class DataMinimization(object):
    def art_data_minimizing(self):
        geo = copy.deepcopy(self.geo)
        for j in range(len(self.angleblocks)):
            if self.blocksize == 1:
                angle = np.array([self.angleblocks[j]], dtype=np.float32)
            else:
                angle = self.angleblocks[j]
            if geo.offOrigin.shape[0] ==self.angles.shape[0]:
               geo.offOrigin = self.geo.offOrigin[j]
            if geo.offDetector.shape[0] == self.angles.shape[0]:
                geo.offOrin = self.geo.offDetector[j]
            if geo.rotDetector.shape[0] ==self.angles.shape[0]:
                geo.rotDetector=self.geo.rotDetector[j]
            if hasattr(geo.DSD,'shape'):
                if geo.DSD.shape[0] ==self.angles.shape[0]:
                    geo.DSD = self.geo.DSD[j]
            if hasattr(geo.DSO,'shape'):
                if geo.DSO.shape[0] ==self.angles.shape[0]:
                    geo.DSO = self.geo.DSO[j]
            self.res += self.lmbda * 1/self.third_dim_sum(self.V[:,:,self.angle_index[j]]) * Atb(self.W[self.angle_index[j]] * (self.proj[self.angle_index[j]]
                                     - Ax(self.res, geo, angle, 'interpolated')),geo, angle, 'FDK')

    def third_dim_sum(self,V):
        if V.ndim == 3:
            return np.sum(V, axis=2, dtype=np.float32)
        else:
            return V












    