import numpy as np
from tigre.Atb import Atb
from tigre.Ax import Ax

class DataMinimization(object):
    def default_data_minimizing(self):
        for j in range(len(self.angleblocks)):
            if self.blocksize == 1:
                angle = np.array([self.angleblocks[j]], dtype=np.float32)
            else:
                angle = self.angleblocks[j]

            self.res += self.lmbda * 1/self.third_dim_sum(self.V[:,:,self.angle_index[j]]) * Atb(self.W[self.angle_index[j]] * (self.proj[self.angle_index[j]]
                                     - Ax(self.res, self.geo, angle, 'interpolated')),self.geo, angle, 'FDK')

    def third_dim_sum(self,V):
        if V.ndim == 3:
            return np.sum(V, axis=2, dtype=np.float32)
        else:
            return V












    