import numpy as np
from tigre.Ax import Ax
from tigre.Atb import Atb
from data_minimization import DataMinimization
from regularisation import Regularisation
from tigre.Utilities.order_subsets import order_subsets
from tigre.Utilities.init_multigrid import init_multigrid
from tigre.Algorithms.FDK import FDK
import copy

class IterativeReconAlg(Regularisation, DataMinimization):

    def __init__(self,proj,geo,angles,niter,**kwargs):

        self.proj = proj
        self.angles = angles
        self.geo = geo
        self.niter = niter

        options = dict(blocksize=20, lmbda=1, lmbda_red=0.99, OrderStrategy=None, Quameasopts=None, init=None,
                       verbose=True, noneg=True, computel2=False, dataminimizing = 'default_data_minimizing')
        for kw in kwargs.keys()
            if not options.has_key(kw):
                #Note: might not want this warning (typo checking).
                print("Warning: " + kw + "not recognised as default parameter for instance of IterativeReconAlg.")
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)
        if not hasattr(self,'W'):
            self.set_w()
        if not hasattr(self,'V'):
            self.set_v()
        if not hasattr(self,'res'):
            self.set_init_parameters()
    def set_w(self):

        geox = copy.deepcopy(self.geo)
        geox.sVoxel[0:] = self.geo.DSD - self.geo.DSO
        geox.sVoxel[2] = max(geox.sDetector[1], geox.sVoxel[2])
        geox.nVoxel = np.array([2, 2, 2])
        geox.dVoxel = geox.sVoxel / geox.nVoxel
        W = Ax(np.ones(geox.nVoxel, dtype=np.float32), geox, self.angles, "ray-voxel")
        W[W < min(self.geo.dVoxel / 4)] = np.inf
        W = 1 / W
        setattr(self,'W',W)

    def set_v(self):
        geo = self.geo
        if geo.mode != 'parallel':

            start = geo.sVoxel[1] / 2 - geo.dVoxel[1] / 2 + geo.offOrigin[1]
            stop = -geo.sVoxel[1] / 2 + geo.dVoxel[1] / 2 + geo.offOrigin[1]
            step = -geo.dVoxel[1]
            xv = np.arange(start, stop + step, step)

            start = geo.sVoxel[2] / 2 - geo.dVoxel[2] / 2 + geo.offOrigin[2]
            stop = -geo.sVoxel[2] / 2 + geo.dVoxel[2] / 2 + geo.offOrigin[2]
            step = -geo.dVoxel[2]
            yv = -1 * np.arange(start, stop + step, step)

            (yy, xx) = np.meshgrid(yv, xv)
            xx = np.expand_dims(xx, axis=2)
            yy = np.expand_dims(yy, axis=2)

            A = (self.angles + np.pi / 2)
            V = (geo.DSO / (geo.DSO + (yy * np.sin(-A)) - (xx * np.cos(-A)))) ** 2
            V = np.array(V, dtype=np.float32)


        else:
            V = np.ones([self.angles.shape[0], geo.nVoxel[1], geo.nVoxel[0]], dtype=np.float32)

        setattr(self,'V',V)

    def set_init_parameters(self):

        setattr(self,'res', np.zeros(self.geo.nVoxel, dtype=np.float32))
        init = self.init
        verbose = self.verbose
        if init == 'multigrid':
            if verbose:
                print('init multigrid in progress...')
                print('default blocksize=1 for init_multigrid(OS_SART)')
            self.res = init_multigrid(self.proj, self.geo, self.angles, alg='SART')
            if verbose:
                print('init multigrid complete.')
        if init == 'FDK':
            self.res = FDK(self.proj, self.geo, self.angles).transpose()

        if type(init) == np.ndarray:
            if (self.geo.nVoxel == init.shape).all():

                self.res = init

            else:
                raise ValueError('wrong dimension of array for initialisation')


    def main_iter(self):







if __name__=='__main__':
    new_class = IterativeReconAlg('proj', 'geo', 'angles', 'niter', **{'e': 1, 'f': 2, 'g': 4})
    new_class.change_proj()
    print(new_class.__dict__)
