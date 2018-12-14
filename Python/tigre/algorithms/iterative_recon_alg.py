import numpy as np
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
from tigre.utilities.order_subsets import order_subsets
from tigre.utilities.init_multigrid import init_multigrid
from tigre.utilities.Measure_Quality import Measure_Quality as MQ
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.algorithms.single_pass_algorithms import FDK
from _minTV import minTV
from _AwminTV import AwminTV
import time
import copy
"""
This module is where the umbrella class IterativeReconAlg is located
which is the umbrella class to all the other algorithms apart from 
the single pass type algorithms. 
"""

# coding: utf8
class DataMinimization(object):
    """
    Class used to define the methods used in run_main_iter in IterativeReconAlg.
    """
    def art_data_minimizing(self):
        """
        VERBOSE:
        for j in range(angleblocks):
            angle = np.array([alpha[j]], dtype=np.float32)
            proj_err = proj[angle_index[j]] - Ax(res, geo, angle, 'ray-voxel')
            weighted_err = W[angle_index[j]] * proj_err
            backprj = Atb(weighted_err, geo, angle, 'FDK')
            weighted_backprj = 1 / V[angle_index[j]] * backprj
            res += weighted_backprj
            res[res<0]=0

        :return: None
        """
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
            if self.noneg:
                self.res = self.res.clip(min=0)

    def third_dim_sum(self,V):
        if V.ndim == 3:
            return np.sum(V, axis=2, dtype=np.float32)
        else:
            return V

class Regularisation(object):

    def minimizeTV(self,res_prev,dtvg):
        return minTV(res_prev,dtvg,self.numiter_tv)

    def minimizeAwTV(self,res_prev,dtvg):
        return AwminTV(res_prev,dtvg,self.numiter_tv,self.delta)


class IterativeReconAlg(Regularisation, DataMinimization):
    """
    Parameters
    ----------
    :param proj: (np.ndarray, dtype=np.float32)
    Input data, shape = (geo.nDector, nangles)

    :param geo: (tigre.geometry)
    Geometry of detector and image (see examples/Demo code)

    :param angles: (np.ndarray , dtype=np.float32)
    angles of projection, shape = (nangles,3)

    :param niter: (int)
    number of iterations for reconstruction algorithm

    :param kwargs: (dict)
    optional parameters

    Keyword Arguments
    -----------------
    :keyword blocksize: (int)
        number of angles to be included in each iteration
        of proj and backproj for OS_SART
    :keyword lmbda: (np.float64)
        Sets the value of the hyperparameter.

    :keyword lmbda_red: (np.float64)
        Reduction of lambda every iteration
        lambda=lambdared*lambda. Default is 0.99

    :keyword init: (str)
        Describes different initialization techniques.
              "none"     : Initializes the image to zeros (default)
              "FDK"      : intializes image to FDK reconstrucition
              "multigrid": Initializes image by solving the problem in
                           small scale and increasing it when relative
                           convergence is reached.
              "image"    : Initialization using a user specified
                           image. Not recommended unless you really
                           know what you are doing.

    :keyword InitImg: (np.ndarray)
        Not yet implemented. Image for the "image" initialization.

    :keyword verbose:  (Boolean)
        Feedback print statements for algorithm progress
        default=True

    :keyword Quameasopts: (list)
        Asks the algorithm for a set of quality measurement
        parameters. Input should contain a list or tuple of strings of
        quality measurement names. Examples:
            RMSE, CC, UQI, MSSIM

    :keyword OrderStrategy : (str)
        Chooses the subset ordering strategy. Options are:
                 "ordered"        : uses them in the input order, but divided
                 "random"         : orders them randomply
                 "angularDistance": chooses the next subset with the
                                    biggest angular distance with the ones used
    Examples
    --------
    tigre.demos.run() to launch ipython notebook file with examples.

     """

    def __init__(self, proj, geo, angles, niter, **kwargs):

        self.proj = proj
        self.angles = angles
        self.geo = geo
        self.niter = niter

        options = dict(blocksize=20, lmbda=1, lmbda_red=0.99,
                       OrderStrategy=None, Quameasopts=None,
                       init=None,verbose=True, noneg=True,
                       computel2=False, dataminimizing='art_data_minimizing',
                       regularisation='minimizeTV')
        allowed_keywords = ['V','W','log_parameters','angleblocks','angle_index','delta']
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)
        for kw in kwargs.keys():
            if not options.has_key(kw) and (kw not in allowed_keywords):
                if self.verbose:
                    # Note: might not want this warning (typo checking).
                    print("Warning: " + kw + " not recognised as default parameter for instance of IterativeReconAlg.")
        if self.angles.ndim == 1:
            a1 = self.angles
            a2 = np.zeros(self.angles.shape[0], dtype=np.float32)
            setattr(self, 'angles', np.vstack((a1, a2, a2)).T)
        if not hasattr(self, 'W'):
            self.set_w()
        if not hasattr(self, 'V'):
            self.set_v()
        if not hasattr(self, 'res'):
            self.set_res()
        if not all([hasattr(self, 'angleindex'), hasattr(self, 'angleblocks')]):
            self.set_angle_index()
        setattr(self, 'lq', [])  # quameasoptslist
        setattr(self, 'l2l', [])  # l2list

    def set_w(self):
        """
        Calculates value of W if this is not given.
        :return: None
        """
        geox = copy.deepcopy(self.geo)
        geox.sVoxel[0:] = self.geo.DSD - self.geo.DSO
        geox.sVoxel[2] = max(geox.sDetector[1], geox.sVoxel[2])
        geox.nVoxel = np.array([2, 2, 2])
        geox.dVoxel = geox.sVoxel / geox.nVoxel
        W = Ax(np.ones(geox.nVoxel, dtype=np.float32), geox, self.angles, "ray-voxel")
        W[W < min(self.geo.dVoxel / 4)] = np.inf
        W = 1 / W
        setattr(self, 'W', W)

    def set_v(self):
        """
        Computes value of V parameter if this is not given.
        :return: None
        """
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
            A = (self.angles[:, 0] + np.pi / 2)
            V = (geo.DSO / (geo.DSO + (yy * np.sin(-A)) - (xx * np.cos(-A)))) ** 2
            V = np.array(V, dtype=np.float32)
            setattr(self, 'V', V)

        else:
            V = np.ones([ geo.nVoxel[1], geo.nVoxel[2], self.angles.shape[0]], dtype=np.float32)
            setattr(self, 'V', V)

    def set_res(self):
        """
        Calulates initial value for res if this is not given.
        :return: None
        """
        setattr(self, 'res', np.zeros(self.geo.nVoxel, dtype=np.float32))
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
            self.res = FDK(self.proj, self.geo, self.angles)

        if type(init) == np.ndarray:
            if (self.geo.nVoxel == init.shape).all():

                self.res = init

            else:
                raise ValueError('wrong dimension of array for initialisation')

    def set_angle_index(self):
        """
        sets angle_index and angleblock if this is not given.
        :return: None
        """
        angleblocks, angle_index = order_subsets(self.angles, self.blocksize, self.OrderStrategy)
        setattr(self, 'angleblocks', angleblocks)
        setattr(self, 'angle_index', angle_index)

    def run_main_iter(self):
        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        Quameasopts = self.Quameasopts

        for i in range(self.niter):
            res_prev = None
            if Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            if self.verbose:
                if i == 0:
                    print("Algorithm in progress.")
                    toc = time.clock()
                if i == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            getattr(self, self.dataminimizing)()
            self.error_measurement(res_prev, i)

    def error_measurement(self, res_prev, iter):
        if self.Quameasopts is not None and iter > 0:
            self.lq.append(MQ(self.res, res_prev, self.Quameasopts))
        if self.computel2:
            # compute l2 borm for b-Ax
            errornow = im3DNORM(self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel'), 2)
            self.l2l.append(errornow)

    def getres(self):
        return self.res

    def getl2(self):
        return self.l2l


def decorator(IterativeReconAlg, name=None, docstring=None):
    """
    Calls run_main_iter when parameters are given to it.

    :param IterativeReconAlg: obj, class
        instance of IterativeReconAlg
    :param name: str
        for name of func
    :param docstring: str
        other documentation that may need to be included from external source.
    :return: func

    Examples
    --------
    import tigre
    from tigre.demos.Test_data.data_loader import load_head_phantom
    geo = tigre.geometry_defaut(high_quality=False)
    src = load_head_phantom(number_of_voxels=geo.nVoxel)
    proj = Ax(src,geo,angles)
    angles = np.linspace(0,2*np.pi,100)
    iterativereconalg = decorator(IterativeReconAlg)
    output = iterativereconalg(proj,geo,angles, niter=50)

    """

    def iterativereconalg(proj, geo, angles, niter, **kwargs):
        alg = IterativeReconAlg(proj, geo, angles, niter, **kwargs)
        alg.run_main_iter()
        if alg.computel2:
            return alg.getres(), alg.getl2()
        else:
            return alg.getres()

    if docstring is not None:
        setattr(iterativereconalg, '__doc__', docstring + IterativeReconAlg.__doc__)
    else:
        setattr(iterativereconalg, '__doc__', IterativeReconAlg.__doc__)
    if name is not None:
        setattr(iterativereconalg, '__name__', name)
    return iterativereconalg


iterativereconalg = decorator(IterativeReconAlg)
