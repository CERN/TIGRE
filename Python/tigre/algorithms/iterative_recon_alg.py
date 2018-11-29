import numpy as np
from tigre.Ax import Ax
from tigre.algorithms.data_minimization import DataMinimization
from tigre.algorithms.regularisation import Regularisation
from tigre.utilities.order_subsets import order_subsets
from tigre.utilities.init_multigrid import init_multigrid
from tigre.utilities.Measure_Quality import Measure_Quality as MQ
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.algorithms.FDK import FDK
import copy
import time
import copy


# coding: utf8

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

    :keyword QualMeas: (list)
        Asks the algorithm for a set of quality measurement
        parameters. Input should contain a list or tuple of strings of
        quality measurement names. These will be computed in each iteration.
        Example:
          ["CC","RMSE","MSSIM"]


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

        options = dict(blocksize=20, lmbda=1, lmbda_red=0.99, OrderStrategy=None, Quameasopts=None, init=None,
                       verbose=True, noneg=True, computel2=False, dataminimizing='default_data_minimizing')
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)
        for kw in kwargs.keys():
            if not options.has_key(kw):
                if self.verbose:
                    # Note: might not want this warning (typo checking).
                    print("Warning: " + kw + " not recognised as default parameter for instance of IterativeReconAlg.")
        if self.angles.ndim == 1:
            a1 = self.angles
            a2 = np.zeros(self.angles.shape[0],dtype=np.float32)
            setattr(self,'angles',np.vstack((a1,a2,a2)).T)
        if not hasattr(self, 'W'):
            self.set_w()
        if not hasattr(self, 'V'):
            self.set_v()
        if not hasattr(self, 'res'):
            self.set_res()
        if not all([hasattr(self, 'angleindex'), hasattr(self, 'angleblocks')]):
            self.set_angle_index()

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
            V = np.ones([self.angles.shape[0], geo.nVoxel[1], geo.nVoxel[0]], dtype=np.float32)

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
        computel2 = self.computel2

        for i in range(self.niter):
            if Quameasopts is not None:
                res_prev = self.res
                setattr(self, 'lq', [])
            if computel2:
                setattr(self, 'l2l', [])
            if self.verbose:
                if i == 0:
                    print("Algorithm in progress.")
                    toc = time.clock()
                if i == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            getattr(self, self.dataminimizing)()
            if Quameasopts is not None:
                self.lq.append(MQ(self.res, res_prev, Quameasopts))
            if computel2:
                # compute l2 borm for b-Ax

                errornow = im3DNORM(self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel'), 2)
                self.l2l.append(errornow)

        setattr(self, 'prev_geo_config', copy.deepcopy(self.geo.__dict__))
        setattr(self, 'prev_angles', copy.deepcopy(self.angles))
        setattr(self, 'prev_niter', copy.deepcopy(self.niter))

    def getres(self):
        """
        Checks if the main iteration has been performed
        on the given configuration. If it hasn't it calls
        run_main_iter() and returns the required result.

        :return res: (np.ndarray, dtype=np.float32)
        Result of performing the main iteration fo the algorithm.

        if :keyword Quameasopts is not none:
            :return res, lq: (np.ndarray, list)
                where lq is the result of tigre.utilities.MQ for
                every iteration

        if :keyword computel2:
            :return res, l2l: (np.ndarray, list)
                where l2l is the result of computing the l2 norm
                for every iteration.
        """
        prev_config = ['prev_angles', 'prev_niter', 'prev_geo_config']
        if all([hasattr(self, attrib) for attrib in prev_config]):
            if all([self.geo.issame(self.prev_geo_config), (self.angles == self.prev_angles).all(),
                    self.niter == self.prev_niter]):
                pass
        else:
            self.run_main_iter()
        if self.computel2:
            return self.res, self.l2l
        if self.Quameasopts is not None:
            return self.res, self.lq
        else:
            return self.res


def decorator(IterativeReconAlg, name=None, docstring=None):
    """
    :param IterativeReconAlg: obj, class
        instance of IterativeReconAlg
    :param name: str
        for name of func
    :param docstring: str
        other documentation that may need to be included from external source.
    :return: func
    """

    def iterativereconalg(proj, geo, angles, niter, **kwargs):
        alg = IterativeReconAlg(proj, geo, angles, niter, **kwargs)
        return alg.getres()

    if docstring is not None:
        setattr(iterativereconalg, '__doc__', docstring + IterativeReconAlg.__doc__)
    else:
        setattr(iterativereconalg, '__doc__', IterativeReconAlg.__doc__)
    if name is not None:
        setattr(iterativereconalg, '__name__', name)
    return iterativereconalg


iterativereconalg = decorator(IterativeReconAlg)
