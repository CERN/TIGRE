from __future__ import print_function
import tigre
import numpy as np
from testforlog import testandlog
import sys
import traceback
nVoxel = np.array([256, 256, 256])


def configuration1(alg):
    try:

        geo = tigre.geometry(mode='cone', nVoxel=nVoxel, default=True)
        niter = 10
        nangles = 100
        angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
        angles_2 = np.zeros((nangles), dtype=np.float32)
        angles = np.vstack((angles_1, angles_2, angles_2)).T
        testandlog([alg], geo, angles, niter, saveresult=True, subdirectoryname='configuration1', sup_kw_warning=True)


    except Exception as e:
        print(traceback.format_exc())
        raise SystemExit()


def configuration2(alg):
    try:
        geo = tigre.geometry(mode='parallel', nVoxel=nVoxel, default=True)
        niter = 10
        nangles = 100
        angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
        angles_2 = np.zeros((nangles), dtype=np.float32)
        angles = np.vstack((angles_1, angles_2, angles_2)).T
        testandlog([alg], geo, angles, niter, saveresult=True, subdirectoryname='configuration2', sup_kw_warning=True)
    except Exception as e:
        print(traceback.format_exc())
        raise SystemExit()


def configuration3(alg):
    try:
        geo = tigre.geometry(mode='cone', nVoxel=np.array([32, 64, 128]), default=True)
        niter = 10
        nangles = 100
        angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
        angles_2 = np.zeros((nangles), dtype=np.float32)
        angles = np.vstack((angles_1, angles_2, angles_2)).T
        testandlog([alg], geo, angles, niter, saveresult=True, subdirectoryname='configuration3', sup_kw_warning=True)
    except Exception as e:
        print(traceback.format_exc())
        raise SystemExit()


def configuration4(alg):
    try:
        geo = tigre.geometry(mode='parallel', nVoxel=np.array([32, 64, 128]), default=True)
        niter = 10
        nangles = 100
        angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
        angles_2 = np.zeros((nangles), dtype=np.float32)
        angles = np.vstack((angles_1, angles_2, angles_2)).T
        testandlog([alg], geo, angles, niter, saveresult=True, subdirectoryname='configuration4', sup_kw_warning=True)
    except Exception as e:
        print(traceback.format_exc())
        raise SystemExit()


def configuration5(alg):
    try:
        geo = tigre.geometry(mode='cone', nVoxel=np.array([32, 64, 128]), default=True)
        niter = 10
        nangles = 100
        angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
        angles_2 = np.ones((nangles), dtype=np.float32) * np.pi
        angles = np.vstack((angles_1, angles_2, angles_2)).T
        testandlog([alg], geo, angles, niter, saveresult=True, subdirectoryname='configuration5', sup_kw_warning=True)
    except Exception as e:
        print(traceback.format_exc())
        raise SystemExit()


if __name__ == '__main__':
    method_name = sys.argv[1]
    parameter_name = sys.argv[2]
    getattr(sys.modules[__name__], method_name)(parameter_name)
