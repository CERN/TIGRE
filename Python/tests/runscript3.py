import tigre
import tigre.algorithms as algs
import numpy as np
from testforlog import testandlog
try:
    nVoxel=np.array([512, 512, 512])
    geo = tigre.geometry(mode='parallel',nVoxel=nVoxel, default_geo=True)
    niter = 10
    nangles = 100
    alglist = algs.__all__
    alglist.remove('sart')
    angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
    angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
    angles_3 = np.ones((nangles), dtype=np.float32) * 0.0001
    angles = np.vstack((angles_1, angles_2, angles_2)).T
    testandlog([alglist[1]],geo,angles,niter, saveresult=True,blocksize=20)
except Exception as e:
    print(e)
    raise SystemExit()
raise SystemExit(__file__ + ' completed successfully')