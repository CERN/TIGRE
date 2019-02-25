import tigre
import numpy as np
from testforlog import testandlog
try:
    nVoxel =np.array([512,512,512])
    geo = tigre.geometry(mode='cone', nVoxel = nVoxel, default_geo=True)
    niter = 10
    nangles = 100
    angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
    angles_2 = np.zeros((nangles), dtype=np.float32)
    angles = np.vstack((angles_1, angles_2, angles_2)).T
    testandlog(['cgls'], geo, angles,niter,saveresult=True)
except Exception as e:
    print(e)
    raise SystemExit()
raise SystemExit(__file__ + ' completed successfully')
