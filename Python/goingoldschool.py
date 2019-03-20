from __future__ import division
import numpy as np
import tigre
import copy
from _tvdenoising import tvdenoise
from tigre.demos.Test_data import data_loader


def fista(y, geo, angles, niter, trueimg):
    lmbda = 0.1
    y_rec = np.zeros((geo.nVoxel), dtype=np.float32)
    x_rec = copy.deepcopy(y_rec)
    t = 1
    L = 2.e5
    bm = 1 / L
    relativeError = []
    for i in range(niter):
        y_rec += bm * 2 * (tigre.Atb((y - tigre.Ax(y_rec, geo, angles, 'interpolated')),
                                     geo, angles, 'FDK'))


        #lambdaForTv = 2 * bm * lmbda
        find_nan(y_rec,22,i)
        y_rec = tvdenoise(y_rec, 20, 12.0)
        find_nan(y_rec,22,i)
        """
        x_rec_old = copy.deepcopy(x_rec)
        x_rec = tvdenoise(y_rec, 0, lambdaForTv)
        t_old = t
        t = (1 + np.sqrt(1 + 4 * t ** 2)) / 2
        y_rec = x_rec + (t_old - 1) / t * (x_rec - x_rec_old)
        """
        relativeError.append(np.linalg.norm((y_rec-trueimg).ravel(),2)/np.linalg.norm(trueimg.ravel(),2))

    return y_rec ,relativeError
def find_nan(img,line,iteration):
    for i in img.ravel():
        if str(i) == 'nan':
            print(img)
            print('NAN DETECTED at line ' +str(line))
            print('iteration ' +str(iteration))

            raise SystemExit
angles = np.linspace(0,2*np.pi,200)


geo = tigre.geometry(mode='cone', nVoxel=np.array([64, 64, 64]), default_geo=True)

src_img = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(src_img,geo,angles)
output, relativeError = fista(proj,geo,angles,50,src_img)
from matplotlib import pyplot as plt
plt.plot(relativeError)
plt.show()
