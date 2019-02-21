from __future__ import print_function
import tigre
import sys, traceback
import os
import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
import time
import tigre.algorithms as algs
import warnings
warnings.filterwarnings("error")

def testandlog(alglist, geo, angles, niter, logfilename=None,saveresult=False):
    nangles = angles.shape[0]
    # createlogfile

    dirname = os.path.dirname(__file__)
    logfolder = os.path.join(dirname, 'logTIGRE')
    if 'logTIGRE' not in os.listdir(dirname):
        os.system('mkdir ' + logfolder)

    if logfilename is not None:
        logf = open(os.path.join(logfolder, logfilename) + '.log', 'w')
    else:
        logf = open(os.path.join(logfolder, str(geo.mode) + '_' + str(time.asctime()).replace(' ', '_') + '.log'), 'w')
    source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
    proj = tigre.Ax(source_img, geo, angles)
    logf.write('GEOMETRY used for instance of testandlog: \n')
    for item in geo.__dict__:
        logf.write(item + ': ' + str(getattr(geo, item)) + '\n')
    logf.write('------------------------------------------------\n')
    for alg in alglist:
        logf.write(str(alg).upper() + ' ' + time.asctime() + '\n')
        logf.write('nproj: ' + str(nangles) + ' niter: ' + str(niter) + '\n')
        if np.sum(angles[1])!=0:
            spherical = True
        else:
            spherical = False
        logf.write('spherical projection: ' + str(spherical) + '\n')
        algpassed = False
        try:
            tic = time.clock()
            res = getattr(algs, alg)(proj, geo, angles, niter=niter)
            algpassed = True
        except Exception:
            formatedlines = traceback.format_exc()
            logf.write(''.join(formatedlines) + '\n')
        finally:
            toc = time.clock()
            pass
        if algpassed:
            logf.write('total time taken: ' + str(abs(tic - toc)) + '\n')
        logf.write('------------------------------------------------\n')

        if saveresult and algpassed:
            npsavefile = os.path.join(logfolder, 'numpysavefile')
            if 'numpysavefile' not in os.listdir(logfolder):
                os.system('mkdir '+npsavefile)
            np.save(os.path.join(npsavefile,alg+'_'+os.path.splitext(os.path.split(logf.name)[-1])[0]),res)



alglist = ['sirt','cgls','fbp']
#alglist.remove('sart')

geo = tigre.geometry(mode='cone', nVoxel=np.array([32, 64, 128]), default_geo=True)
niter = 10
nangles = 100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.ones((nangles), dtype=np.float32) * 0.0001
angles = np.vstack((angles_1, angles_3, angles_3)).T


testandlog(alglist, geo, angles,niter,saveresult=False)

geo = tigre.geometry(mode='parallel', nVoxel=np.array([32, 64, 128]), default_geo=True)
niter = 10
nangles = 100

testandlog(alglist, geo, angles,niter,saveresult=True)