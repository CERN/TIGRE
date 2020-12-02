import time
import numpy as np
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb

def mlem(proj,geo,angles,niter,verbose=True):
    name="mlem"
#MLEM solves the tomographic problem by using Maximum Likelihood Expection
# Maximitation algorithm. 
#
#   MLEM(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
#   using the projection data PROJ taken over ALPHA angles, corresponding
#   to the geometry descrived in GEO, using NITER iterations.
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
# 
# Copyright (c) 2015, University of Bath and 
#                     CERN-European Organization for Nuclear Research
#                     All rights reserved.
#
# License:            Open Source under BSD. 
#                     See the full license at
#                     https://github.com/CERN/TIGRE/blob/master/LICENSE
#
# Contact:            tigre.toolbox@gmail.com
# Codes:              https://github.com/CERN/TIGRE/
# Coded by:           Ander Biguri 
#--------------------------------------------------------------------------
    res=np.ones(geo.nVoxel, np.float32)
    W=Atb(np.ones((angles.shape[0], geo.nDetector[0], geo.nDetector[1]), np.float32), geo, angles)

    for ii in range(niter):
        if verbose:
            if ii == 0:
                print(str(name).upper() +
                        ' ' + "algorithm in progress.")
                toc = time.clock()
            if ii == 1:
                tic = time.clock()
                print('Esitmated time until completetion (s): ' +
                        str((niter - 1) * (tic - toc)))
        auxMLEM=proj/Ax(res,geo,angles)
        auxMLEM[np.isnan(auxMLEM)] = 0
        auxMLEM[np.isinf(auxMLEM)] = 0
        
        imgupdate = Atb(auxMLEM, geo,angles)/W
        imgupdate[np.isnan(imgupdate)] = 0
        imgupdate[np.isinf(imgupdate)] = 0
        
        res *= imgupdate
        res[res<0] = 0

    return res
