import tigre
import tigre.algorithms as algs
import numpy as np

geo = tigre.geometry(mode = 'cone', nVoxel=np.array([32,64,128]), default_geo=True)
from tigre.demos.Test_data import data_loader

img = data_loader.load_head_phantom(geo.nVoxel)
angles = np.linspace(0,np.pi*2,100,dtype=np.float32)
proj = tigre.Ax(img,geo,angles)

algs.fdk(proj,geo,angles)

