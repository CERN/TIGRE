# For simplicity I divide the A into M=2 N=2 
# I selected 1 row blcoks and 4 columen blocks each time
import tigre
#import copy
import numpy as np
from tigre.demos.Test_data import data_loader
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
rank=0
geo = tigre.geometry(mode='cone', nVoxel=np.array([512,512,512]),default=True)
geo.dDetector = np.array([0.8, 0.8])*2               # size of each pixel            (mm)
geo.sDetector = geo.dDetector * geo.nDetector
nangles = 30
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
head = data_loader.load_head_phantom(geo.nVoxel)
y = tigre.Ax(head,geo,angles,rank)
#print('it runs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')


