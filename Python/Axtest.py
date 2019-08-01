import time
import scipy.io as sio
import tigre
import copy
import numpy as np
from tigre.demos.Test_data import data_loader
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
geo = tigre.geometry(mode='cone', nVoxel=np.array([512,512,512]),default=True)
Ingeo = copy.deepcopy(geo)
print(geo)
print('===============================')
print(Ingeo)
head = data_loader.load_head_phantom(geo.nVoxel)
if rank ==0:
  print('head.shape',head.shape)
  print(geo)
nangles = 500
geo.dDetector = np.array([0.8, 0.8])*2               # size of each pixel            (mm)
geo.sDetector = geo.dDetector * geo.nDetector
maxIte = 1
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
comm.Barrier()
for i in range(maxIte):
  comm.Barrier()
  print('rank',rank,'iteration,',i)
  Start_time = time.time()
  y = tigre.Ax(head,geo,angles,rank)  
  End_time = time.time()
  duration_time = End_time - Start_time
  print('rank,',rank,'Ax geo form time =',duration_time)   
  print('rank,',rank,'Ax geo y.sum', np.sum(y))
  Start_time = time.time()
  y2 = tigre.Ax(head,Ingeo,angles,rank)
  End_time = time.time()
  duration_time = End_time - Start_time
  print('rank,',rank,'Ax INGEO form time =',duration_time)   
  print('rank,',rank,'Ax INGEO y2.sum', np.sum(y2))
if rank ==0:
  tigre.plotimg(y-y2,'x')
