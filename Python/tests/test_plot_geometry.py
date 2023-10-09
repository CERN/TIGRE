
## test    

import numpy as np
from matplotlib import pyplot as plt
import tigre

# # plot 1, check zero angle position, scan rotation direction, detector offset

# geo = tigre.geometry(nVoxel=np.array([256,128,256]),default=True)
# geo.sVoxel = geo.nVoxel * np.array([1,1,1]) # (z,y,x)
# geo.nDetector = np.array([256,128])
# geo.dDetector = np.array([0.8, 0.8])*2               
# geo.sDetector = geo.dDetector * geo.nDetector 
# geo.offDetector=np.array([0,0]) # viewing from S, D move (up, right) => (v, u)
# geo.rotDetector=np.array([30,0,0])/180*np.pi # [roll, pitch, yaw] viewing from S to D
# geo.offOrigin = np.array([0,0,0]) # (z,y,x)
# geo.COR=0
# angles=np.linspace(0,np.pi,100)
# ani1=tigre.plot_geometry(geo,angles,10,animate=True)  # angle=0, S is at (x=DSO, y=0, z=0)
# ani1
# # confirm the plot with projection
# from scipy.io import loadmat
# head=loadmat('head.mat')['img'].transpose(2,1,0).copy()
# head=head[:,:128,:].copy()
# proj = tigre.Ax(head,geo,angles)
# plt.figure()
# plt.subplot(1,2,1)
# plt.imshow(head[:,:,128],origin='lower')
# plt.title('dim2=128')
# plt.ylabel('dim0 ->')
# plt.xlabel('dim1 ->')
# plt.subplot(1,2,2)
# plt.imshow(proj[0,:,:],origin='lower')
# plt.ylabel('v ->')
# plt.xlabel('u ->')


# plot 2, check staticDetectorGeo() for tomosymthesis setup

geo = tigre.geometry_default()
angles=np.linspace(-60,60,31)/180*np.pi
geo = tigre.staticDetectorGeo(geo,angles,60)

ani2=tigre.plot_geometry(geo,angles,0,animate=True,fname='Tomosynthesis')
ani2

 
## plot 3, fixed target object and detector positions and orientations, source moving linearly

geo = tigre.geometry_default()
df = np.linspace(-510,510,64)  # source position on 
geo.DSO = 750
geo.DSD = 1000

geo, angles = tigre.staticDetLinearSourceGeo(geo,df,10,60)

ani3=tigre.plot_geometry(geo, angles, 0, animate=True, fname='Linear_Tomosynthesis')
ani3


## plot 4, helical CT
geo = tigre.geometry_default(high_resolution=False)

angles = np.linspace(0, 2 * np.pi, 100)
angles = np.hstack([angles, angles, angles])  # loop 3 times

# This makes it helical, axis order (z,y,x) for python
geo.offOrigin = np.zeros((angles.shape[0], 3))
geo.offOrigin[:, 0] = np.linspace(
    -1024 / 2 + 128, 1024 / 2 - 128, angles.shape[0])

ani4 = tigre.plot_geometry(geo, angles, 0, animate=True, fname='Helical_CT')
ani4
