
import numpy as np
from matplotlib import pyplot as plt

from numpy import sin, cos
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import animation
import tigre

   
def plot_geometry(geo,angles=np.linspace(0,2*np.pi,100),pos=0,animate=False,fname=None):
#PLOT_GEOMETRY(GEO,ANGLES,POS,Animate,fname) plots a simplified version of the CBCT geometry with the
# given geomerty GEO and scanning angles ANGLES at angle POS. If angles is 
# not given, [0,2pi] at 100 steps will be given, and pos=0 will be chosen. If 
# Animate=TRUE, an animation file will be generated. Animation file name can be
# specified by fname.
# 
# h=PLOT_GEOMETRY(...) will return the figure handle 
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
# Coded by:           Ander Biguri, modified by Yi Liu
#--------------------------------------------------------------------------

    # check geo, which makes DSD, DSO etc all matrices
    try:
        geo.check_geo(angles)
    except:
        print('Please check that geometry and angles are consistent')
        
    if pos < angles[0] or pos > angles[-1]:
        raise ValueError('Pos should be within the range of [{}, {}]'.format(angles[0],angles[-1]))
        
    pos = abs(angles-pos).argmin()
    
    ## Figure stuff
    fig=plt.figure(figsize=(10,8))
    ax=plt.axes(projection='3d')
    ax.view_init(azim=52,elev=26)
           
    ln = geo.sVoxel.min()   
    
    ## source trajectory, cordinates in (z,y,x)
    thz,thx,thy = geo.angles[:,0], geo.angles[:,1], geo.angles[:,2]
    Rs = rotmat3d(geo.angles,order='ZYZ')  # source and detector rotation are in opposite direction to that of object
    stj = np.zeros_like(geo.angles)
    # Note: angles = (0,0,0) <==> source at (x=DSO, y=0, z=0)
    scent = np.array([geo.DSO,0*thy,0*thz]).T # source centre (x,y,z) before scan rotation
    for j in range(len(thx)): 
        stj[j,:] = np.matmul(Rs[:,:,j], scent[j,:])  # no offset for source
    # displacement in y for geo.COR
    if hasattr(geo, 'COR'):
        stj[:,1] += geo.COR
    if np.ptp(stj,axis=0).any():
        ax.plot3D(stj[:,0],stj[:,1],stj[:,2],color='grey',ls='',marker=".",markersize=2.5,mfc='grey',mec="grey")       
    # source centre at pos
    source=ax.scatter(stj[pos,0],stj[pos,1],stj[pos,2],color='r',s=5)
    stext=ax.text(stj[pos,0],stj[pos,1],stj[pos,2]+15,'S',None)
    
    ## detector trajectory: scan -> offset -> rotate
    # detector centre (x,y,z) scan rotation. 
    dcent = np.array([-geo.DSD+geo.DSO, geo.offDetector[:,1], geo.offDetector[:,0]]).T
    dtj = np.zeros_like(geo.angles)
    for j in range(len(thx)):
        dtj[j,:] = np.matmul(Rs[:,:,j], dcent[j,:]) 
    # displacement in y for geo.COR
    if hasattr(geo, 'COR'):
        dtj[:,1] += geo.COR
    # # detector offset
    if np.ptp(dtj,axis=0).any():
        ax.plot3D(dtj[:,0],dtj[:,1],dtj[:,2],color='grey',ls='',marker=".",markersize=2.5,mfc='grey',mec="grey")
    # detector centre at pos
    det=ax.scatter(dtj[pos,0],dtj[pos,1],dtj[pos,2],color='brown',s=5)
    dtext=ax.text(dtj[pos,0],dtj[pos,1],dtj[pos,2]+15,'D',None)
    # after scan to pos, then rotate the detector (roll,pitch,yaw) <==> (x,y,z)
    R = np.zeros_like(Rs)
    for j in range(Rs.shape[2]):
        R[:,:,j] = np.matmul( Rs[:,:,j], rotmat3d(geo.rotDetector[j,:],order='XYZ'))
    # Detector, cp returns four cordinates of corners closest to source
    ddp = 5 # detector depth
    dsz = np.array([ddp,geo.sDetector[1],geo.sDetector[0]])  # at angles (0,0,0)
    dverts = calCube(dtj,dsz,R)
    dcube = art3d.Poly3DCollection(dverts[pos],color='brown',alpha=0.3)
    ax.add_collection3d(dcube)
        
    ## origin trajectory. NOTE: offOrigin is in (z,y,x)
    otj = np.fliplr(geo.offOrigin).astype(np.float32)
    # displacement in y for geo.COR
    if hasattr(geo, 'COR'):
        otj[:,1] += geo.COR   
    if np.ptp(otj,axis=0).any():
        ax.plot3D(otj[:,0],otj[:,1],otj[:,2],color='grey',ls='',marker=".",markersize=2.5,mfc='grey',mec="grey")        
    # Cordinates Arrows from origin
    ax.quiver(0,0,0,1,0,0,length=ln,color='r')
    ax.quiver(0,0,0,0,1,0,length=ln,color='b')
    ax.quiver(0,0,0,0,0,1,length=ln,color='g')
    ax.text(-10,-10,-10,'O',None)
    # CUBE/Image, sVoxel dims order (z,y,x) in python
    overts = calCube(otj,geo.sVoxel[[2,1,0]])                              
    ocube = art3d.Poly3DCollection(overts[pos],color='c',alpha=0.1)
    ax.add_collection3d(ocube)
    
    # effective beam profile
    beampf = [0,0,0,0]
    for i in range(4):
        beampf[i]=ax.plot3D(*zip(stj[pos,:],dverts[pos][0][i]),color='y')
    # cenrtal beam
    cbeam=ax.plot3D(*zip(stj[pos,:],dtj[pos,:]),color='pink')
    
    # set tight limits and aspect
    limX = ax.get_xlim3d()
    limY = ax.get_ylim3d()
    limZ = ax.get_zlim3d()
    ax.set_box_aspect((limX[1]-limX[0], limY[1]-limY[0], limZ[1]-limZ[0]))

    # set up plot
    roll, pitch, yaw = geo.rotDetector[:,0]/np.pi*180, geo.rotDetector[:,1]/np.pi*180, geo.rotDetector[:,1]/np.pi*180
    plt.title('Current CBCT geometry at angles {} (deg), in scale'.format(list(map('{:.1f}'.format,geo.angles[pos]/np.pi*180))) \
             +'\n\nDetector rotation [{:.1f} - {:.1f}, {:.1f} - {:.1f}, {:.1f} - {:.1f}] (deg)'.format(roll[0],roll[-1],pitch[0],pitch[-1],yaw[0],yaw[-1]))       
    ax.set_xlabel('X');
    ax.set_ylabel('Y');
    ax.set_zlabel('Z');
          
    ax.view_init(azim=52,elev=26)
    plt.tight_layout()

    def update(pos,stj,dtj,otj,dverts,overts):
        source.set_offsets(stj[pos,:2])
        source.set_3d_properties(stj[pos,2], 'z')
        stext.set_position(stj[pos,:2]+15)
        det.set_offsets(dtj[pos,:2])
        det.set_3d_properties(dtj[pos,2], 'z')
        dtext.set_position(dtj[pos,:2]+15)
        cbeam[0].set_data(*zip(stj[pos,:2],dtj[pos,:2]))
        cbeam[0].set_3d_properties((stj[pos,2],dtj[pos,2]),'z')
        for i in range(4):
            beampf[i][0].set_data(*zip(stj[pos,:2],dverts[pos][0][i][:2]))
            beampf[i][0].set_3d_properties((stj[pos,2],dverts[pos][0][i][2]),'z')
        dcube.set_verts(dverts[pos])
        ocube.set_verts(overts[pos])
        ax.set_title('Current CBCT geometry at angles {} (deg), in scale'.format(list(map('{:.1f}'.format,geo.angles[pos]/np.pi*180))) \
                 +'\n\nDetector rotation [{:.1f} - {:.1f}, {:.1f} - {:.1f}, {:.1f} - {:.1f}] (deg)'.format(roll[0],roll[-1],pitch[0],pitch[-1],yaw[0],yaw[-1]))
    

    if animate:
        ani = animation.FuncAnimation(fig, update, len(angles), fargs=(stj,dtj,otj,dverts,overts), interval=100)
        if isinstance(fname, str):
            fname += '_geometry'
            try:
                ani.save('%s.mp4'%(fname),writer='ffmpeg',fps=30)
            except ValueError:
                print('Movie writer "ffmpeg" unavailable, try "Pillow"')
                try:
                    ani.save('%s.gif'%(fname),writer='pillow',fps=30)
                except ValueError:
                    print('Movie writer "Pillow" unavailable, try "ImageMagick"')
                    try: 
                        ani.save('%s.gif'%(fname),writer='imagemagick',fps=30)
                    except ValueError:
                        print('Movie writer unavailable, animation not saved.')
        return ani
    else:
        return ax

# other useful functions
def calCube(centre,size,R='eye'):
    # Calculate vertices of a cuboid, centred at "centre" with "size" and rotated by R, 3D rotation matrix 
    # centre is either (n,3) or (3,) array and R is (3,3,n) or (3,3) array  
    ## 3D cube corners (x,y,z)
    corner = np.array([[-1, -1, -1],
                       [1, -1, -1],
                       [1, 1, -1],
                       [-1, 1, -1],
                       [-1, -1, 1],
                       [1, -1, 1],
                       [1, 1, 1],
                       [-1, 1, 1]])   
        
    if centre.ndim>1:
        n = centre.shape[0]
        verts = []
        for i in range(n):
            # scaling to cuboid and rotation
            if isinstance(R, str):
                c = np.matmul( corner*size/2, np.eye(3) ) + centre[i,:]
            else:
                c = np.matmul( corner*size/2, R[:,:,i].T ) + centre[i,:]
        
            # vertices of the cube
            verts.append( [[c[1], c[2], c[6], c[5]],   # front
                           [c[0], c[1], c[2], c[3]],   # bottom
                           [c[4], c[5], c[6], c[7]],   # top
                           [c[0], c[1], c[5], c[4]],   # left
                           [c[2], c[3], c[7], c[6]],   # right
                           [c[0], c[3], c[7], c[4]]] ) # back           
    else:
        # scaling and rotation
        if isinstance(R, str):
            c = np.matmul( corner*size/2, np.eye(3) ) + centre
        else:    
            c = np.matmul( corner*size/2, R.T ) + centre
    
        # vertices of the cube
        verts = [[c[1], c[2], c[6], c[5]],   # front
                 [c[0], c[1], c[2], c[3]],   # bottom
                 [c[4], c[5], c[6], c[7]],   # top
                 [c[0], c[1], c[5], c[4]],   # left
                 [c[2], c[3], c[7], c[6]],   # right
                 [c[0], c[3], c[7], c[4]]]   # back
   
    return verts


def rotmat3d(ang,order='ZYX'):
    # Calculate 3D rotation matrices of any Eular configuation (right-hand)
    # ang: can be (3,) or (n,3) angles in rads
    # order: str of any combination of "X", "Y", "Z", in any order
    # cordinates in (x,y,z)
    
    if ang.ndim == 1:
        an = np.expand_dims(ang,0)
    else:
        an = ang.copy()
        
    r = np.eye(3)
    R = np.stack([r for _  in range(an.shape[0])], axis=2)

    for j in range(an.shape[0]):
        for i in range(len(order)):
            if order[i] in 'X':
                R[:,:,j] = np.matmul(rotX(an[j,i]), R[:,:,j])
            elif order[i] in 'Y':
                R[:,:,j] = np.matmul(rotY(an[j,i]), R[:,:,j])
            elif order[i] in 'Z':
                R[:,:,j] = np.matmul(rotZ(an[j,i]), R[:,:,j])          
    return np.squeeze(R)

    def do_3d_projection(self, render=None):
        from mpl_toolkits.mplot3d import proj3d
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def rotZ(a):
    # rotation around z axis of a rads (yaw)
    Rz = np.array([[cos(a),-sin(a),0*a],
                   [sin(a),cos(a),0*a],
                   [0*a,0*a,0*a+1]])
    return Rz

def rotY(b):
    # rotation around y axis of b rads (pitch)
    Ry = np.array([[cos(b),0*b,sin(b)],
                   [0*b,0*b+1,0*b],
                   [-sin(b),0*b,cos(b)]])
    return Ry

def rotX(c):
    # rotation aroud x axis of c rads (roll)
    Rx = np.array([[0*c+1,0*c,0*c],
                   [0*c,cos(c),-sin(c)],
                   [0*c,sin(c),cos(c)]])
    return Rx
 
if __name__ == "__main__":
    n_voxel_z = 128
    n_voxel_y = 256
    n_voxel_x = 512
    n_detector_u = 400
    n_detector_v = 300
    off_detector_u = 0  # 500
    off_detector_v = 0  # 500
    off_origin_x = 0  # 300
    off_origin_y = 0  # 100
    off_origin_z = 0  # 100

    geo = tigre.geometry(mode="cone", default=True)
    geo.nVoxel = np.array([n_voxel_z, n_voxel_y, n_voxel_x])
    geo.sVoxel = geo.nVoxel
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.nDetector = np.array([n_detector_v, n_detector_u])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offDetector = np.array([off_detector_v, off_detector_u])
    geo.offOrigin = np.array([off_origin_z, off_origin_y, off_origin_x])
    print(geo)
    angles = np.linspace(0.0,np.pi,180)
    pos = np.pi / 6
    ani=plot_geometry(geo, angles, pos, animate=True)
    plt.show()
