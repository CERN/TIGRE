import numpy as np
import copy
from scipy.spatial.transform import Rotation

# Tomosymthesis
def staticDetectorGeo(geo,angles,rot=0):
    """
    # angles: angle off the aixs between the source and detector centre (on x-y plane)
    #         when angles=0, the source is perpendicular to the detector
    # rot: rotation of both source and detector around origin
    """
    ngeo = copy.deepcopy(geo)
    R = ngeo.DSD-ngeo.DSO
    rot = rot/180*np.pi
    ngeo.DSD = ngeo.DSO + R*np.cos(angles)
    ngeo.offDetector = np.vstack([0*angles, R*np.sin(angles)]).T
    ngeo.rotDetector = np.vstack([0*angles, 0*angles, -angles]).T
    angles += rot
    ngeo.angles = np.vstack([angles, 0*angles, 0*angles]).T
    
    return ngeo

# Linear scan of source
def staticDetLinearSourceGeo(geo,s_pos,s_rot=0,rot=0):
    """
    # s_pos: distance along source scanning linear trajectry 
    #        when s_pos = 0, source is aligned to the origin and detector centre on x-axis
    # s_rot: rotation angle between the source linear trajectry and detector on (looking from top on x-y plane, anticlock-wise) 
    # rot: source and detector rotation angle around the origin (looking from top on x-y plane, anit-clockwise) 
    """
    ngeo = copy.deepcopy(geo)
    if np.isscalar(s_rot):
        s_rot = s_rot * np.ones_like(s_pos)
    s_rot = s_rot / 180 *np.pi
    if np.isscalar(rot):
        rot = rot * np.ones_like(s_pos)
    if len(s_pos) != len(rot) or len(s_pos) != len(s_rot):
        raise("Inputs length do not match")
    
    rot = rot/180*np.pi 
               
    ang  = np.arctan2(s_pos*np.cos(s_rot), ngeo.DSO + s_pos*np.sin(s_rot)) 
    R = ngeo.DSD - ngeo.DSO
    if hasattr(ngeo,"offDetector"):
        ngeo.offDetector = ngeo.offDetector.astype(np.float32) + np.array([ang*0, R*np.sin(ang)]).T 
    else:
        ngeo.offDetector = np.array([ang*0, R*np.sin(ang)]).T 
    if hasattr(ngeo,"rotDetector"):    
        ngeo.rotDetector += np.array([0*ang, 0*ang, -ang]).T   
    else:
        ngeo.rotDetector = np.array([0*ang, 0*ang, -ang]).T 
    ngeo.DSO = np.sqrt((s_pos*np.cos(s_rot))**2 + (ngeo.DSO + s_pos*np.sin(s_rot))**2)
    ngeo.DSD = np.sqrt((s_pos*np.cos(s_rot))**2 + (ngeo.DSD + s_pos*np.sin(s_rot))**2)
    ang += rot    
    ngeo.angles = np.vstack([ang, 0*ang, 0*ang]).T
    return ngeo

def ArbitrarySourceDetMoveGeo(geo,s_pos,d_pos=None,d_rot=None):
    """
    # Source and Detector can move arbitrarily while the object is fixed
    #
    #   Parameters
    #   ----------
    #   geo:   standard cone beam geometry
    #   s_pos: nx3 array, source movement coordinates (x,y,z), in mm
    #           default s_pos = (DSO, 0, 0), source is on x-axis (+) 
    #   d_pos: nx3 array, detector centre movement coordinates (x,y,z), in mm 
    #           default d_pos = (DSO-DSD, 0, 0), detector centre on x-axis (-)
    #   d_rot: nx3 array, detector rotation angles (roll, pitch, yaw) in degrees, 
    #           default - no rotation, detector facing origin
    #   Note:
    #          a point source rotation has no effect, ignored
    #
    #   Returns
    #   -------
    #   geometry with arbitrarily specified movements of source and detector
    #
    """
    s_pos = np.array(s_pos) if isinstance(s_pos,list) else s_pos
    ngeo = copy.deepcopy(geo)
    if s_pos.ndim != 2:
        raise("Input s_pos should be an n x 3 array")
    else:
        n = s_pos.shape[0]
    if s_pos.shape[1] != 3:
        raise("Input s_pos should be an n x 3 array")
    if d_pos is None:
        d_pos = np.repeat([[ngeo.DSO-ngeo.DSD, 0, 0]], n, axis=0)
    elif isinstance(d_pos,list):
        d_pos = np.array(d_pos)
    if d_rot is None:
        d_rot = np.zeros((n,3))
    elif isinstance(d_rot,list):
        d_rot = np.array(d_rot)
    if s_pos.shape != d_pos.shape or s_pos.shape != d_rot.shape:
        raise("Inputs dimensions do not match")        
    
    d_rot = d_rot / 180*np.pi
    
    # source and detector vector lengths and directions
    rs = np.linalg.norm(s_pos, axis=1)
    ns = s_pos / np.tile(rs, (3,1)).T
    rd = np.linalg.norm(d_pos, axis=1)
    nd = d_pos / np.tile(rd, (3,1)).T
    
    # source euler angles (intrinsic) away from x-axis in "ZYZ" order
    s0 = np.zeros_like(s_pos)
    s0[:,0] = s_pos[:,0]
    s_ang = euler_from_vecs(s0, s_pos, order="zyz")   
    
    # check angles between OS and OD
    sd_cross = np.cross(ns, nd)
    # note: A x B = ||A||.||B||sin(theta)
    sd_ang = np.arcsin(sd_cross)
    if (abs(sd_ang) >= np.pi/2).any():  
        raise RuntimeError("Source and detector are not always on different sides (some angle<=90 degree)")        
    
    # detector angles from ray
    d_cross = np.cross(-ns, nd) 
    d_ang = np.arcsin(d_cross)
    pd = d_cross * np.tile(rd, (3,1)).T  
    
    if hasattr(ngeo,"offDetector"):
        ngeo.offDetector = ngeo.offDetector.astype(np.float64) + np.vstack([pd[:,1], -pd[:,2]]).T
    else:
        ngeo.offDetector = np.vstack([pd[:,1], -pd[:,2]]).T
        
    if hasattr(ngeo,"rotDetector"):
        ngeo.rotDetector += d_ang + d_rot
    else:
        ngeo.rotDetector = d_ang + d_rot
    ngeo.DSO = np.linalg.norm(s_pos, axis=1) 
    ngeo.DSD = np.linalg.norm(s_pos-d_pos, axis=1) 
    ngeo.angles = s_ang 
    return ngeo

    
def euler_from_vecs(a_vec,b_vec,order="xyz"):
    """
    Calculate Euler angles from two vectors

    Parameters
    ----------
    a_vec : np.array of (n,3) or (3,)
        A vector(s).
    b_vec : np.array of (n,3) or (3,)
        B vector(s).
    order : string, optional
        Order of the euler angles in rotations. The default is "xyz".

    Returns
    -------
    euler : np.array of (n,3) or (3,)
        Euler angles of rotations.

    """
    na = 1 if a_vec.ndim<2 else a_vec.shape[0]
    nb = 1 if b_vec.ndim<2 else b_vec.shape[0]
    n = max(na,nb)
    a_vec = np.repeat([a_vec], n, axis=0) if a_vec.ndim<2 else a_vec
    b_vec = np.repeat([b_vec], n, axis=0) if b_vec.ndim<2 else b_vec
    euler = np.zeros_like(a_vec,dtype=np.float64)
    for i in range(n):
        R = rotation_from_vecs(a_vec[i,:],b_vec[i,:])
        euler[i,:] = Rotation.from_matrix(R).as_euler(order)
        for j in range(3):
            if abs(euler[i,j]) < 2e-8: # very small angle
                euler[i,j] = 0 
    return euler


def rotation_from_vecs(v1, v2):
    """
    Compute a matrix R that rotates v1 to align with v2.
    v1 and v2 must be length-3 1d numpy arrays.
    """
    # unit vectors
    u = v1 / np.linalg.norm(v1)
    Ru = v2 / np.linalg.norm(v2)
    # dimension of the space and identity
    I = np.identity(u.size)
    # the cos angle between the vectors (dot product)
    c = np.dot(u, Ru)
    # a small number
    eps = 1.0e-9
    if np.abs(c - 1.0) < eps:
        # same direction
        return I
    elif np.abs(c + 1.0) < eps:
        # opposite direction
        return -I
    else:
        # the cross product matrix of a vector to rotate around
        v = np.cross(u, Ru)
        K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        # Rodrigues' formula
        return I + K + (K @ K) / (1 + c)    

