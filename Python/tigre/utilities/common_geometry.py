import numpy as np

# Tomosymthesis
def staticDetectorGeo(geo,angles,rot=0):
    # angles: angle off the aixs between the source and detector centre (on x-y plane)
    #         when angles=0, the source is perpendicular to the detector
    # rot: rotation of both source and detector around origin
    R = geo.DSD-geo.DSO
    rot = rot/180*np.pi
    geo.DSD = geo.DSO + R*np.cos(angles)
    geo.offDetector = np.vstack([0*angles, R*np.sin(angles)]).T
    geo.rotDetector = np.vstack([0*angles, 0*angles, -angles]).T
    angles += rot
    geo.angles = np.vstack([angles, 0*angles, 0*angles]).T
    
    return geo

# Linear scan of source
def staticDetLinearSourceGeo(geo,s_pos,s_rot=0,rot=0):
    # s_pos: distance along source scanning linear trajectry 
    #        when s_pos = 0, source is aligned to the origin and detector centre
    # s_rot: rotation angle between the source linear trajectry and detector on (x-y) plane 
    # rot: source and detector rotation angle around the origin (on x-y plane)
    if np.isscalar(s_rot):
        s_rot = s_rot * np.ones_like(s_pos)
    s_rot = s_rot / 180 *np.pi
    if np.isscalar(rot):
        rot = rot * np.ones_like(s_pos)
    rot = rot/180*np.pi 
               
    ang  = np.arctan2(s_pos*np.cos(s_rot), geo.DSO + s_pos*np.sin(s_rot))
    R = geo.DSD - geo.DSO
    geo.offDetector = np.array([ang*0, R*np.sin(ang)]).T
    geo.rotDetector = np.array([0*ang, 0*ang, -ang]).T    
    geo.DSO = np.sqrt((s_pos*np.cos(s_rot))**2 + (geo.DSO + s_pos*np.sin(s_rot))**2)
    geo.DSD = geo.DSO + R*np.cos(ang)
    ang += rot        
    geo.angles = np.vstack([ang, 0*ang, 0*ang]).T
    return geo, ang
