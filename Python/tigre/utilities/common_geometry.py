import numpy as np

# Tomosymthesis
def staticDetectorGeo(geo,angles):
    R = geo.DSD-geo.DSO
    geo.DSD -= R*(1-np.cos(angles))
    geo.offDetector = np.vstack([0*angles, R*np.sin(angles)]).T
    geo.rotDetector = np.vstack([0*angles, 0*angles, -angles]).T
    geo.angles = np.vstack([angles, 0*angles, 0*angles]).T
    
    return geo

# Linear scan of source
def staticDetLinearSourceGeo(geo,s_pos,s_rot=0,d_loc='x',d_rot=0):
    # s_pos: distance along source scanning linear trajectry 
    #        when s_pos = 0, source is aligned origin and detector centre
    # s_rot:   rotation angle of the linear trajectry on (x-y) plane 
    # d_loc: the source location on (x-y) plane when s_pos=0,
    if isinstance(d_loc, str):
        if d_loc not in ['x','X','y','Y','+x','+X','+y','+Y','-x','-X','-y','-Y']: 
            raise Exception('Wrong d_loc option, valid choices are {"x" "y", "-x", "-y"}')
    if np.isscalar(s_rot):
        s_rot = s_rot * np.ones_like(s_pos)
    if (abs(s_rot)>90).any():       
        raise Exception('Source scanning rotation angle should be within (-90, 90) degree')
    else: 
        s_rot = s_rot/180*np.pi 
    if np.isscalar(d_rot):
        d_rot = d_rot * np.ones_like(s_pos)
    if (abs(d_rot)>90).any():       
        raise Exception('Detector rotation angle should be within (-90, 90) degree')
    else: 
        d_rot = d_rot/180*np.pi 
               
    ang  = np.arctan2(s_pos*np.cos(s_rot), geo.DSD + s_pos*np.sin(s_rot))
    dshift = s_pos*np.cos(s_rot)*(geo.DSD-geo.DSO)/geo.DSO
    geo.offDetector = np.array([ang*0, dshift]).T
    geo.rotDetector = np.array([0*ang, 0*ang, -ang]).T    
    geo.DSO = np.sqrt((s_pos*np.cos(s_rot))**2 + (geo.DSO + s_pos*np.sin(s_rot))**2)
    geo.DSD = np.sqrt((s_pos*np.cos(s_rot))**2 + (geo.DSD + s_pos*np.sin(s_rot))**2)
    for i in range(len(s_pos)):
        if d_loc[i] in ['y','Y','+y','+Y']:
            ang[i] += np.pi/2 + d_rot[i]
        elif d_loc[i] in ['-y','-Y']: 
            ang[i] -= np.pi/2 + d_rot[i]
        elif d_loc[i] in ['-x','-X']:
            ang[i] += np.pi + d_rot[i]         
    geo.angles = np.vstack([ang, 0*ang, 0*ang]).T
    return geo, ang

