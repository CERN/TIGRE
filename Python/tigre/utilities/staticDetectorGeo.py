import numpy as np

def staticDetectorGeo(geo,angles):
        R = geo.DSD-geo.DSO
        geo.DSD -= R*(1-np.cos(angles))
        geo.offDetector = np.vstack([0*angles, R*np.sin(angles)]).T
        geo.rotDetector = np.vstack([0*angles, 0*angles, -angles]).T
        geo.angles = np.vstack([angles, 0*angles, 0*angles]).T
        return geo
