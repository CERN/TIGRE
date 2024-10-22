import math
import numpy as np

def flatten_detector(proj, geo, oversample=1,arclength=None):
    
    if arclength is None:
        arclength=geo.sDetector[1]
    (num_angles, num_detector_rows, orig_num_detectors) = np.shape(proj)
    pixel_arclength = arclength / orig_num_detectors

    # Calculate the size of the detector projected from the arc onto the plane
    # This is slightly different for odd and even detectors
    if orig_num_detectors % 2:
        # Odd number of detectors - one on the centreline
        detector_size = math.tan(pixel_arclength) * geo.DSD
    else:
        detector_size = math.tan(pixel_arclength / 2) * geo.DSD * 2

    detector_size = detector_size / oversample
    total_detector_size = math.tan(arclength / 2) * geo.DSD * 2

    # Calculate the equal distance detector positions on the plane detector
    if orig_num_detectors % 2:
        # Odd number of detectors - one on the centreline
        detectors = np.arange(detector_size, total_detector_size / 2, detector_size)
        detectors = np.concatenate((np.flip(-detectors), [0], detectors))
    else:
        detectors = np.arange(detector_size / 2, total_detector_size / 2, detector_size)
        detectors = np.concatenate((np.flip(-detectors), detectors))

    num_detector_columns = len(detectors)

    # Calculate the angle from the source to detector positions on the plane detector
    angles = np.arctan2(detectors, geo.DSD)

    # Normailse to give this angle relative to the projection detector
    normalised_angles = angles / pixel_arclength + orig_num_detectors / 2

    # Interpolate for projected detectors
    flattened_proj = np.zeros((num_angles, num_detector_rows, num_detector_columns), dtype=np.float32)
    for i in range(num_angles):
        for r in range(num_detector_rows):
            flattened_proj[i, r, :] = np.interp(normalised_angles, np.arange(orig_num_detectors), proj[i, r, :])

    # Update geometry
    geo.nDetector = np.array([num_detector_rows, num_detector_columns])
    geo.dDetector = np.array([geo.dVoxel[0], detector_size])
    geo.sDetector = geo.nDetector * geo.dDetector

    return flattened_proj

