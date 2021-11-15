from __future__ import division
from __future__ import print_function
from numpy.core.arrayprint import dtype_is_implied
from tigre.utilities.parkerweight import parkerweight
import numpy as np
from scipy.fft  import fft, ifft

import warnings

import numpy as np
from tigre.utilities.parkerweight import parkerweight


# TODO: Fix parker
def filtering(proj, geo, angles, parker, verbose=False):
    if parker:
        proj=parkerweight(proj.transpose(0,2,1),geo,angles,parker).transpose(0,2,1)

    filt_len=max(64,2**nextpow2(2*max(geo.nDetector)))
    ramp_kernel=ramp_flat(filt_len)

    d=1
    filt=filter(geo.filter,ramp_kernel[0],filt_len,d,verbose=verbose)
    filt=np.kron(np.ones((np.int64(geo.nDetector[0]),1),dtype=np.float32),filt)

    padding = int((filt_len-geo.nDetector[1])//2 )
    scale_factor = (geo.DSD[0]/geo.DSO[0]) * (2 * np.pi/ len(angles)) / ( 4 * geo.dDetector[0] ) 

    #filter 2 projection at a time packing in to complex container
    fproj=np.empty((geo.nDetector[0],filt_len),dtype=np.complex64)
    for i in range(0,angles.shape[0]-1,2):
        fproj.fill(0)
        fproj.real[:,padding:padding+geo.nDetector[1]]=proj[i]
        fproj.imag[:,padding:padding+geo.nDetector[1]]=proj[i+1]

        fproj=fft(fproj,axis=1)
        fproj=fproj*filt
        fproj=ifft(fproj,axis=1)

        proj[i]=fproj.real[:,padding:padding+geo.nDetector[1]] * scale_factor
        proj[i+1]=fproj.imag[:,padding:padding+geo.nDetector[1]] * scale_factor

    #if odd number of projections filter last solo
    if angles.shape[0] % 2:
        fproj.fill(0)
        fproj.real[:,padding:padding+geo.nDetector[1]]=proj[angles.shape[0]-1]

        fproj=fft(fproj,axis=1)
        fproj=fproj*filt
        fproj=np.real(ifft(fproj,axis=1))     
        proj[angles.shape[0]-1]=fproj[:,padding:padding+geo.nDetector[1]] * scale_factor

    return proj


def ramp_flat(n, verbose=False):
    nn = np.arange(-n / 2, n / 2)
    h = np.zeros(nn.shape, dtype=np.float32)
    h[int(n / 2)] = 1 / 4
    odd = nn % 2 == 1
    h[odd] = -1 / (np.pi * nn[odd]) ** 2
    return h, nn


def filter(filter, kernel, order, d, verbose=False):
    f_kernel = abs(np.fft.fft(kernel)) * 2

    filt = f_kernel[: int((order / 2) + 1)]
    w = 2 * np.pi * np.arange(len(filt)) / order

    if filter in {"ram_lak", None}:
        if filter is None and verbose:
            warnings.warn("no filter selected, using default ram_lak")
    elif filter == "shepp_logan":
        filt[1:] *= np.sin(w[1:] / (2 * d)) / (w[1:] / (2 * d))
    elif filter == "cosine":
        filt[1:] *= np.cos(w[1:] / (2 * d))
    elif filter == "hamming":
        filt[1:] *= 0.54 + 0.46 * np.cos(w[1:] / d)
    elif filter == "hann":
        filt[1:] *= (1 + np.cos(w[1:]) / d) / 2
    else:
        raise ValueError("filter not recognised: " + str(filter))

    filt[w > np.pi * d] = 0
    filt = np.hstack((filt, filt[1:-1][::-1]))
    return filt.astype(np.float32)


def nextpow2(n):
    i = 1
    while (2 ** i) < n:
        i += 1
    return i
