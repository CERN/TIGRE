from __future__ import division
from __future__ import print_function
from numpy.core.arrayprint import dtype_is_implied
from tigre.utilities.parkerweight import parkerweight
import numpy as np
from scipy.fft  import fft, ifft

import os
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from tigre.utilities.parkerweight import parkerweight

# bytes of padded complex data a thread filters at a time, small enough to stay in cache
BLOCK_BYTES = 2**21
# threads used when the caller does not choose; more than this was not faster where measured
MAX_WORKERS = 16


# TODO: Fix parker
def filtering(proj, geo, angles, parker, verbose=False, workers=None):
    if parker:
        proj=parkerweight(proj.transpose(0,2,1),geo,angles,parker).transpose(0,2,1)

    filt_len=max(64,2**nextpow2(2*geo.nDetector[1])) 
    ramp_kernel=ramp_flat(filt_len)

    d=1
    filt=filter(geo.filter,ramp_kernel[0],filt_len,d,verbose=verbose)

    padding = int((filt_len-geo.nDetector[1])//2 )
    scale_factor = (geo.DSD[0]/geo.DSO[0]) * (2 * np.pi/ len(angles)) / ( 4 * geo.dDetector[1] ) 

    nangles=angles.shape[0]
    nrows=int(geo.nDetector[0])
    cols=slice(padding,padding+geo.nDetector[1])
    if workers is None:
        workers=default_workers()

    #detector rows are filtered independently of each other: cut them in blocks that fit
    #in cache and share the blocks among threads (numpy and scipy.fft release the GIL).
    #A multiple of 16 rows, so that the FFT packs the same rows in a SIMD vector as it
    #does on a whole projection, and rounds the same way
    block=max(16,BLOCK_BYTES//(8*filt_len)//16*16)
    local=threading.local()

    def filter_blocks(tasks):
        if not hasattr(local,"fproj"):
            local.fproj=np.empty((block,filt_len),dtype=np.complex64)

        for i,r0,r1 in tasks:
            #filter 2 projection at a time packing in to complex container
            fproj=local.fproj[:r1-r0]
            fproj.fill(0)
            fproj.real[:,cols]=proj[i,r0:r1]
            #if odd number of projections filter last solo
            if i+1<nangles:
                fproj.imag[:,cols]=proj[i+1,r0:r1]

            fproj=fft(fproj,axis=1,overwrite_x=True)
            fproj*=filt
            fproj=ifft(fproj,axis=1,overwrite_x=True)

            proj[i,r0:r1]=fproj.real[:,cols] * scale_factor
            if i+1<nangles:
                proj[i+1,r0:r1]=fproj.imag[:,cols] * scale_factor

    tasks=[(i,r0,min(r0+block,nrows)) for i in range(0,nangles,2) for r0 in range(0,nrows,block)]
    #a detector with few rows gives small tasks: hand them to the threads several at a time
    batch=max(1,block//nrows)
    batches=[tasks[k:k+batch] for k in range(0,len(tasks),batch)]
    workers=min(workers,len(batches))
    if workers>1:
        with ThreadPoolExecutor(workers) as pool:
            #list() so that an exception raised in a thread is raised here
            list(pool.map(filter_blocks,batches))
    else:
        filter_blocks(tasks)

    return proj


def default_workers():
    try:
        ncpu = len(os.sched_getaffinity(0))  # the cores this process is allowed to use
    except AttributeError:  # not available on Windows and macOS
        ncpu = os.cpu_count() or 1
    return min(ncpu, MAX_WORKERS)


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
        filt[1:] *= (1 + np.cos(w[1:] / d)) / 2
    else:
        raise ValueError("filter not recognized: " + str(filter))

    filt[w > np.pi * d] = 0
    filt = np.hstack((filt, filt[1:-1][::-1]))
    return filt.astype(np.float32)


def nextpow2(n):
    i = 1
    while (2 ** i) <= n:
        i += 1
    return i
