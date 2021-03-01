#   This file is part of the TIGRE Toolbox

#   Copyright (c) 2015, University of Bath and
#                       CERN-European Organization for Nuclear Research
#                       All rights reserved.

#   License:            Open Source under BSD.
#                       See the full license at
#                       https://github.com/CERN/TIGRE/license.txt

#   Contact:            tigre.toolbox@gmail.com
#   Codes:              https://github.com/CERN/TIGRE/
# --------------------------------------------------------------------------
#   Coded by:          MATLAB (original code): Ander Biguri
#                      PYTHON : Sam Loescher, Reuben Lindroos

cimport numpy as np 
import numpy as np
from tigre.utilities.errors import TigreCudaCallError
from tigre.utilities.cuda_interface._gpuUtils cimport GpuIds as c_GpuIds, convert_to_c_gpuids, free_c_gpuids

np.import_array()

from libc.stdlib cimport malloc, free 

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
    void PyArray_CLEARFLAGS(np.ndarray arr, int flags)

cdef extern from "tvdenoising.hpp":
    cdef void tvdenoising(float* src, float* dst, float lamda, float* spacing, long* image_size, int maxiter, c_GpuIds gpuids)


def cuda_raise_errors(error_code):
    if error_code:
        raise TigreCudaCallError('tvdenoising:',error_code)


def tvdenoise(np.ndarray[np.float32_t, ndim=3] src, int maxiter = 100, float lamda = 15.0, gpuids = None):
    cdef c_GpuIds* c_gpuids = convert_to_c_gpuids(gpuids)
    if not c_gpuids:
        raise MemoryError()
    
    
    cdef np.npy_intp size_img[3]
    size_img[0]= <np.npy_intp> src.shape[0]
    size_img[1]= <np.npy_intp> src.shape[1]
    size_img[2]= <np.npy_intp> src.shape[2]

    cdef float* c_imgout = <float*> malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float))

    cdef float spacing[3]
    spacing[0]=<float> 1
    spacing[1]=<float> 1
    spacing[2]=<float> 1

    cdef long imgsize[3]
    imgsize[0] = <long> size_img[2]
    imgsize[1] = <long> size_img[1]
    imgsize[2] = <long> size_img[0]

    cdef float* c_src = <float*> src.data
    cdef np.npy_intp c_maxiter = <np.npy_intp> maxiter
    cuda_raise_errors(tvdenoising(c_src, c_imgout, lamda, spacing, imgsize, c_maxiter, c_gpuids[0]))

    imgout = np.PyArray_SimpleNewFromData(3, size_img, np.NPY_FLOAT32, c_imgout)
    PyArray_ENABLEFLAGS(imgout, np.NPY_OWNDATA)

    return imgout
