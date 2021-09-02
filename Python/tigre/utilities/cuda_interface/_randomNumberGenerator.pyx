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
#   Coded by:          Tomoyuki SADAKANE

cimport numpy as np 
import numpy as np
from tigre.utilities.errors import TigreCudaCallError
from tigre.utilities.cuda_interface._gpuUtils cimport GpuIds as c_GpuIds, convert_to_c_gpuids, free_c_gpuids

np.import_array()

from libc.stdlib cimport malloc, free 

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
    void PyArray_CLEARFLAGS(np.ndarray arr, int flags)

cdef extern from "RandomNumberGenerator.hpp":
    cdef void poisson_1d(float* img, size_t uiLen, float* dst, c_GpuIds gpuids)
    cdef void poisson_gaussian_1d(float* img, size_t uiLen, float fGaussMu, float fGaussSigma, float* dst, c_GpuIds gpuids)

def cuda_raise_errors(error_code):
    if error_code:
        raise TigreCudaCallError('RandomNumberGenerator::', error_code)

def add_poisson(np.ndarray[np.float32_t, ndim=3] src, gpuids=None):
    # print("add_poisson()")
    cdef c_GpuIds* c_gpuids = convert_to_c_gpuids(gpuids)
    if not c_gpuids:
        raise MemoryError()

    cdef np.npy_intp size_img[3]
    size_img[0]= <np.npy_intp> src.shape[0]
    size_img[1]= <np.npy_intp> src.shape[1]
    size_img[2]= <np.npy_intp> src.shape[2]

    cdef float* c_imgout = <float*> malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float))

    cdef long imgsize[3]
    imgsize[0] = <long> size_img[2]
    imgsize[1] = <long> size_img[1]
    imgsize[2] = <long> size_img[0]

    cdef float* c_src = <float*> src.data
    cdef np.npy_intp c_uiSigLen = <np.npy_intp> (size_img[0] *size_img[1] *size_img[2])
    cuda_raise_errors(poisson_1d(c_src, c_uiSigLen, c_imgout, c_gpuids[0]))
    imgout = np.PyArray_SimpleNewFromData(3, size_img, np.NPY_FLOAT32, c_imgout)
    PyArray_ENABLEFLAGS(imgout, np.NPY_OWNDATA)
    
    return imgout

def add_noise(np.ndarray[np.float32_t, ndim=3] poisson_lambda,
              np.float32_t gaussian_mu,
              np.float32_t gaussian_sigma,
              gpuids=None):
    # print("add_noise()")
    cdef c_GpuIds* c_gpuids = convert_to_c_gpuids(gpuids)
    if not c_gpuids:
        raise MemoryError()

    cdef np.npy_intp size_img[3]
    size_img[0]= <np.npy_intp> poisson_lambda.shape[0]
    size_img[1]= <np.npy_intp> poisson_lambda.shape[1]
    size_img[2]= <np.npy_intp> poisson_lambda.shape[2]

    cdef float* c_imgout = <float*> malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float))

    cdef long imgsize[3]
    imgsize[0] = <long> size_img[2]
    imgsize[1] = <long> size_img[1]
    imgsize[2] = <long> size_img[0]

    cdef float* c_src = <float*> poisson_lambda.data
    cdef np.npy_intp c_uiSigLen = <np.npy_intp> (size_img[0] *size_img[1] *size_img[2])
    cdef np.float32_t c_fGaussMu = gaussian_mu
    cdef np.float32_t c_fGaussSigma = gaussian_sigma
    cuda_raise_errors(poisson_gaussian_1d(c_src, c_uiSigLen, c_fGaussMu, c_fGaussSigma, c_imgout, c_gpuids[0]))
    imgout = np.PyArray_SimpleNewFromData(3, size_img, np.NPY_FLOAT32, c_imgout)
    PyArray_ENABLEFLAGS(imgout, np.NPY_OWNDATA)
    
    return imgout
