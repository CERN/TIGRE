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

cdef extern from "FbpFiltration.hpp":
    cdef void apply_filtration(const float* pf32In, size_t uiXLen, size_t uiYLen, float* pf32Filter, float* pfOut, const c_GpuIds gpuids)

def cuda_raise_errors(error_code):
    if error_code:
        raise TigreCudaCallError('FFT::', error_code)

def apply(np.ndarray[np.float32_t, ndim=2] src, np.ndarray[np.float32_t, ndim=1] flt, gpuids=None):
    cdef c_GpuIds* c_gpuids = convert_to_c_gpuids(gpuids)
    if not c_gpuids:
        raise MemoryError()

    cdef np.npy_intp size_img[2]
    size_img[0]= <np.npy_intp> src.shape[0]
    size_img[1]= <np.npy_intp> src.shape[1]

    cdef float* c_imgout = <float*> malloc(size_img[0] *size_img[1]* sizeof(float))

    src = np.ascontiguousarray(src)
    flt = np.ascontiguousarray(flt)

    cdef float* c_src = <float*> src.data
    cdef float* c_flt = <float*> flt.data
    cdef size_t uiXLen = size_img[1]
    cdef size_t uiYLen = size_img[0]
    cuda_raise_errors(apply_filtration(c_src, uiXLen, uiYLen, c_flt, c_imgout, c_gpuids[0]))
    imgout = np.PyArray_SimpleNewFromData(2, size_img, np.NPY_FLOAT32, c_imgout)
    PyArray_ENABLEFLAGS(imgout, np.NPY_OWNDATA)
    return imgout
