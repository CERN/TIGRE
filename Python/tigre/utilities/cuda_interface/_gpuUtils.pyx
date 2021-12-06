
from tigre.utilities.cuda_interface._types cimport Geometry as c_Geometry, convert_to_c_geometry, free_c_geometry

cimport numpy as np
import numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

from libc.stdlib cimport malloc, free
# from libc.stdio cimport printf, sprintf
# from tigre.utilities.errors import TigreCudaCallError
# cdef extern from "numpy/arrayobject.h":
#     void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef extern from "gpuUtils.hpp":
    cdef int GetGpuIdArray(char* kacGPUName, int* piDeviceIds, int iIdCountMax, char* pcMessage)
cdef extern from "gpuUtils.hpp":
    cdef void GetGpuName(int iDeviceId, char* pcName)
cdef extern from "gpuUtils.hpp":
    cdef int GetGpuCount()
# cdef extern from "gpuUtils.hpp":
#     cdef int GetMaxGpuCount()

# def cuda_raise_errors(error_code):
#     if error_code:
#         raise TigreCudaCallError('Atb:',error_code)

def getGpuIdList(nameGPU = ""):
    # print("getGpuIdList(nameGPU = {})".format(nameGPU))
    cdef int nGPUMax = GetGpuCount()
    cdef int* paIndexes    = <int*> malloc(nGPUMax* sizeof(int))
    cdef char* pcMessage    = <char*> malloc(65535* sizeof(char))
    cdef int iGPUCount = GetGpuIdArray(nameGPU.encode('utf-8'), paIndexes, nGPUMax, pcMessage)
    strMessage = (<bytes>pcMessage).decode('utf-8')
    # print("iGPUCount = {}".format(iGPUCount))
    # print(strMessage)
    if iGPUCount > 0:
        listDevices = [x for x in paIndexes[:iGPUCount]]
    else:
        listDevices = []
    free(paIndexes)
    free(pcMessage)
    return listDevices

def getGpuName(iDeviceId = 0):
    # print("iDeviceId = {}".format(iDeviceId))
    cdef char* pcDeviceName    = <char*> malloc(128* sizeof(char))
    GetGpuName(iDeviceId, pcDeviceName)
    return (<bytes>pcDeviceName).decode('utf-8')

def getGpuNames():
    nDevices = GetGpuCount()
    return [getGpuName(x) for x in range(nDevices)]

        

# deviceProp.GetMultiProcessorCount() * deviceProp.GetSMCores() * deviceProp.GetClockRate()

# multiProcessorCount
