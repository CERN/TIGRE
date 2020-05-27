from libc.stdlib cimport malloc, free
cdef extern from "GpuIds.hpp":
    ctypedef struct GpuIds:
        # char* m_strName;
        int*  m_piDeviceIds; 
        int m_iCount; 

#TODO: Change from inline to stop duplication. Discussed on slack "python-questions" 28-Mar-2017
cdef inline void free_c_gpuids(GpuIds* c_gpuids):
    # free(c_gpuids.m_strName)
    if c_gpuids.m_piDeviceIds != NULL:
        free(c_gpuids.m_piDeviceIds)
    free(c_gpuids)

#TODO: Change from inline to stop duplication. Discussed on slack "python-questions" 28-Mar-2017
cdef inline GpuIds* convert_to_c_gpuids(p_gpuids):
    cdef GpuIds* c_gpuids =<GpuIds *>malloc(sizeof(GpuIds))
    if not c_gpuids:
        MemoryError()
    if p_gpuids is not None:
        c_gpuids.m_iCount = len(p_gpuids)
    else:
        c_gpuids.m_iCount = 0
    
    if c_gpuids.m_iCount > 0:
        c_gpuids.m_piDeviceIds = <int*>malloc(c_gpuids.m_iCount * sizeof(int))
        if not c_gpuids.m_piDeviceIds:
            raise MemoryError()
        for iI in range(c_gpuids.m_iCount):
            c_gpuids.m_piDeviceIds[iI] = p_gpuids.devices[iI]
    else:
        c_gpuids.m_iCount = 0
        c_gpuids.m_piDeviceIds = NULL

    return c_gpuids
