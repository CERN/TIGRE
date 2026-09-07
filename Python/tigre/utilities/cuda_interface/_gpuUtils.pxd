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
        raise MemoryError("Error allocating memory for GPU IDs")
    if p_gpuids is None:
        # Default to EVERY visible GPU, not to no GPU.
        #
        # This used to build a valid pointer with m_iCount = 0 and an empty
        # device list, and hand that to CUDA. The algorithm classes hide it -
        # they all do `if self.gpuids is None: self.gpuids = GpuIds()` first -
        # so tigre.Ax/Atb and every algs.* entry point look fine. But the
        # kernels reachable directly, minTV and AwminTV among them, take
        # gpuids=None as their own default signature and went straight through
        # here: minTV(img, alpha, iters, None) ran with zero devices and killed
        # the interpreter outright, no exception, no message.
        #
        # GpuIds() with no name means "all of them", which is the only sane
        # reading of an unspecified device, and matches what the algorithm
        # classes were already doing by hand.
        from tigre.utilities.gpu import GpuIds as _AllGpuIds
        p_gpuids = _AllGpuIds()
    c_gpuids.m_iCount = len(p_gpuids)
    
    if c_gpuids.m_iCount > 0:
        c_gpuids.m_piDeviceIds = <int*>malloc(c_gpuids.m_iCount * sizeof(int))
        if not c_gpuids.m_piDeviceIds:
            raise MemoryError("Error allocating memory for GPU IDs")
        for iI in range(c_gpuids.m_iCount):
            c_gpuids.m_piDeviceIds[iI] = p_gpuids.devices[iI]
    else:
        c_gpuids.m_iCount = 0
        c_gpuids.m_piDeviceIds = NULL

    return c_gpuids
