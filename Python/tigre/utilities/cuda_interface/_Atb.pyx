
from tigre.utilities.cuda_interface._types cimport Geometry as c_Geometry, convert_to_c_geometry, free_c_geometry
from tigre.utilities.cuda_interface._gpuUtils cimport GpuIds as c_GpuIds, convert_to_c_gpuids, free_c_gpuids

cimport numpy as np
import numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

from libc.stdlib cimport malloc, free
from tigre.utilities.errors import TigreCudaCallError
cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef extern from "voxel_backprojection.hpp":
    cdef int voxel_backprojection(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha, const c_GpuIds& gpuids)
cdef extern from "voxel_backprojection2.hpp":
    cdef int voxel_backprojection2(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha, const c_GpuIds& gpuids)
cdef extern from "voxel_backprojection_parallel.hpp":
    cdef int voxel_backprojection_parallel(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha, const c_GpuIds& gpuids)


def cuda_raise_errors(error_code):
    if error_code:
        raise TigreCudaCallError('Atb:',error_code)



def _Atb_ext(np.ndarray[np.float32_t, ndim=3] projections, geometry, np.ndarray[np.float32_t, ndim=2] angles, backprojection_type="FDK", mode="cone", gpuids=None):

    cdef c_GpuIds* c_gpuids = convert_to_c_gpuids(gpuids)
    if not c_gpuids:
        raise MemoryError()
    
    cdef int total_projections = angles.shape[0]

    cdef c_Geometry* c_geometry = convert_to_c_geometry(geometry, total_projections)

    cdef float* c_model = <float*> malloc(geometry.nVoxel[0] * geometry.nVoxel[1] * geometry.nVoxel[2] * sizeof(float))
    cdef float* c_angles = <float*> angles.data

    # TODO: Error if backprojection_type isn't FDK or matched
    if backprojection_type == "matched":
        krylov_proj = True
    elif backprojection_type == "FDK":
        krylov_proj = False
    else:
        print("Warning: Unknown backprojector, using default matched")
        krylov_proj = True

    if mode == "parallel":
        cone_beam = False
    elif mode == "cone":
        cone_beam = True
    else:
        print("Warning: Unknown mode, using default cone beam")
        cone_beam = True

    cdef float* c_projections = <float*> projections.data

    if cone_beam:
        if krylov_proj:
            cuda_raise_errors(voxel_backprojection2(c_projections, c_geometry[0], c_model, c_angles, total_projections, c_gpuids[0]))
        else:
            cuda_raise_errors(voxel_backprojection(c_projections, c_geometry[0], c_model, c_angles, total_projections, c_gpuids[0]))


    else:
        cuda_raise_errors(voxel_backprojection_parallel(c_projections, c_geometry[0], c_model, c_angles, total_projections, c_gpuids[0]))

    cdef np.npy_intp shape[3]
    shape[0] = <np.npy_intp> geometry.nVoxel[0]
    shape[1] = <np.npy_intp> geometry.nVoxel[1]
    shape[2] = <np.npy_intp> geometry.nVoxel[2]

    # TODO: Swap axis here could be making a copy
    model = np.PyArray_SimpleNewFromData(3, shape, np.NPY_FLOAT32, c_model)
    PyArray_ENABLEFLAGS(model, np.NPY_OWNDATA) # Attribute new memory owner

    free_c_geometry(c_geometry)


    return model
