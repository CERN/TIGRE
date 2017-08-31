from _types cimport Geometry as c_Geometry, convert_to_c_geometry, free_c_geometry

cimport numpy as np
import numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

from libc.stdlib cimport malloc, free

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef extern from "voxel_backprojection.hpp":
    cdef int voxel_backprojection(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha)
cdef extern from "voxel_backprojection2.hpp":
    cdef int voxel_backprojection2(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha)
cdef extern from "voxel_backprojection_parallel.hpp":
    cdef int voxel_backprojection_parallel(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha)

def Atb(np.ndarray[np.float32_t, ndim=3] projections, geometry, np.ndarray[np.float32_t, ndim=1] angles, krylov="matched", mode="cone"):
    cdef int total_projections = angles.size
    cdef c_Geometry* c_geometry = convert_to_c_geometry(geometry, total_projections)
    cdef float* c_model = <float*> malloc(geometry.nVoxel[0] * geometry.nVoxel[1] * geometry.nVoxel[2] * sizeof(float))
    cdef float* c_angles = <float*> angles.data

    # TODO: Error if krylov isn't FDK or matched
    if krylov == "matched":
        krylov_proj = True
    elif krylov == "FDK":
        krylov_proj = False
    else:
        print("Error: Unknown krylov, using default matched")
        krylov_proj = True

    if mode == "parallel":
        cone_beam = False
    elif mode == "cone":
        cone_beam = True
    else:
        print("Error: Unknown mode, using default cone beam")
        cone_beam = True

    projections = projections.swapaxes(0,1).copy(order='F')
    cdef float* c_projections = <float*> projections.data

    if cone_beam:
        if krylov_proj:
            voxel_backprojection2(c_projections, c_geometry[0], c_model, c_angles, total_projections)
        else:
            voxel_backprojection(c_projections, c_geometry[0], c_model, c_angles, total_projections)
    else:
        voxel_backprojection_parallel(c_projections, c_geometry[0], c_model, c_angles, total_projections)

    projections = projections.swapaxes(0,1).copy(order='C')

    cdef np.npy_intp shape[3]
    shape[0] = <np.npy_intp> geometry.nVoxel[2]
    shape[1] = <np.npy_intp> geometry.nVoxel[1]
    shape[2] = <np.npy_intp> geometry.nVoxel[0]

    # TODO: Swap axis here could be making a copy
    model = np.PyArray_SimpleNewFromData(3, shape, np.NPY_FLOAT32, c_model).swapaxes(0,2)
    PyArray_ENABLEFLAGS(model, np.NPY_OWNDATA) # Attribute new memory owner

    free_c_geometry(c_geometry)

    return model