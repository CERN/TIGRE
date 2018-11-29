from tigre.Source._types cimport Geometry as c_Geometry, convert_to_c_geometry, free_c_geometry

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
cdef extern from "voxel_backprojection_spherical.hpp":
    cdef int voxel_backprojection_spherical(float* projections, c_Geometry geo, float* result,float * angles,int nalpha)
cdef extern from "voxel_backprojection2.hpp":
    cdef int voxel_backprojection2(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha)
cdef extern from "voxel_backprojection2_spherical.hpp":
    cdef int voxel_backprojection2_spherical(float* projections, c_Geometry geo, float* result,float * angles,int nalpha)
cdef extern from "voxel_backprojection_parallel.hpp":
    cdef int voxel_backprojection_parallel(float* projections, c_Geometry geo, float* result,float * alphas,int nalpha)
cdef extern from "voxel_backprojection_parallel_spherical.hpp":
    cdef int voxel_backprojection_parallel_spherical(float* projections, c_Geometry geo, float* result,float * angles,int nalpha)

def _Atb_ext(np.ndarray[np.float32_t, ndim=3] projections, geometry, np.ndarray[np.float32_t, ndim=2] angles, krylov="matched", mode="cone"):
    cdef int total_projections = angles.shape[0]
    geometry.convert_contig_mode()
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

    projections = projections.swapaxes(1,2).swapaxes(0,2).copy(order='F')
    cdef float* c_projections = <float*> projections.data

    cdef float theta,psi;
    theta=0;
    psi=0;
    for i in range(total_projections):
        theta+=abs(c_angles[i*3+1])
        psi  +=abs(c_angles[i*3+2])
    
    standard_rotation = True;
    
    if psi == 0.0 and theta == 0.0:
        standard_rotation=True
    else:
        standard_rotation=False

    if cone_beam:
        if krylov_proj:
            if standard_rotation:
                #print("vox_back2 being called")
                voxel_backprojection2(c_projections, c_geometry[0], c_model, c_angles, total_projections)
            else:
                #print("vox_back2_spherical being called")
                voxel_backprojection2_spherical(c_projections, c_geometry[0], c_model, c_angles, total_projections)
        else:
            if standard_rotation:
                #print("vox_back being called")
                voxel_backprojection(c_projections, c_geometry[0], c_model, c_angles, total_projections)
            else:
                #print("vox_back_spherical being called")
                voxel_backprojection_spherical(c_projections, c_geometry[0], c_model, c_angles, total_projections)


    else:
        if standard_rotation:
            #print("voxel_backprojection_parallel being called")
            voxel_backprojection_parallel(c_projections, c_geometry[0], c_model, c_angles, total_projections)
        else:
            #print("voxel_backprojection_parallel_spherical being called")
            voxel_backprojection_parallel_spherical(c_projections, c_geometry[0], c_model, c_angles, total_projections)
    
    projections = projections.swapaxes(0,2).swapaxes(1,2).copy(order='C')

    cdef np.npy_intp shape[3]
    shape[0] = <np.npy_intp> geometry.nVoxel[2]
    shape[1] = <np.npy_intp> geometry.nVoxel[1]
    shape[2] = <np.npy_intp> geometry.nVoxel[0]

    # TODO: Swap axis here could be making a copy
    model = np.PyArray_SimpleNewFromData(3, shape, np.NPY_FLOAT32, c_model)
    PyArray_ENABLEFLAGS(model, np.NPY_OWNDATA) # Attribute new memory owner

    free_c_geometry(c_geometry)
    geometry.convert_contig_mode()
    return model