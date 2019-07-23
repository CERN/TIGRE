cimport numpy as np
import numpy as np
from tigre.Source._types cimport Geometry as c_Geometry, convert_to_c_geometry, free_c_geometry

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

from libc.stdlib cimport malloc, free

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
    void PyArray_CLEARFLAGS(np.ndarray arr, int flags)


cdef extern from "Siddon_projection.hpp":
    cdef int siddon_ray_projection(float* img, c_Geometry geo, float** result, float* alphas, int nalpha)
cdef extern from "Siddon_projection_parallel.hpp":
    cdef int siddon_ray_projection_parallel(float* img, c_Geometry geo, float** result, float* alphas, int nalpha)
cdef extern from "ray_interpolated_projection.hpp":
    cdef int interpolation_projection(float* img, c_Geometry geo, float** result, float* alphas, int nalpha)
cdef extern from "ray_interpolated_projection_parallel.hpp":
    cdef int interpolation_projection_parallel(float* img, c_Geometry geo, float** result, float* alphas, int nalpha)

def cuda_raise_errors(error_code):
    if error_code:
        raise ValueError('TIGRE: Call to Ax failed')

def _Ax_ext(np.ndarray[np.float32_t, ndim=3] img, geometry, np.ndarray[np.float32_t, ndim=2] angles, projection_type="ray-voxel", mode="cone"):

    cdef int total_projections = angles.shape[0]

    #PERMUTE INPUT: convert_contig_mode(C_CONTIG) -> F_CONTIG
    # (N, V, U) -> (U, V, N) 
    geometry.convert_contig_mode()
    cdef c_Geometry* c_geometry = convert_to_c_geometry(geometry, total_projections)
    if not c_geometry:
        raise MemoryError()
    cdef float** c_projections = <float**> malloc(total_projections * sizeof(float*))
    if not c_projections:
        raise MemoryError()
    cdef int i = 0
    for i in range(total_projections):
       c_projections[i] = <float*> malloc(geometry.nDetector[0] * geometry.nDetector[1] * sizeof(float))

    cdef float* c_angles = <float*> angles.data
    if not c_angles:
        raise MemoryError()
    if projection_type == "ray-voxel":
        interpolated = False
    elif projection_type == "interpolated":
        interpolated = True
    else:
        print("Warning: Unknown projection_type, using default ray-voxel")
        interpolated = False

    if mode == "parallel":
        cone_beam = False
    elif mode == "cone":
        cone_beam = True
    else:
        print("Warning: Unknown mode, using default cone beam")
        cone_beam = True
    #PERMUTE INPUT: (Z, Y, X) -> (X ,Y ,Z)
    img = img.transpose().copy(order='F')

    cdef float* c_img = <float*> img.data
    if cone_beam:
        if not interpolated:
            cuda_raise_errors(siddon_ray_projection(c_img, c_geometry[0], c_projections, c_angles, total_projections))
        else:
            cuda_raise_errors(interpolation_projection(c_img, c_geometry[0], c_projections, c_angles, total_projections))
    else:
        if not interpolated:
            cuda_raise_errors(siddon_ray_projection_parallel(c_img, c_geometry[0], c_projections, c_angles, total_projections))
        else:
            cuda_raise_errors(interpolation_projection_parallel(c_img, c_geometry[0], c_projections, c_angles, total_projections))
    img = img.copy(order='C')

    cdef np.npy_intp shape[2]
    shape[0] = <np.npy_intp> geometry.nDetector[0]
    shape[1] = <np.npy_intp> geometry.nDetector[1]
    projections = [None] * total_projections
    for i in range(total_projections):
        projections[i] = np.PyArray_SimpleNewFromData(2, shape, np.NPY_FLOAT32, c_projections[i])
        PyArray_ENABLEFLAGS(projections[i], np.NPY_OWNDATA)  # Attribute new memory owner


    free(c_projections)  # Free pointer array, not actual data
    free_c_geometry(c_geometry)
    geometry.convert_contig_mode()

    # PERMUTE OUTPUT: (U, V, N) -> (N, V, U)
    return np.stack(projections,0).swapaxes(1,2).copy(order='C')