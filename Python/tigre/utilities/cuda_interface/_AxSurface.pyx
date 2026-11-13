cimport numpy as np
import numpy as np

np.import_array()

cdef extern from "Siddon_surface_projection.hpp":
    cdef int siddon_surface_projection(const float* img,
            int nx, int ny, int nz,
            float dx, float dy, float dz,
            float offx, float offy, float offz,
            int nviews, int rows, int cols,
            const float* source,
            const float* pixels, int pixel_views,
            const unsigned char* mask, int mask_views,
            float* result, int device) nogil
    cdef int siddon_surface_backprojection(const float* proj,
            int nx, int ny, int nz,
            float dx, float dy, float dz,
            float offx, float offy, float offz,
            int nviews, int rows, int cols,
            const float* source,
            const float* pixels, int pixel_views,
            const unsigned char* mask, int mask_views,
            float* volume, int device) nogil



def _validate_inputs(source, pixels, mask, voxel_size, offset_origin, device):
    if source.shape[0] < 1 or source.shape[1] != 3 or not np.isfinite(source).all():
        raise ValueError("source must have finite shape (views, 3)")
    if (pixels.shape[0] not in (1, len(source)) or pixels.shape[-1] != 3
            or min(pixels.shape[1:3]) < 1 or not np.isfinite(pixels).all()):
        raise ValueError("pixels must have finite shape (one or views, rows, columns, 3)")
    if mask is not None:
        m = np.asarray(mask)
        if (m.ndim != 3 or m.shape[0] not in (1, len(source))
                or m.shape[1:] != pixels.shape[1:3] or np.any((m != 0) & (m != 1))):
            raise ValueError("mask must be binary (one or views, rows, columns)")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        spacing = np.asarray(voxel_size, dtype=np.float32)
        offset = np.asarray(offset_origin, dtype=np.float32)
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("voxel_size must contain three positive finite values")
    if offset.shape != (3,) or not np.isfinite(offset).all():
        raise ValueError("offset_origin must contain three finite values")
    if device < 0:
        raise ValueError("device must be nonnegative")


def _Atb_surface_ext(np.ndarray[np.float32_t, ndim=3] proj,
                     shape_zyx, voxel_size, offset_origin,
                     np.ndarray[np.float32_t, ndim=2] source,
                     np.ndarray[np.float32_t, ndim=4] pixels,
                     mask=None, int device=0):
    """Discrete adjoint of _Ax_surface_ext, up to float32 accumulation error."""
    _validate_inputs(np.asarray(source), np.asarray(pixels), mask, voxel_size, offset_origin, device)
    shape = np.asarray(shape_zyx)
    if shape.shape != (3,) or shape.dtype.kind not in "iu" or np.any(shape < 1) or np.any(shape > np.iinfo(np.int32).max):
        raise ValueError("shape_zyx must contain three positive int32 dimensions")
    if tuple(np.asarray(proj).shape) != (source.shape[0], pixels.shape[1], pixels.shape[2]):
        raise ValueError("projection shape does not match the detector")
    if not np.isfinite(proj).all():
        raise ValueError("projection values must be finite after masking")
    proj = np.ascontiguousarray(proj)
    source = np.ascontiguousarray(source)
    pixels = np.ascontiguousarray(pixels)
    cdef int nviews = source.shape[0]
    cdef int pixel_views = pixels.shape[0]
    cdef int rows = pixels.shape[1]
    cdef int cols = pixels.shape[2]
    cdef int nz = shape_zyx[0]
    cdef int ny = shape_zyx[1]
    cdef int nx = shape_zyx[2]
    cdef float dx = voxel_size[0]
    cdef float dy = voxel_size[1]
    cdef float dz = voxel_size[2]
    cdef float ox = offset_origin[0]
    cdef float oy = offset_origin[1]
    cdef float oz = offset_origin[2]

    cdef np.ndarray[np.uint8_t, ndim=3] mask_arr
    cdef const unsigned char* c_mask = NULL
    cdef int mask_views = 0
    if mask is not None:
        mask_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        mask_views = mask_arr.shape[0]
        c_mask = <const unsigned char*> mask_arr.data

    volume = np.zeros((nz, ny, nx), dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=3] volume_c = volume
    cdef int status
    with nogil:
        status = siddon_surface_backprojection(<const float*> proj.data, nx, ny, nz, dx, dy, dz, ox, oy, oz,
                                               nviews, rows, cols,
                                               <const float*> source.data,
                                               <const float*> pixels.data, pixel_views,
                                               c_mask, mask_views,
                                               <float*> volume_c.data, device)
    if status != 0:
        raise ValueError('TIGRE: Call to Atb_surface failed')
    return volume


def _Ax_surface_ext(np.ndarray[np.float32_t, ndim=3] img,
                    voxel_size, offset_origin,
                    np.ndarray[np.float32_t, ndim=2] source,
                    np.ndarray[np.float32_t, ndim=4] pixels,
                    mask=None, int device=0):
    """Siddon forward projection with explicit source/pixel coordinates.

    img: (nz, ny, nx); source: (views, 3); pixels: (1 or views, rows, cols, 3);
    mask: None or uint8 (1 or views, rows, cols). Returns (views, rows, cols).
    """
    _validate_inputs(np.asarray(source), np.asarray(pixels), mask, voxel_size, offset_origin, device)
    if min(np.asarray(img).shape) < 1 or not np.isfinite(img).all():
        raise ValueError("volume must be nonempty and finite")
    img = np.ascontiguousarray(img)
    source = np.ascontiguousarray(source)
    pixels = np.ascontiguousarray(pixels)
    cdef int nviews = source.shape[0]
    cdef int pixel_views = pixels.shape[0]
    cdef int rows = pixels.shape[1]
    cdef int cols = pixels.shape[2]
    cdef int nz = img.shape[0]
    cdef int ny = img.shape[1]
    cdef int nx = img.shape[2]
    cdef float dx = voxel_size[0]
    cdef float dy = voxel_size[1]
    cdef float dz = voxel_size[2]
    cdef float ox = offset_origin[0]
    cdef float oy = offset_origin[1]
    cdef float oz = offset_origin[2]

    cdef np.ndarray[np.uint8_t, ndim=3] mask_arr
    cdef const unsigned char* c_mask = NULL
    cdef int mask_views = 0
    if mask is not None:
        mask_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        mask_views = mask_arr.shape[0]
        c_mask = <const unsigned char*> mask_arr.data

    result = np.zeros((nviews, rows, cols), dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=3] result_c = result
    cdef int status
    with nogil:
        status = siddon_surface_projection(<const float*> img.data, nx, ny, nz, dx, dy, dz, ox, oy, oz,
                                           nviews, rows, cols,
                                           <const float*> source.data,
                                           <const float*> pixels.data, pixel_views,
                                           c_mask, mask_views,
                                           <float*> result_c.data, device)
    if status != 0:
        raise ValueError('TIGRE: Call to Ax_surface failed')
    return result
