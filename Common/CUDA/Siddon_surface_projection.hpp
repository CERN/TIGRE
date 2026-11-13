#ifndef SIDDON_SURFACE_PROJECTION_HPP
#define SIDDON_SURFACE_PROJECTION_HPP

/* Ray-voxel (Siddon) forward projection for an arbitrary detector surface.
 *
 * All positions are in mm, in the volume-fixed lab frame (x,y,z); the volume is
 * centred at offOrig (usually 0). The image is laid out with x fastest
 * (numpy shape (nz,ny,nx)).
 *
 * source : nviews*3 floats (x,y,z)
 * pixels : pixel_views*rows*cols*3 floats; pixel_views is 1 (shared by all views)
 *          or nviews
 * mask   : NULL or mask_views*rows*cols bytes (0 = missing sample), mask_views 1 or nviews
 * result : nviews*rows*cols floats; masked samples are set to 0
 *
 * Returns 0 on success.
 */
int siddon_surface_projection(const float* img,
        int nx, int ny, int nz,
        float dx, float dy, float dz,
        float offx, float offy, float offz,
        int nviews, int rows, int cols,
        const float* source,
        const float* pixels, int pixel_views,
        const unsigned char* mask, int mask_views,
        float* result, int device);

/* Exact adjoint of siddon_surface_projection: scatters proj (nviews*rows*cols)
 * into volume (nx*ny*nz, overwritten) with the same ray weights. */
int siddon_surface_backprojection(const float* proj,
        int nx, int ny, int nz,
        float dx, float dy, float dz,
        float offx, float offy, float offz,
        int nviews, int rows, int cols,
        const float* source,
        const float* pixels, int pixel_views,
        const unsigned char* mask, int mask_views,
        float* volume, int device);


#endif
