#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Siddon_surface_projection.hpp"

#define surfCall(call, msg) \
    do { \
        cudaError_t __err = (call); \
        if (__err != cudaSuccess) { \
            printf("TIGRE surface projection: %s (%s)\n", msg, cudaGetErrorString(__err)); \
            cudaGetLastError(); \
            return 1; \
        } \
    } while (0)

struct SurfaceParams {
    int nx, ny, nz;
    float dx, dy, dz;
    float offx, offy, offz;
    int rows, cols;
};

// A parallel ray on the lower face belongs to the first voxel; the upper
// face is excluded. Handle it explicitly instead of relying on 0/0 in slabs.
__device__ __forceinline__ bool surfaceAxisInterval(float source, float ray,
        int count, float& enter, float& leave)
{
    if (ray == 0.0f) {
        enter = -INFINITY;
        leave = INFINITY;
        return source >= 0.0f && source < count;
    }
    float a = __fdividef(-source, ray);
    float b = __fdividef(count - source, ray);
    enter = fminf(a, b);
    leave = fmaxf(a, b);
    return true;
}

// Releases every device resource it holds, so early error returns do not leak.
struct DeviceResources {
    cudaArray* arr = NULL;
    cudaTextureObject_t tex = 0;
    float* ptrs[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    unsigned char* mask = NULL;
    ~DeviceResources() {
        if (tex) cudaDestroyTextureObject(tex);
        if (arr) cudaFreeArray(arr);
        for (float* q : ptrs) if (q) cudaFree(q);
        if (mask) cudaFree(mask);
        cudaGetLastError();
    }
};

// Siddon's ray-voxel traversal in voxel-index space (voxel i spans [i,i+1)).
// BACK=false: detector is the output, the volume is read through tex.
// BACK=true : detector is the input, its values are scattered into volume
//             with the same Siddon weights (exact adjoint of the forward).
__device__ __forceinline__ void scatterVoxel(float* volume, const SurfaceParams& p,
        float fi, float fj, float fk, float value)
{
    int i = (int)floorf(fi), j = (int)floorf(fj), k = (int)floorf(fk);
    if (i < 0 || j < 0 || k < 0 || i >= p.nx || j >= p.ny || k >= p.nz)
        return;
    atomicAdd(&volume[((size_t)k * p.ny + j) * p.nx + i], value);
}

template <bool BACK>
__global__ void kernelSurface(SurfaceParams p,
        const float* source,       // chunkViews*3 (already in the chunk)
        const float* pixels,       // pixel_stride*chunkViews*rows*cols*3
        const unsigned char* mask, // NULL or mask_stride*chunkViews*rows*cols
        size_t pixelViewStride,    // floats between views in pixels (0 = shared)
        size_t maskViewStride,
        float* detector,
        cudaTextureObject_t tex,
        float* volume)
{
    unsigned int u = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int v = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int view = blockIdx.z;
    if (u >= (unsigned)p.cols || v >= (unsigned)p.rows)
        return;

    size_t pix = (size_t)v * p.cols + u;
    size_t idx = (size_t)view * p.rows * p.cols + pix;

    if (mask != NULL && mask[view * maskViewStride + pix] == 0) {
        if (!BACK) detector[idx] = 0;
        return;
    }

    const float* sw = source + (size_t)view * 3;
    const float* pw = pixels + view * pixelViewStride + pix * 3;

    // world (mm) -> voxel-index coordinates
    float3 source_;
    source_.x = (sw[0] - p.offx + p.nx * p.dx * 0.5f) / p.dx;
    source_.y = (sw[1] - p.offy + p.ny * p.dy * 0.5f) / p.dy;
    source_.z = (sw[2] - p.offz + p.nz * p.dz * 0.5f) / p.dz;
    float3 pixel1D;
    pixel1D.x = (pw[0] - p.offx + p.nx * p.dx * 0.5f) / p.dx;
    pixel1D.y = (pw[1] - p.offy + p.ny * p.dy * 0.5f) / p.dy;
    pixel1D.z = (pw[2] - p.offz + p.nz * p.dz * 0.5f) / p.dz;

    float3 ray;
    ray.x = pixel1D.x - source_.x;
    ray.y = pixel1D.y - source_.y;
    ray.z = pixel1D.z - source_.z;
    float eps = 0.001f;
    ray.x = (fabsf(ray.x) < eps) ? 0 : ray.x;
    ray.y = (fabsf(ray.y) < eps) ? 0 : ray.y;
    ray.z = (fabsf(ray.z) < eps) ? 0 : ray.z;

    float axm, aym, azm, axM, ayM, azM;
    if ((ray.x == 0.0f && ray.y == 0.0f && ray.z == 0.0f) ||
            !surfaceAxisInterval(source_.x, ray.x, p.nx, axm, axM) ||
            !surfaceAxisInterval(source_.y, ray.y, p.ny, aym, ayM) ||
            !surfaceAxisInterval(source_.z, ray.z, p.nz, azm, azM)) {
        if (!BACK) detector[idx] = 0;
        return;
    }

    // The ray is the segment source -> pixel (alpha in [0, 1]); clip so that a pixel or
    // source lying inside the volume box does not extend the path past the segment.
    float am = fmaxf(fmaxf(fmaxf(axm, aym), azm), 0.0f);
    float aM = fminf(fminf(axM, ayM), azM);
    bool clipEnd = aM > 1.0f;
    if (clipEnd) aM = 1.0f;

    // Ray misses the volume (also guards the NaN/inf case of axis-aligned rays).
    if (!(am < aM)) {
        if (!BACK) detector[idx] = 0;
        return;
    }

    float imin, imax, jmin, jmax, kmin, kmax;
    if (source_.x < pixel1D.x) {
        imin = (am == axm) ? 1.0f : ceilf(source_.x + am * ray.x);
        imax = (aM == axM) ? p.nx : floorf(source_.x + aM * ray.x);
    } else {
        imax = (am == axm) ? p.nx - 1.0f : floorf(source_.x + am * ray.x);
        imin = (aM == axM) ? 0.0f : ceilf(source_.x + aM * ray.x);
    }
    if (source_.y < pixel1D.y) {
        jmin = (am == aym) ? 1.0f : ceilf(source_.y + am * ray.y);
        jmax = (aM == ayM) ? p.ny : floorf(source_.y + aM * ray.y);
    } else {
        jmax = (am == aym) ? p.ny - 1.0f : floorf(source_.y + am * ray.y);
        jmin = (aM == ayM) ? 0.0f : ceilf(source_.y + aM * ray.y);
    }
    if (source_.z < pixel1D.z) {
        kmin = (am == azm) ? 1.0f : ceilf(source_.z + am * ray.z);
        kmax = (aM == azM) ? p.nz : floorf(source_.z + aM * ray.z);
    } else {
        kmax = (am == azm) ? p.nz - 1.0f : floorf(source_.z + am * ray.z);
        kmin = (aM == azM) ? 0.0f : ceilf(source_.z + aM * ray.z);
    }

    float ax, ay, az;
    ax = (source_.x < pixel1D.x) ? __fdividef(imin - source_.x, ray.x) : __fdividef(imax - source_.x, ray.x);
    ay = (source_.y < pixel1D.y) ? __fdividef(jmin - source_.y, ray.y) : __fdividef(jmax - source_.y, ray.y);
    az = (source_.z < pixel1D.z) ? __fdividef(kmin - source_.z, ray.z) : __fdividef(kmax - source_.z, ray.z);
    ax = isinf(ax) ? fabsf(ax) : ax;
    ay = isinf(ay) ? fabsf(ay) : ay;
    az = isinf(az) ? fabsf(az) : az;
    if (ray.x == 0.0f) ax = INFINITY;
    if (ray.y == 0.0f) ay = INFINITY;
    if (ray.z == 0.0f) az = INFINITY;

    float aminc = fminf(fminf(ax, ay), az);
    unsigned long i, j, k;
    i = (unsigned long)floorf(source_.x + (aminc + am) * 0.5f * ray.x);
    j = (unsigned long)floorf(source_.y + (aminc + am) * 0.5f * ray.y);
    k = (unsigned long)floorf(source_.z + (aminc + am) * 0.5f * ray.z);
    float ac = am;
    float axu = __frcp_rd(fabsf(ray.x));
    float ayu = __frcp_rd(fabsf(ray.y));
    float azu = __frcp_rd(fabsf(ray.z));
    float iu = (source_.x < pixel1D.x) ? 1.0f : -1.0f;
    float ju = (source_.y < pixel1D.y) ? 1.0f : -1.0f;
    float ku = (source_.z < pixel1D.z) ? 1.0f : -1.0f;

    float maxlength = __fsqrt_rd(ray.x * ray.x * p.dx * p.dx + ray.y * ray.y * p.dy * p.dy + ray.z * ray.z * p.dz * p.dz);
    float sum = 0.0f;
    float proj = BACK ? detector[idx] * maxlength : 0.0f;
    unsigned long Np = (imax - imin + 1) + (jmax - jmin + 1) + (kmax - kmin + 1);
    float fi = i + 0.5f, fj = j + 0.5f, fk = k + 0.5f;
    for (unsigned long ii = 0; ii < Np; ii++) {
        if (clipEnd && aminc >= aM) break;
        if (ax == aminc) {
            if (BACK) scatterVoxel(volume, p, fi, fj, fk, (ax - ac) * proj);
            else sum += (ax - ac) * tex3D<float>(tex, fi, fj, fk);
            fi += iu; ac = ax; ax += axu;
        } else if (ay == aminc) {
            if (BACK) scatterVoxel(volume, p, fi, fj, fk, (ay - ac) * proj);
            else sum += (ay - ac) * tex3D<float>(tex, fi, fj, fk);
            fj += ju; ac = ay; ay += ayu;
        } else if (az == aminc) {
            if (BACK) scatterVoxel(volume, p, fi, fj, fk, (az - ac) * proj);
            else sum += (az - ac) * tex3D<float>(tex, fi, fj, fk);
            fk += ku; ac = az; az += azu;
        }
        aminc = fminf(fminf(ax, ay), az);
    }
    if (clipEnd && aM > ac) {
        // last partial voxel, ending at the pixel inside the volume
        if (BACK) scatterVoxel(volume, p, fi, fj, fk, (aM - ac) * proj);
        else sum += (aM - ac) * tex3D<float>(tex, fi, fj, fk);
    }
    if (!BACK) detector[idx] = sum * maxlength;
}

static int surfaceDriver(bool back,
        const float* img, float* volumeOut, const float* projIn,
        int nx, int ny, int nz,
        float dx, float dy, float dz,
        float offx, float offy, float offz,
        int nviews, int rows, int cols,
        const float* source,
        const float* pixels, int pixel_views,
        const unsigned char* mask, int mask_views,
        float* result, int device)
{
    if (nviews <= 0 || rows <= 0 || cols <= 0)
        return 1;
    if (pixel_views != 1 && pixel_views != nviews)
        return 1;
    if (mask != NULL && mask_views != 1 && mask_views != nviews)
        return 1;

    surfCall(cudaSetDevice(device), "set device");

    cudaExtent extent = make_cudaExtent(nx, ny, nz);
    DeviceResources res;
    cudaArray*& arr = res.arr;
    cudaTextureObject_t& tex = res.tex;
    float*& d_vol = res.ptrs[0];
    float*& d_source = res.ptrs[1];
    float*& d_pixels = res.ptrs[2];
    float*& d_det = res.ptrs[3];
    unsigned char*& d_mask = res.mask;
    if (!back) {
        // Volume as a texture (the whole volume must fit on the device).
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        surfCall(cudaMalloc3DArray(&arr, &channelDesc, extent), "volume allocation (volume must fit in GPU memory)");
        cudaMemcpy3DParms copyParams = {0};
        copyParams.srcPtr = make_cudaPitchedPtr((void*)img, extent.width * sizeof(float), extent.width, extent.height);
        copyParams.dstArray = arr;
        copyParams.extent = extent;
        copyParams.kind = cudaMemcpyHostToDevice;
        surfCall(cudaMemcpy3D(&copyParams), "volume copy");

        cudaResourceDesc texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array = arr;
        cudaTextureDesc texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModePoint;
        texDescr.addressMode[0] = cudaAddressModeBorder;
        texDescr.addressMode[1] = cudaAddressModeBorder;
        texDescr.addressMode[2] = cudaAddressModeBorder;
        texDescr.readMode = cudaReadModeElementType;
        surfCall(cudaCreateTextureObject(&tex, &texRes, &texDescr, NULL), "texture creation");
    } else {
        surfCall(cudaMalloc((void**)&d_vol, sizeof(float) * (size_t)nx * ny * nz), "volume allocation (volume must fit in GPU memory)");
        surfCall(cudaMemset(d_vol, 0, sizeof(float) * (size_t)nx * ny * nz), "volume initialization");
    }

    size_t npix = (size_t)rows * cols;
    // Bound the per-chunk device footprint to ~256 MB of pixel table.
    size_t maxFloats = (size_t)64 * 1024 * 1024;
    int chunk = (pixel_views == 1) ? nviews : (int)std::max<size_t>(1, maxFloats / (npix * 3));
    chunk = std::min(chunk, nviews);
    chunk = std::min(chunk, 65535);

    surfCall(cudaMalloc((void**)&d_source, sizeof(float) * 3 * chunk), "source allocation");
    surfCall(cudaMalloc((void**)&d_det, sizeof(float) * npix * chunk), "detector allocation");
    size_t pixelFloats = (size_t)(pixel_views == 1 ? 1 : chunk) * npix * 3;
    surfCall(cudaMalloc((void**)&d_pixels, sizeof(float) * pixelFloats), "pixel allocation");
    if (mask != NULL)
        surfCall(cudaMalloc((void**)&d_mask, (size_t)(mask_views == 1 ? 1 : chunk) * npix), "mask allocation");

    if (pixel_views == 1)
        surfCall(cudaMemcpy(d_pixels, pixels, sizeof(float) * npix * 3, cudaMemcpyHostToDevice), "shared pixel copy");
    if (mask != NULL && mask_views == 1)
        surfCall(cudaMemcpy(d_mask, mask, npix, cudaMemcpyHostToDevice), "shared mask copy");

    SurfaceParams p = {nx, ny, nz, dx, dy, dz, offx, offy, offz, rows, cols};
    dim3 block(16, 16, 1);
    for (int start = 0; start < nviews; start += chunk) {
        int n = std::min(chunk, nviews - start);
        surfCall(cudaMemcpy(d_source, source + (size_t)start * 3, sizeof(float) * 3 * n, cudaMemcpyHostToDevice), "source copy");
        if (pixel_views != 1)
            surfCall(cudaMemcpy(d_pixels, pixels + (size_t)start * npix * 3, sizeof(float) * npix * 3 * n, cudaMemcpyHostToDevice), "pixel copy");
        if (mask != NULL && mask_views != 1)
            surfCall(cudaMemcpy(d_mask, mask + (size_t)start * npix, npix * n, cudaMemcpyHostToDevice), "mask copy");
        dim3 grid((cols + block.x - 1) / block.x, (rows + block.y - 1) / block.y, n);
        if (back) {
            surfCall(cudaMemcpy(d_det, projIn + (size_t)start * npix, sizeof(float) * npix * n, cudaMemcpyHostToDevice), "projection copy");
            kernelSurface<true><<<grid, block>>>(p, d_source, d_pixels, d_mask,
                    pixel_views == 1 ? 0 : npix * 3, mask_views == 1 ? 0 : npix, d_det, tex, d_vol);
        } else {
            kernelSurface<false><<<grid, block>>>(p, d_source, d_pixels, d_mask,
                    pixel_views == 1 ? 0 : npix * 3, mask_views == 1 ? 0 : npix, d_det, tex, NULL);
        }
        surfCall(cudaGetLastError(), "kernel launch");
        if (!back)
            surfCall(cudaMemcpy(result + (size_t)start * npix, d_det, sizeof(float) * npix * n, cudaMemcpyDeviceToHost), "forward result copy");
    }
    if (back) {
        surfCall(cudaMemcpy(volumeOut, d_vol, sizeof(float) * (size_t)nx * ny * nz, cudaMemcpyDeviceToHost), "adjoint result copy");
    }
    return 0;
}

int siddon_surface_projection(const float* img,
        int nx, int ny, int nz,
        float dx, float dy, float dz,
        float offx, float offy, float offz,
        int nviews, int rows, int cols,
        const float* source,
        const float* pixels, int pixel_views,
        const unsigned char* mask, int mask_views,
        float* result, int device)
{
    return surfaceDriver(false, img, NULL, NULL, nx, ny, nz, dx, dy, dz, offx, offy, offz,
            nviews, rows, cols, source, pixels, pixel_views, mask, mask_views, result, device);
}

int siddon_surface_backprojection(const float* proj,
        int nx, int ny, int nz,
        float dx, float dy, float dz,
        float offx, float offy, float offz,
        int nviews, int rows, int cols,
        const float* source,
        const float* pixels, int pixel_views,
        const unsigned char* mask, int mask_views,
        float* volume, int device)
{
    return surfaceDriver(true, NULL, volume, proj, nx, ny, nz, dx, dy, dz, offx, offy, offz,
            nviews, rows, cols, source, pixels, pixel_views, mask, mask_views, NULL, device);
}
