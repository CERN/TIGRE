
#include "cuda_to_hip.h"
#include "gpuUtils.hpp"
#include <string.h>
#include <stdio.h>

#if defined(USE_HIP) && (defined(__CUDACC__) || defined(__HIPCC__))
/*
 * Verified self-test for fp32 element-read hardware linear texture filtering.
 *
 * Some ROCm/hardware combinations reject creation of a cudaFilterModeLinear
 * texture over an element-read float array outright; others accept it but
 * silently point-sample, which would ship wrong results. So we do not trust
 * creation success alone: we build a tiny 4x1x1 ramp texture (texel x has
 * value x), create it Linear, then sample at x=2.0 (texel centre 1.5..2.5).
 * Hardware linear returns 1.5; point sampling returns 2.0. The result is
 * cached for the process and the verdict is logged once to stderr.
 */
__global__ static void tigre_hw_linear_probe(cudaTextureObject_t tex, float* out) {
    /* x=2.0 is the midpoint between texel centres 1.5 (value 1) and 2.5 (value
     * 2): hardware linear returns 1.5, point sampling returns 2.0. y,z sit on
     * interior texel centres so border addressing does not pollute the lerp. */
    *out = tex3D<float>(tex, 2.0f, 2.5f, 2.5f);
}

static bool tigre_run_hw_linear_test() {
    /* Mirror the real projection/backprojection textures: a true 3D array with
     * cudaAddressModeBorder and cudaReadModeElementType, created Linear. Some
     * ROCm/hardware combinations accept a 1D Linear texture but reject this 3D
     * border-addressed one, so the probe must match the real configuration. */
    const int nx = 4, ny = 4, nz = 4;
    float ramp[nz * ny * nx];
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                ramp[(z * ny + y) * nx + x] = (float)x;

    cudaChannelFormatDesc ch = cudaCreateChannelDesc<float>();
    cudaArray_t arr = 0;
    cudaExtent ext = make_cudaExtent(nx, ny, nz);
    if (cudaMalloc3DArray(&arr, &ch, ext) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }

    cudaMemcpy3DParms cp = {0};
    cp.srcPtr = make_cudaPitchedPtr((void*)ramp, nx * sizeof(float), nx, ny);
    cp.dstArray = arr;
    cp.extent = ext;
    cp.kind = cudaMemcpyHostToDevice;
    if (cudaMemcpy3D(&cp) != cudaSuccess) {
        cudaGetLastError();
        cudaFreeArray(arr);
        return false;
    }

    cudaResourceDesc res;
    memset(&res, 0, sizeof(res));
    res.resType = cudaResourceTypeArray;
    res.res.array.array = arr;

    cudaTextureDesc td;
    memset(&td, 0, sizeof(td));
    td.normalizedCoords = false;
    td.filterMode = cudaFilterModeLinear;
    td.addressMode[0] = cudaAddressModeBorder;
    td.addressMode[1] = cudaAddressModeBorder;
    td.addressMode[2] = cudaAddressModeBorder;
    td.readMode = cudaReadModeElementType;

    cudaTextureObject_t tex = 0;
    if (cudaCreateTextureObject(&tex, &res, &td, NULL) != cudaSuccess) {
        cudaGetLastError();
        cudaFreeArray(arr);
        return false;
    }

    float* dout = 0;
    if (cudaMalloc(&dout, sizeof(float)) != cudaSuccess) {
        cudaGetLastError();
        cudaDestroyTextureObject(tex);
        cudaFreeArray(arr);
        return false;
    }

    tigre_hw_linear_probe<<<1, 1>>>(tex, dout);
    bool ok = (cudaDeviceSynchronize() == cudaSuccess);
    float v = 0.0f;
    if (ok && cudaMemcpy(&v, dout, sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
        ok = false;
    }
    cudaGetLastError();

    cudaFree(dout);
    cudaDestroyTextureObject(tex);
    cudaFreeArray(arr);

    return ok && (fabsf(v - 1.5f) < 0.05f);
}

bool tigre_hw_linear_supported() {
    static int cached = -1;
    if (cached < 0) {
        bool ok = tigre_run_hw_linear_test();
        char name[128];
        int dev = 0;
        cudaGetDevice(&dev);
        GetGpuName(dev, name);
        fprintf(stderr,
                "TIGRE: fp32 hardware linear texture filtering on %s: %s\n",
                name,
                ok ? "supported (hardware trilinear)"
                   : "unsupported (software trilinear fallback)");
        cached = ok ? 1 : 0;
    }
    return cached != 0;
}
#endif

int GetGpuIdArray(const char* kacGPUName, int* piDeviceIds, int iIdCountMax, char* pcMessage) {
    if (pcMessage) {
        for (int iI = 0; iI < 65535; ++iI) {
            pcMessage[iI] = '\0';
        }
    }
    if (piDeviceIds == 0 || iIdCountMax == 0) {
        return 0;
    }
    int iMessagePos = 0;
    // Count installed GPUs.
    int iCudaDeviceCount = GetGpuCount();
    iMessagePos += sprintf(pcMessage + iMessagePos, "Found GPUs: %d\n", iCudaDeviceCount);
    if (iCudaDeviceCount == 0) {
        // printf("No GPU found\n");
        return 0;
    }

    iCudaDeviceCount = min(iCudaDeviceCount, iIdCountMax);
    iMessagePos += sprintf(pcMessage + iMessagePos, "Max GPUs: %d\n", iCudaDeviceCount);
    if (strlen(kacGPUName) == 0) {
        // Semi-compatible mode:
        //    Return all GPUs
        for (int iI = 0; iI < iCudaDeviceCount; ++iI) {
            piDeviceIds[iI] = iI;
        }
        return iCudaDeviceCount;
    }

    cudaError_t err;
    cudaDeviceProp propDevice;
    int nMatch = 0;
    for (int iId = 0; iId < iCudaDeviceCount; ++iId) {
        err = cudaGetDeviceProperties(&propDevice, iId);
        iMessagePos += sprintf(pcMessage + iMessagePos, "propDevice.name = %s\n", propDevice.name);
        if (strcmp(propDevice.name, kacGPUName) == 0) {
            piDeviceIds[nMatch] = iId;
            ++nMatch;
        }
    }

    for (int iI = 0; iI < nMatch; ++iI) {
        iMessagePos += sprintf(pcMessage + iMessagePos, "%d, ", piDeviceIds[iI]);
    }  
    return nMatch;

}

void GetGpuName(int iDeviceId, char* pcName) {
    memset(pcName, 0, 128);
    cudaError_t err;
    cudaDeviceProp propDevice;
    int id = iDeviceId;
    err = cudaGetDeviceProperties(&propDevice, id);
    memcpy(pcName, propDevice.name, strlen(propDevice.name)*sizeof(char));
}


int GetGpuCount() {
    int iCudaDeviceCount = 0;
    cudaGetDeviceCount(&iCudaDeviceCount);
    return iCudaDeviceCount;
}
