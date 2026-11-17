
#include "gpuUtils.hpp"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <string.h>
#include <stdio.h>

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


bool UsePeerHaloExchange(int iDevA, int iDevB, bool bManySmallCopies) {
    if (iDevA == iDevB) {
        return false;
    }
#if defined(__HIP_PLATFORM_AMD__)
    // Measured on MI355X: peer copies win for a single array per boundary (GD_TV, GD_AwTV), but several
    // concurrent small peer copies (tv_proximal) serialize on ROCm and end up slower than staging.
    if (bManySmallCopies) {
        return false;
    }
#else
    // Measured on V100 with NVLink only. On PCIe, peer access can be reported while the copies are slow,
    // so require a direct link. Native peer atomics are supported over NVLink but not over PCIe, which
    // makes them a usable proxy for "direct link" without linking against NVML.
    int iAtomicAB = 0, iAtomicBA = 0;
    if (cudaDeviceGetP2PAttribute(&iAtomicAB, cudaDevP2PAttrNativeAtomicSupported, iDevA, iDevB) != cudaSuccess ||
        cudaDeviceGetP2PAttribute(&iAtomicBA, cudaDevP2PAttrNativeAtomicSupported, iDevB, iDevA) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    if (!iAtomicAB || !iAtomicBA) {
        return false;
    }
#endif
    int iAccessAB = 0, iAccessBA = 0;
    cudaDeviceCanAccessPeer(&iAccessAB, iDevA, iDevB);
    cudaDeviceCanAccessPeer(&iAccessBA, iDevB, iDevA);
    if (!iAccessAB || !iAccessBA) {
        return false;
    }
    // Peer access persists in the context, so an earlier call may have enabled it already.
    int iCurrent;
    cudaGetDevice(&iCurrent);
    cudaSetDevice(iDevA);
    cudaError_t errAB = cudaDeviceEnablePeerAccess(iDevB, 0);
    cudaSetDevice(iDevB);
    cudaError_t errBA = cudaDeviceEnablePeerAccess(iDevA, 0);
    cudaSetDevice(iCurrent);
    if (errAB != cudaSuccess || errBA != cudaSuccess) {
        cudaGetLastError();
    }
    return (errAB == cudaSuccess || errAB == cudaErrorPeerAccessAlreadyEnabled) &&
           (errBA == cudaSuccess || errBA == cudaErrorPeerAccessAlreadyEnabled);
}


int GetGpuCount() {
    int iCudaDeviceCount = 0;
    cudaGetDeviceCount(&iCudaDeviceCount);
    return iCudaDeviceCount;
}
