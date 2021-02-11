
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


int GetGpuCount() {
    int iCudaDeviceCount = 0;
    cudaGetDeviceCount(&iCudaDeviceCount);
    return iCudaDeviceCount;
}
