
#include "gpuUtils.hpp"
#include "tigre_sycl.hpp"
#include <algorithm>
#include <string>
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

    iCudaDeviceCount = std::min(iCudaDeviceCount, iIdCountMax);
    iMessagePos += sprintf(pcMessage + iMessagePos, "Max GPUs: %d\n", iCudaDeviceCount);
    if (strlen(kacGPUName) == 0) {
        // Semi-compatible mode:
        //    Return all GPUs
        for (int iI = 0; iI < iCudaDeviceCount; ++iI) {
            piDeviceIds[iI] = iI;
        }
        return iCudaDeviceCount;
    }

    const auto& devices = tigreSyclDevices();
    int nMatch = 0;
    for (int iId = 0; iId < iCudaDeviceCount; ++iId) {
        std::string name = devices[iId].get_info<sycl::info::device::name>();
        iMessagePos += sprintf(pcMessage + iMessagePos, "propDevice.name = %s\n", name.c_str());
        if (strcmp(name.c_str(), kacGPUName) == 0) {
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
    std::string name = tigreSyclDevices()[iDeviceId].get_info<sycl::info::device::name>();
    memcpy(pcName, name.c_str(), std::min(name.size(), (size_t)127)*sizeof(char));
}


int GetGpuCount() {
    return (int)tigreSyclDevices().size();
}
