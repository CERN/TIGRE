
#ifndef GPUUTILS_HPP
#define GPUUTILS_HPP
//! @brief # of installed GPUs
int GetGpuCount();

//! @brief IDs of GPUs whose name is kacGPUName.
//! @note Call GetGpuCount and allocate sufficient memory for piDeviceIds.
//! @param [in] kacGPUName
//! @param [in, out] piDeviceIds. 
//! @param [in] iIdCountMax. Return value of GetGpuCount() 
int GetGpuIdArray(const char* kacGPUName, int* piDeviceIds, int iIdCountMax, char* pcMessage);

//! @brief GPU name of index iDeviceId. Allocate 128bytes for pcName before call.
void GetGpuName(int iDeviceId, char* pcName);

#endif  // GPUUTILS_HPP

