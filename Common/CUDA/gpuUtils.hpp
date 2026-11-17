
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

//! @brief Whether the multi-GPU halo exchange between iDevA and iDevB should copy peer-to-peer instead of
//! staging through pinned host memory. Only says yes where peer copies were measured to be faster, and
//! enables peer access (both ways) only in that case, since enabling it slows later allocations.
//! @param [in] bManySmallCopies true if the caller issues several concurrent copies per pair (tv_proximal:
//! u, px, py, pz), which serialize on ROCm.
bool UsePeerHaloExchange(int iDevA, int iDevB, bool bManySmallCopies);

#endif  // GPUUTILS_HPP

