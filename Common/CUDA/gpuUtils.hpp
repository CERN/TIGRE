
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

//! @brief Whether the multi-GPU halo exchange between iDevA and iDevB should copy peer-to-peer instead of staging
//! through pinned host memory. Enables peer access (both ways) only when it says yes; on HIP it always says no.
bool UsePeerHaloExchange(int iDevA, int iDevB);

//! @brief How the TV functions (minTV, AwminTV, tvdenoise) share the z slices of an image between deviceCount GPUs
//! and splits rounds: deviceCount*splits pieces of slices_per_split slices, the last one possibly shorter. Pieces
//! must not be empty: the code computes a piece's size as nz minus where it starts, which wraps around (unsigned)
//! past the end of the image, and launches kernels with no blocks for an empty piece. A plain ceil(nz/pieces) can
//! leave the last pieces empty (e.g. 5 slices on 4 GPUs: 2, 2, 1, -1). Given the largest piece that fits
//! (slices_per_split) and at most deviceCount GPUs, this picks the fewest rounds, then the most GPUs, then the
//! largest pieces, with every piece non-empty. With deviceCount=1 and slices_per_split=1 it always succeeds.
inline void tvPiecesWithoutEmpty(unsigned long nz, int* deviceCount, unsigned int* splits, unsigned int* slices_per_split) {
    int bestD = 1;
    unsigned int bestS = 1;
    unsigned long bestRounds = nz;
    for (int d = *deviceCount; d >= 1; d--) {
        for (unsigned int s = *slices_per_split; s >= 1; s--) {
            const unsigned long pieces = (nz + s - 1) / s;  // all non-empty, the last one possibly shorter
            if (pieces % d != 0) continue;
            const unsigned long rounds = pieces / d;
            if (rounds < bestRounds || (rounds == bestRounds && (d > bestD || (d == bestD && s > bestS)))) {
                bestRounds = rounds; bestD = d; bestS = s;
            }
        }
    }
    *deviceCount = bestD;
    *splits = (unsigned int)bestRounds;
    *slices_per_split = bestS;
}

#endif  // GPUUTILS_HPP

