#include "GpuIds.hpp"
#include <stdlib.h>
#include <string.h>
#include <string>
#include "tigre_sycl.hpp"

GpuIds::~GpuIds() {
    free(m_piDeviceIds); m_piDeviceIds = nullptr;
    m_iCount = 0;
}
GpuIds::GpuIds() : m_piDeviceIds (nullptr), m_iCount(0) {

}
void GpuIds::SetIds(int iCount, int* piDeviceIds) {
    if (iCount > 0 && piDeviceIds != 0) {
        if (m_piDeviceIds) {
            free(m_piDeviceIds); m_piDeviceIds = nullptr;
            m_iCount = 0;
        }
        m_piDeviceIds = (int*)malloc(iCount * sizeof(int));
        if (m_piDeviceIds) {
            for (int iI = 0; iI < iCount; ++iI) {
                m_piDeviceIds[iI] = piDeviceIds[iI];
            }
            m_iCount = iCount;
        }
    }
}

int GpuIds::GetLength() const {
    return m_iCount;
}
int& GpuIds::operator[](int iIndex){
    return m_piDeviceIds[iIndex];
}
int GpuIds::operator[](int iIndex) const {
    return m_piDeviceIds[iIndex];
}

void GpuIds::SetAllGpus(int iTotalDeviceCount) {
    // Set all GPUs for compatibility
    // Makeup valid GpuIds.
    int* aiIds = nullptr;
    if (iTotalDeviceCount == 0) {
        (int*)malloc(iTotalDeviceCount*sizeof(int));
        for (int iI = 0; iI < iTotalDeviceCount; ++iI) {
            aiIds[iI] = iI;
        }
    }
    SetIds(iTotalDeviceCount, aiIds);
    free(aiIds); aiIds = 0;    
}

bool GpuIds::AreEqualDevices() const {
    int deviceCount = this->GetLength();
    const auto& devices = tigreSyclDevices();
    std::string devicename;
    for (int dev = 0; dev < deviceCount; dev++) {
        std::string name = devices[m_piDeviceIds[dev]].get_info<sycl::info::device::name>();
        if (dev>0) {
            if (devicename != name) {
                return false;
            }
        }
        devicename = name;
    }
    return true;
}
