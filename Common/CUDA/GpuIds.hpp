
#ifndef GPUIDS_H
#define GPUIDS_H
struct GpuIds {
    int* m_piDeviceIds;
    int m_iCount;
    ~GpuIds();
    GpuIds();
    void SetIds(int iCount, int* piDeviceIds);
    int GetLength() const;
    void SetAllGpus(int iTotalDeviceCount);
    int& operator[](int iIndex);
    int operator[](int iIndex) const;
    bool AreEqualDevices() const;
};
#endif

