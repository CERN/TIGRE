/*--------------------------------------------------------------------------
 This file is part of the TIGRE Toolbox

 Copyright (c) 2015, University of Bath and
                     CERN-European Organization for Nuclear Research
                     All rights reserved.

 License:            Open Source under BSD.
                     See the full license at
                     https://github.com/CERN/TIGRE/blob/master/LICENSE

 Contact:            tigre.toolbox@gmail.com
 Codes:              https://github.com/CERN/TIGRE/
 --------------------------------------------------------------------------
 */

#include "tigre_sycl.hpp"
#include "TIGRE_common.hpp"
#include <memory>
#include <mutex>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

static void printCpuFallbackOnce(const std::vector<sycl::device>& cpus){
    // Each TIGRE Python extension links tigre_sycl.cpp separately, so record
    // the current PID to print once per process.
    char pid[32];
#ifdef _WIN32
    snprintf(pid, sizeof(pid), "%d", _getpid());
    const char* seen = getenv("TIGRE_SYCL_CPU_NOTIFIED");
    if (seen && strcmp(seen, pid) == 0)
        return;
    _putenv_s("TIGRE_SYCL_CPU_NOTIFIED", pid);
#else
    snprintf(pid, sizeof(pid), "%ld", (long)getpid());
    const char* seen = getenv("TIGRE_SYCL_CPU_NOTIFIED");
    if (seen && strcmp(seen, pid) == 0)
        return;
    setenv("TIGRE_SYCL_CPU_NOTIFIED", pid, 1);
#endif
    fprintf(stderr, "TIGRE: no SYCL GPU found; running on the CPU (%s).\n",
            cpus[0].get_info<sycl::info::device::name>().c_str());
    fflush(stderr);
}

static std::vector<sycl::device> enumerateDevices(){
    std::vector<sycl::device> l0Gpus, otherGpus, cpus;
    for (const auto& platform : sycl::platform::get_platforms()){
        for (const auto& device : platform.get_devices()){
            if (device.is_gpu()){
                if (platform.get_backend()==sycl::backend::ext_oneapi_level_zero)
                    l0Gpus.push_back(device);
                else
                    otherGpus.push_back(device);
            } else if (device.is_cpu()){
                cpus.push_back(device);
            }
        }
    }
    std::vector<sycl::device> gpus = l0Gpus.empty() ? otherGpus : l0Gpus;
    if (gpus.empty()){
        if (!cpus.empty())
            printCpuFallbackOnce(cpus);
        return cpus;
    }
    // TIGRE_SYCL_ALL_DEVICES=1: list the CPU devices after the GPUs.
    const char* all = getenv("TIGRE_SYCL_ALL_DEVICES");
    if (all && *all && strcmp(all, "0") != 0)
        gpus.insert(gpus.end(), cpus.begin(), cpus.end());
    return gpus;
}

const std::vector<sycl::device>& tigreSyclDevices(){
    static const std::vector<sycl::device> devices = enumerateDevices();
    return devices;
}

sycl::queue& tigreSyclQueue(int iDeviceId){
    static std::vector<std::unique_ptr<sycl::queue>> queues;
    static std::mutex lock;
    std::lock_guard<std::mutex> guard(lock);
    const auto& devices = tigreSyclDevices();
    if (queues.size()!=devices.size())
        queues.resize(devices.size());
    if (!queues[iDeviceId])
        queues[iDeviceId].reset(new sycl::queue(devices[iDeviceId], sycl::property::queue::in_order()));
    return *queues[iDeviceId];
}

// TIGRE_SYCL_MEMORY_LIMIT_MB caps the memory TIGRE plans with on every device.
// For tests: it makes the split and out-of-memory paths reachable on any device.
static size_t memoryLimit(){
    const char* value = getenv("TIGRE_SYCL_MEMORY_LIMIT_MB");
    if (!value || !*value)
        return (size_t)-1;
    return (size_t)strtoull(value, nullptr, 10) << 20;
}

size_t tigreSyclMaxAlloc(int iDeviceId){
    const sycl::device& device = tigreSyclDevices()[iDeviceId];
    size_t global   = device.get_info<sycl::info::device::global_mem_size>();
    size_t maxAlloc = device.get_info<sycl::info::device::max_mem_alloc_size>();
    return std::min(std::min((size_t)((double)global*0.95), maxAlloc), memoryLimit());
}

// Free memory as reported by the runtime, if it reports it.
static bool reportedFreeMemory(int iDeviceId, size_t* free_bytes){
    const sycl::device& device = tigreSyclDevices()[iDeviceId];
    if (!device.has(sycl::aspect::ext_intel_free_memory))
        return false;
    try {
        *free_bytes = device.get_info<sycl::ext::intel::info::device::free_memory>();
        return true;
    } catch (const sycl::exception&) {
        // Not available after all (e.g. Level Zero Sysman not initialised).
        return false;
    }
}

size_t tigreSyclFreeMemory(int iDeviceId){
    size_t usable = tigreSyclMaxAlloc(iDeviceId);
    size_t free_bytes;
    if (reportedFreeMemory(iDeviceId, &free_bytes))
        usable = std::min(usable, free_bytes);
    return usable;
}

size_t tigreSyclTotalMemory(int iDeviceId){
    return tigreSyclDevices()[iDeviceId].get_info<sycl::info::device::global_mem_size>();
}

size_t tigreSyclCheckFreeMemory(int iDeviceId, const char* tag){
    size_t free_bytes;
    if (reportedFreeMemory(iDeviceId, &free_bytes) && free_bytes < tigreSyclTotalMemory(iDeviceId)/2){
        mexWarnMsgIdAndTxt(tag,"One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
    }
    return tigreSyclFreeMemory(iDeviceId);
}
