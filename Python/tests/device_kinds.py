"""Devices for the tests that compare two or more GPUs with one.

Those tests expect each device to give the same numbers, which holds for identical GPUs. With the
SYCL build and TIGRE_SYCL_ALL_DEVICES=1 the list can hold a GPU and the CPU device, which round
float32 arithmetic differently (and AwminTV amplifies that), so the results differ by device kind,
not by how the work is shared. Then same_kind() lists the first device once per device instead:
the same multi-device code runs, on identical devices.
"""

import functools
import os
import subprocess
import sys


@functools.lru_cache(maxsize=None)
def mixed_device_kinds():
    # As in test_projectors_gpu.py: GPUs and other devices are listed together exactly when the
    # list is longer than without TIGRE_SYCL_ALL_DEVICES.
    if os.environ.get("TIGRE_SYCL_ALL_DEVICES", "0") in ("", "0"):
        return False
    from tigre.utilities.gpu import GpuIds

    env = dict(os.environ)
    del env["TIGRE_SYCL_ALL_DEVICES"]
    code = "from tigre.utilities.gpu import getGpuNames; print(len(getGpuNames()))"
    out = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
    )
    return int(out.stdout.split()[-1]) < len(GpuIds())


def same_kind(devices):
    devices = list(devices)
    if len(devices) > 1 and mixed_device_kinds():
        return [devices[0]] * len(devices)
    return devices
