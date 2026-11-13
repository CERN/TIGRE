# Install

TIGRE runs on NVIDIA GPUs. It is available for Python and MATLAB; both versions
share the same CUDA code and the same set of algorithms.

:::{important}
TIGRE is **not** on PyPI. `pip install tigre` installs an unrelated package
with the same name, and the old `pytigre` package on PyPI is an outdated
v0.1 release.
:::

## Choose your route

| You want to… | Use | Needs CUDA toolkit and compiler? |
|---|---|---|
| Use TIGRE from Python on Linux | [Prebuilt wheel](python.md#prebuilt-wheels-linux) | No, only an NVIDIA driver |
| Use TIGRE from Python on Windows, or a Python version without a wheel | [Install from GitHub with pip](python.md#install-from-source) | Yes |
| Use the conda ecosystem | [CCPi conda package](python.md#conda-unofficial) (unofficial) | No |
| Use TIGRE from MATLAB | [Compile with `Compile.m`](matlab.md) | Yes |
| Change TIGRE's code | [Install from a clone](python.md#from-a-local-clone) | Yes |

## Requirements at a glance

- An NVIDIA GPU. When you compile TIGRE yourself, it is built for every GPU
  architecture your CUDA toolkit supports, so any GPU that your CUDA version
  supports will work.
- Python 3.11 to 3.14, or MATLAB R2016b or newer.
- To compile: the CUDA toolkit and a C++ compiler (MSVC on Windows, gcc on
  Linux).

```{toctree}
:maxdepth: 2

python
matlab
```
