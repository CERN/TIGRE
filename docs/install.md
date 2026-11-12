# Install

TIGRE runs on NVIDIA GPUs and needs the CUDA toolkit and a C++ compiler, as it
compiles its CUDA code on your machine.

:::{important}
TIGRE is **not** on PyPI. `pip install tigre` installs an unrelated package
with the same name.
:::

::::{tab-set}

:::{tab-item} Python
Install straight from GitHub:

```sh
# latest master
pip install git+https://github.com/CERN/TIGRE.git

# a specific release
pip install git+https://github.com/CERN/TIGRE.git@vX.Y.Z
```

Replace `vX.Y.Z` with a release tag. The latest release is **v{{ release }}**.

Prebuilt Linux wheels are attached to each
[GitHub release](https://github.com/CERN/TIGRE/releases).

Full requirements and troubleshooting:
[Python installation guide](https://github.com/CERN/TIGRE/blob/master/Frontispiece/python_installation.md).
:::

:::{tab-item} MATLAB
See the
[MATLAB installation guide](https://github.com/CERN/TIGRE/blob/master/Frontispiece/MATLAB_installation.md).
:::
::::

:::{note}
The full installation guides are moving into this site. Until then, the links
above point to the current guides in the repository.
:::
