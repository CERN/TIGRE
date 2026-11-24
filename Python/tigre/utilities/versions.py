"""TIGRE's version, and the versions of what it runs on (for bug reports)."""
import platform
import sys


def _tigre_version():
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover
        return "unknown"
    try:
        return version("tigre")
    except PackageNotFoundError:
        # run from a source tree without installing it
        return "unknown"


def _gpu_names():
    from .gpu import getGpuNames

    return getGpuNames()


def show_versions():
    """Print the versions of TIGRE, Python, NumPy and the GPUs.

    Paste the output into bug reports: https://github.com/CERN/TIGRE/issues

    Examples
    --------
    >>> import tigre
    >>> tigre.show_versions()  # doctest: +SKIP
    TIGRE    3.3.1
    Python   3.12.3 (Linux-6.8.0-x86_64-with-glibc2.39)
    NumPy    2.1.3
    GPUs     NVIDIA GeForce RTX 4090
    """
    import numpy

    lines = [
        ("TIGRE", _tigre_version()),
        ("Python", f"{platform.python_version()} ({platform.platform()})"),
        ("NumPy", numpy.__version__),
    ]
    try:
        names = _gpu_names()
        lines.append(("GPUs", ", ".join(names) if names else "none found"))
    except Exception as err:  # the GPU query itself may be what is broken
        lines.append(("GPUs", f"could not query them ({err})"))
    torch = sys.modules.get("torch")
    if torch is not None:
        lines.append(("PyTorch", torch.__version__))
    for name, value in lines:
        print(f"{name:<8} {value}")
