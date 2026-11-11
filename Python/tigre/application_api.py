"""Headless numerical services for desktop, notebook and batch applications.

Geometry recipes and persistence live in ``design_api``. These services bind
that geometry to the public operator algorithms; they contain no UI or solver
recurrences. Progress/cancellation occur at call boundaries, not inside CUDA
kernels or solver iterations.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import copy
import importlib.util
import time

import numpy as np


@dataclass(frozen=True)
class AlgorithmInfo:
    name: str
    family: str
    requires_tolerance: bool = False


@dataclass(frozen=True)
class DeviceInfo:
    index: int
    name: str


@dataclass(frozen=True)
class ReconstructionProgress:
    """A lifecycle stage, not an iteration counter or fractional estimate."""
    stage: str
    algorithm: str
    requested_iterations: int


@dataclass(frozen=True)
class ReconstructionResult:
    """Owned float32 volume and optional algorithm-specific diagnostics.

    ``diagnostics`` is None, an error/residual array, or a POCS history mapping
    according to the chosen algorithm/options. It is never fabricated into a
    universal residual history. ``elapsed_seconds`` measures the solver call.
    """
    volume: np.ndarray
    diagnostics: object
    algorithm: str
    requested_iterations: int
    elapsed_seconds: float


class ReconstructionCancelled(RuntimeError):
    """Cancellation observed before work or after a completed backend call."""


class BackendExecutionError(RuntimeError):
    """Backend execution failed; the original exception is chained as cause."""


_FAMILIES = {
    "algebraic": ("sart", "sirt", "os_sart", "fast_os_sart", "sart_tv", "os_sart_tv", "iterativereconalg"),
    "krylov": ("cgls", "lsqr", "lsmr", "hybrid_lsqr", "irn_tv_cgls", "hybrid_flsqr_tv", "ab_gmres", "ba_gmres"),
    "proximal_statistical": ("ista", "fista", "mlem"),
    "pocs": ("asd_pocs", "os_asd_pocs", "awasd_pocs", "os_awasd_pocs", "pcsd", "os_pcsd", "aw_pcsd", "os_aw_pcsd"),
}
_ALGORITHMS = tuple(AlgorithmInfo(name, family, family == "pocs")
                    for family, names in _FAMILIES.items() for name in names)
_ALIASES = {"ossart": "os_sart", "ossart_tv": "os_sart_tv"}


def available_algorithms():
    """Registered operator methods; availability does not test CUDA/hardware."""
    from . import algorithms
    return tuple(info for info in _ALGORITHMS if callable(getattr(algorithms, info.name + "_operator", None)))


def backend_capabilities():
    """Describe installed features without enumerating or initializing GPUs."""
    return {
        "algorithms": tuple(info.name for info in available_algorithms()),
        "explicit_surfaces": True,
        "surface_fdk": importlib.util.find_spec("tigre.algorithms.fdk_surface") is not None,
        "iteration_progress": False,
        "in_flight_cancellation": False,
        "cuda_runtime_checked": False,
    }


def available_devices():
    """Enumerate zero-based CUDA device indices on explicit request only."""
    from .utilities.gpu import getGpuNames
    try:
        return tuple(DeviceInfo(index, name) for index, name in enumerate(getGpuNames()))
    except Exception as error:
        raise BackendExecutionError("CUDA device discovery failed: " + str(error)) from error


def create_operator(geometry, *, surface=None, angles=None, gpuids=None, n_slabs=None):
    """Bind an owned conventional or explicit-surface acquisition.

    Supply exactly one of surface/angles. Geometry is not silently moved to an
    object frame, resized, rebinned or otherwise changed. Use the corresponding
    design_api helper explicitly before binding moving-object acquisitions.
    """
    from .utilities.projection_operator import SurfaceOperator, TigreOperator
    if (surface is None) == (angles is None):
        raise ValueError("supply exactly one of surface or angles")
    if surface is not None:
        return SurfaceOperator(geometry, surface, gpuids=gpuids, n_slabs=n_slabs)
    if n_slabs is not None:
        raise ValueError("n_slabs belongs to explicit-surface operators")
    return TigreOperator(geometry, angles, gpuids=gpuids)


def forward_project(volume, geometry=None, *, operator=None, surface=None, angles=None,
                    gpuids=None, n_slabs=None):
    """Project through the selected operator, preserving caller inputs."""
    from .utilities.projection_operator import _array
    op = _bind(geometry, operator, surface, angles, gpuids, n_slabs)
    image = _array(volume, op.volume_shape, "volume")
    try:
        projected = op.forward(image)
    except Exception as error:
        raise BackendExecutionError("Forward projection failed: " + str(error)) from error
    return _array(projected, op.projection_shape, "forward result", op.valid_mask)


def _bind(geometry, operator, surface, angles, gpuids, n_slabs):
    if operator is not None:
        if any(value is not None for value in (geometry, surface, angles, gpuids, n_slabs)):
            raise ValueError("with operator=, omit geometry, acquisition and device settings")
        return operator
    if geometry is None:
        raise ValueError("geometry is required when operator is not supplied")
    return create_operator(geometry, surface=surface, angles=angles, gpuids=gpuids, n_slabs=n_slabs)


def reconstruct(projections, geometry=None, *, algorithm="cgls", niter=20,
                operator=None, surface=None, angles=None, gpuids=None, n_slabs=None,
                options=None, progress=None, cancelled=None):
    """Reconstruct with a registered operator method and return owned results.

    ``progress(event)`` is called synchronously on the caller's thread with
    stages validating, running and completed. UI clients marshal events onto
    their own event loop. Callback exceptions propagate unchanged.

    ``cancelled()`` is checked before validation, before execution and after
    execution. An in-flight solver/native call finishes before cancellation
    is observed; its result is then discarded. There is no per-iteration
    cancellation, solver restart or fabricated iteration progress.

    ``options`` uses the selected operator algorithm's documented parameters.
    POCS requires explicit ``maxl2err``. FDK/FBP are outside this iterative API.
    """
    from . import algorithms
    from .utilities.projection_operator import prepare_projections, _array
    if not isinstance(algorithm, str):
        raise TypeError("algorithm must be a registered method name")
    name = _ALIASES.get(algorithm.lower(), algorithm.lower())
    info = next((item for item in available_algorithms() if item.name == name), None)
    if info is None:
        raise ValueError("unsupported iterative algorithm: " + name)
    if isinstance(niter, (bool, np.bool_)) or not isinstance(niter, (int, np.integer)) or niter < 0:
        raise ValueError("niter must be a nonnegative integer")
    niter = int(niter)
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    if cancelled is not None and not callable(cancelled):
        raise TypeError("cancelled must be callable or None")
    if options is not None and not isinstance(options, Mapping):
        raise TypeError("options must be a mapping or None")
    kwargs = dict(options or {})
    reserved = set(kwargs) & {"operator", "niter", "geo", "geometry", "angles", "surface", "gpuids", "n_slabs"}
    if reserved:
        raise ValueError("pass geometry, device and iteration settings outside options: " + str(sorted(reserved)))
    if info.requires_tolerance:
        epsilon = kwargs.get("maxl2err")
        if epsilon is None or np.ndim(epsilon) != 0 or not np.isrealobj(epsilon) or not np.isfinite(epsilon) or epsilon < 0:
            raise ValueError(name + " requires a finite nonnegative absolute maxl2err")

    def checkpoint(stage):
        if cancelled is not None and cancelled():
            raise ReconstructionCancelled("Reconstruction cancelled " + stage)

    def report(stage):
        if progress is not None:
            progress(ReconstructionProgress(stage, name, niter))

    checkpoint("before validation")
    report("validating")
    op = _bind(geometry, operator, surface, angles, gpuids, n_slabs)
    data = prepare_projections(projections, op)
    checkpoint("before execution")
    report("running")
    checkpoint("before execution")
    started = time.perf_counter()
    try:
        value = getattr(algorithms, name + "_operator")(data, op, niter=niter, **kwargs)
    except Exception as error:
        raise BackendExecutionError(name + " reconstruction failed: " + str(error)) from error
    elapsed = time.perf_counter() - started
    checkpoint("after execution; completed result discarded")
    image, diagnostics = value if isinstance(value, tuple) else (value, None)
    image = _array(image, op.volume_shape, "reconstructed volume")
    result = ReconstructionResult(image, copy.deepcopy(diagnostics), name, niter, elapsed)
    report("completed")
    return result


__all__ = [
    "AlgorithmInfo", "DeviceInfo", "ReconstructionProgress", "ReconstructionResult",
    "ReconstructionCancelled", "BackendExecutionError", "available_algorithms",
    "backend_capabilities", "available_devices", "create_operator", "forward_project", "reconstruct",
]
