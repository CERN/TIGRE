"""Bind existing iterative recurrences to explicit projection operators.

View tokens below are internal indices, never physical angles and never sent
to native projectors. Conventional calls retain their original projector calls.
"""

import copy
from types import SimpleNamespace
import numpy as np
from tigre.utilities.projection_operator import prepare_projections, _array
from .operator_algorithms import _initial, _iterations, _valid_samples


def _validated_options(cls, kwargs, shape):
    defaults = dict(
        blocksize=20,
        lmbda=1.0,
        lmbda_red=1.0,
        OrderStrategy=None,
        Quameasopts=None,
        init=None,
        verbose=False,
        noneg=True,
        mu_max=None,
        support=None,
        computel2=False,
        dataminimizing="art_data_minimizing",
        name=cls.__name__,
        niter_outer=4,
        restart=True,
    )
    allowed = set(defaults) | {
        "W",
        "V",
        "tviter",
        "tvlambda",
        "hyper",
        "fista_p",
        "fista_q",
        "backprojector",
        "sup_kw_warning",
    }
    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(f"Unsupported operator options: {sorted(unknown)}")
    defaults.update(kwargs)
    for key in ("tviter", "niter_outer"):
        if key in defaults:
            value = defaults[key]
            if (
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer))
                or value < 0
            ):
                raise ValueError(f"{key} must be a nonnegative integer")
    if cls.__name__ == "IRN_TV_CGLS" and defaults["niter_outer"] == 0:
        raise ValueError("niter_outer must be positive")
    if "tvlambda" in defaults:
        value = defaults["tvlambda"]
        if (
            not np.isfinite(value)
            or value < 0
            or (value == 0 and cls.__name__ in ("SART_TV", "OS_SART_TV"))
        ):
            raise ValueError("Invalid TV regularization parameter")
    if "hyper" in defaults:
        value = defaults["hyper"]
        if (isinstance(value, str) and value != "auto") or (
            not isinstance(value, str) and (not np.isfinite(value) or value <= 0)
        ):
            raise ValueError("hyper must be auto or finite and positive")
    for key in ("W", "V"):
        if defaults.get(key) is not None:
            raise ValueError("Operator normalization is computed from the selected operator")

    for key in ("lmbda", "lmbda_red", "mu_max"):
        value = defaults.get(key)
        if value is None and key == "mu_max":
            continue
        if np.ndim(value) != 0 or not np.isrealobj(value) or not np.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    size = defaults["blocksize"]
    if isinstance(size, (bool, np.bool_)) or not isinstance(size, (int, np.integer)) or size < 1:
        raise ValueError("blocksize must be a positive integer")
    if defaults["OrderStrategy"] not in (None, "ordered", "random", "random2"):
        raise ValueError("Operator ordering supports ordered or random view indices")
    for key in ("fista_p", "fista_q"):
        if key in defaults:
            value = defaults[key]
            if np.ndim(value) != 0 or not np.isfinite(value) or value < 0:
                raise ValueError(f"{key} must be finite and nonnegative")
    if defaults["support"] is not None:
        support = _array(defaults["support"], shape, "support", copy_data=False)
        if np.any((support != 0) & (support != 1)):
            raise ValueError("support must be binary")
    if defaults["init"] is not None:
        if isinstance(defaults["init"], str):
            raise ValueError("Operator initializers must be volumes, not FDK/multigrid")
        _array(defaults["init"], shape, "init", copy_data=False)
    return defaults


class OperatorReconMixin:
    def _operator_lipschitz(self, seed=0):
        from .operator_algorithms import _squared_norm
        rng=np.random.default_rng(seed)
        x=rng.standard_normal(tuple(self.geo.nVoxel)).astype(np.float32)
        x /= np.float32(np.sqrt(_squared_norm(x)))
        scale=1.0
        for _ in range(8):
            y=self.operator.adjoint(self.operator.forward(x))
            scale=np.sqrt(_squared_norm(y))
            if scale==0:
                return 1.0
            x=y/np.float32(scale)
        return 2.1*scale

    def apply_constraints(self):
        if self.noneg:
            np.maximum(self.res, 0, out=self.res)
        if self.mu_max is not None:
            np.minimum(self.res, self.mu_max, out=self.res)
        if self.support is not None:
            self.res *= self.support

    def apply_terminal_constraints(self):
        if self.mu_max is not None:
            np.minimum(self.res, self.mu_max, out=self.res)
        if self.support is not None:
            self.res *= self.support

    def _initialize_operator(self, proj, operator, niter, kwargs):
        self.operator = operator
        self.proj = prepare_projections(proj, operator)
        self.niter = _iterations(niter)
        shape = tuple(operator.volume_shape)
        views = operator.projection_shape[0]
        spacing = np.asarray(getattr(getattr(operator, "geo", None), "dVoxel", np.ones(3)))
        self.geo = SimpleNamespace(
            nVoxel=np.asarray(shape, dtype=int),
            dVoxel=spacing,
            sVoxel=np.asarray(shape) * spacing,
            nDetector=np.asarray(operator.projection_shape[1:], dtype=int),
        )
        self.angles = np.zeros((views, 3), np.float32)
        self.angles[:, 0] = np.arange(views)
        defaults = _validated_options(type(self), kwargs, shape)
        self.__dict__.update(defaults)
        self.gpuids = copy.deepcopy(getattr(operator, "gpuids", None))
        if self.gpuids is not None:
            from tigre.utilities.Ax_surface import _devices
            from tigre.utilities.gpu import GpuIds

            devices = GpuIds()
            devices.devices = _devices(self.gpuids)
            self.gpuids = devices
        if isinstance(self.init, str):
            raise ValueError("Operator initializers must be volumes, not FDK/multigrid")
        self.res = _initial(self.init, shape)
        if self.support is not None:
            self.support = _array(self.support, shape, "support")
            if np.any((self.support != 0) & (self.support != 1)):
                raise ValueError("support must be binary")
        if self.mu_max is not None and (not np.isfinite(self.mu_max) or self.mu_max < 0):
            raise ValueError("mu_max must be finite and nonnegative")
        self.apply_constraints()
        if self.OrderStrategy not in (None, "ordered", "random", "random2"):
            raise ValueError("Operator ordering supports ordered or random view indices")
        size = self.blocksize
        if (
            isinstance(size, (bool, np.bool_))
            or not isinstance(size, (int, np.integer))
            or size < 1
        ):
            raise ValueError("blocksize must be a positive integer")
        order = np.arange(views)
        if self.OrderStrategy in ("random", "random2"):
            np.random.shuffle(order)
        count = max(1, -(-views // size))
        self.angle_index = [order[k::count] for k in range(count)]
        self.angleblocks = [self.angles[idx] for idx in self.angle_index]
        self.W = None
        self.V = None
        self._operator_cache = {}
        if isinstance(self.Quameasopts, str):
            self.Quameasopts = [self.Quameasopts]
        self.lq = np.zeros((0 if self.Quameasopts is None else len(self.Quameasopts), self.niter))
        self.l2l = np.zeros((1, self.niter), np.float32)
        if hasattr(self, "backproject"):
            self.backproject = self._matched_backproject

    def _view_operator(self, angles):
        tokens = np.asarray(angles).reshape(-1, 3)[:, 0]
        indices = tokens.astype(np.intp)
        if np.any(indices != tokens) or np.any(indices < 0) or np.any(indices >= len(self.angles)):
            raise ValueError("Invalid internal operator view indices")
        if np.array_equal(indices, np.arange(len(self.angles))):
            return self.operator
        key = tuple(indices)
        if key not in self._operator_cache:
            self._operator_cache[key] = self.operator.subset(indices)
        return self._operator_cache[key]

    def _project_forward(self, legacy, image, geo, angles, *args, **kwargs):
        if not hasattr(self, "operator"):
            return legacy(image, geo, angles, *args, **kwargs)
        sub = self._view_operator(angles)
        image = _array(image, sub.volume_shape, "volume", copy_data=False)
        return _array(sub.forward(image), sub.projection_shape, "operator forward", sub.valid_mask)

    def _project_adjoint(self, legacy, data, geo, angles, *args, **kwargs):
        if not hasattr(self, "operator"):
            return legacy(data, geo, angles, *args, **kwargs)
        sub = self._view_operator(angles)
        return _array(
            sub.adjoint(prepare_projections(data, sub)), sub.volume_shape, "operator adjoint"
        )

    def _matched_backproject(self, data, geo, angles, **kwargs):
        return self._project_adjoint(None, data, geo, angles, **kwargs)

    def _volume_denoise(self, legacy, image, iterations=50, lmbda=15.0, gpuids=None):
        if not hasattr(self, "operator"):
            return legacy(image, iterations, lmbda, gpuids)
        if iterations == 0 or np.ptp(image) == 0:
            return image.copy()
        if min(image.shape) < 2:
            raise ValueError("Native volume TV needs at least two voxels along every axis")
        result = np.asarray(legacy(image, iterations, lmbda, gpuids), dtype=np.float32)
        if self.noneg:
            np.maximum(result, 0, out=result)
        if self.mu_max is not None:
            np.minimum(result, self.mu_max, out=result)
        if self.support is not None:
            result *= self.support
        return result

    def _operator_tv_precondition(self, weight, rhs):
        from .operator_algorithms import _squared_norm

        rhs = np.asarray(rhs, np.float32).reshape(self.res.shape).copy()
        rhs -= np.float32(np.mean(rhs, dtype=np.float64))
        x = np.zeros_like(rhs)
        r = rhs.copy()
        p = r.copy()
        gamma = _squared_norm(r)
        initial = gamma
        for _ in range(50):
            if gamma <= 1e-12 * initial or gamma == 0:
                break
            q = self.Ltx(weight, self.Lx(weight, p))
            denominator = float(np.sum(p.astype(np.float64) * q.astype(np.float64)))
            if denominator <= 0:
                break
            alpha = gamma / denominator
            x += np.float32(alpha) * p
            r -= np.float32(alpha) * q
            next_gamma = _squared_norm(r)
            p = r + np.float32(next_gamma / gamma) * p
            gamma = next_gamma
        return x.ravel()

    def _operator_data_step(self):
        for k, angle in enumerate(self.angleblocks):
            self.update_image(self.geo, angle, k)
            self.apply_constraints()

    def _operator_update(self, iteration):
        idx = self.angle_index[iteration]
        sub = self._view_operator(self.angleblocks[iteration])
        key = ("row", tuple(idx))
        if key not in self._operator_cache:
            row = sub.forward(np.ones(sub.volume_shape, np.float32))
            w = np.zeros_like(row)
            np.divide(1.0, row, out=w, where=row > 1e-6)
            self._operator_cache[key] = w
        col = sub.adjoint(_valid_samples(sub))
        v = np.zeros_like(col)
        np.divide(1.0, col, out=v, where=col > 1e-6)
        residual = (self.proj[idx] - sub.forward(self.res)) * self._operator_cache[key]
        self.res += np.float32(self.lmbda) * v * sub.adjoint(residual.astype(np.float32))


def run_existing(cls, proj, operator, niter, **kwargs):
    from .operator_algorithms import _squared_norm

    niter = _iterations(niter)
    if kwargs.get("backprojector", "matched") != "matched":
        raise ValueError(
            "Operator GMRES requires its adjoint; FDK backprojectors are geometry-specific"
        )
    if "nonneg" in kwargs:
        if "noneg" in kwargs:
            raise TypeError("specify only one of noneg and nonneg")
        kwargs["noneg"] = kwargs.pop("nonneg")
    init = kwargs.get("init")
    b = prepare_projections(proj, operator)
    options = _validated_options(cls, kwargs, operator.volume_shape)
    if cls.__name__ == "MLEM" and np.any(b < 0):
        raise ValueError("MLEM requires nonnegative valid data")
    if niter == 0:
        result = _initial(init, operator.volume_shape)
        if cls.__name__ == "MLEM":
            if init is None:
                result += 1
            np.maximum(result, 0, out=result)
        if kwargs.get("computel2") or kwargs.get("Quameasopts") is not None:
            quality = kwargs.get("Quameasopts")
            count = 0 if quality is None else (1 if isinstance(quality, str) else len(quality))
            return result, np.zeros((count + bool(kwargs.get("computel2")), 0))
        return result
    if _squared_norm(b) == 0 and init is None:
        result = np.zeros(operator.volume_shape, np.float32)
        if kwargs.get("computel2") or kwargs.get("Quameasopts") is not None:
            quality = kwargs.get("Quameasopts")
            count = 0 if quality is None else (1 if isinstance(quality, str) else len(quality))
            return result, np.zeros((count + bool(kwargs.get("computel2")), niter))
        return result
    if cls.__name__ in ("LSQR", "LSMR", "hybrid_LSQR", "AB_GMRES", "BA_GMRES", "hybrid_fLSQR_TV"):
        candidate = _initial(init, operator.volume_shape)
        if options["noneg"]:
            np.maximum(candidate, 0, out=candidate)
        if options["mu_max"] is not None:
            np.minimum(candidate, options["mu_max"], out=candidate)
        if options["support"] is not None:
            candidate *= options["support"]
        residual = b if init is None else b - operator.forward(candidate)
        stationary = _squared_norm(residual) == 0
        if not stationary:
            stationary = _squared_norm(operator.adjoint(residual)) == 0
        if stationary:
            if kwargs.get("computel2") or kwargs.get("Quameasopts") is not None:
                quality = kwargs.get("Quameasopts")
                count = 0 if quality is None else (1 if isinstance(quality, str) else len(quality))
                errors = np.zeros((count + bool(kwargs.get("computel2")), niter))
                if kwargs.get("computel2"):
                    errors[0] = np.sqrt(_squared_norm(residual))
                return candidate, errors
            return candidate
    tokens = np.zeros((operator.projection_shape[0], 3), np.float32)
    tokens[:, 0] = np.arange(len(tokens))
    algorithm = cls(b, None, tokens, niter, _projection_operator=operator, **kwargs)
    if cls.__name__ == "hybrid_fLSQR_TV":
        # Reuse the validated finite-difference/transpose pair rather than
        # the historical flexible solver's off-by-one boundary stencil.
        from types import MethodType
        from .krylov_subspace_algorithms import IRN_TV_CGLS

        for method in ("__D__", "__Dt__", "__build_weights__"):
            setattr(algorithm, method, MethodType(getattr(IRN_TV_CGLS, method), algorithm))
        algorithm.Lx = lambda weight, image: IRN_TV_CGLS.Lx(
            algorithm, weight, np.reshape(image, operator.volume_shape)
        )
        algorithm.Ltx = lambda weight, image: IRN_TV_CGLS.Ltx(
            algorithm, weight, np.reshape(image, (3, *operator.volume_shape))
        )
    algorithm.run_main_iter()
    algorithm.apply_terminal_constraints()
    result = np.asarray(algorithm.getres(), dtype=np.float32)
    if not np.isfinite(result).all():
        raise FloatingPointError(f"{cls.__name__} produced nonfinite values")
    if algorithm.computel2 or algorithm.Quameasopts is not None:
        return result, algorithm.geterrors()
    return result
