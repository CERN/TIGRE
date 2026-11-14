"""Docstrings of the iterative algorithms.

The algorithms share most of their arguments, so their docstrings are put
together here from shared pieces (numpydoc format) and attached to each
algorithm function by ``iterative_recon_alg.decorator``.
"""

_PARAMETERS = """\
Parameters
----------
proj : numpy.ndarray
    Projections (line integrals), float32, shape ``(len(angles), *geo.nDetector)``.
geo : tigre.utilities.geometry.Geometry
    The scanner geometry.
angles : numpy.ndarray
    Projection angles in radians: shape ``(N,)`` for rotations around z, or
    ``(N, 3)`` for ZYZ Euler angles.
niter : int
    Number of iterations."""

_RETURNS = """\
Returns
-------
img : numpy.ndarray
    The reconstruction, float32, shape ``geo.nVoxel``.
errors : numpy.ndarray
    Only returned when ``computel2=True`` or ``Quameasopts`` is given: one row
    per measure (the data residual norm first, if requested) and one column per
    iteration."""

_COMMON = """\
init : {None, "FDK", "multigrid"} or numpy.ndarray, optional
    Starting image: zeros (default), an FDK reconstruction, a multigrid
    initialisation (OS-SART on coarser grids), or an array of shape
    ``geo.nVoxel``.
verbose : bool, optional
    Print progress and an estimate of the time left. Default True.
Quameasopts : str or list of str, optional
    Image quality measures to record at every iteration, comparing each
    iterate with the previous one, for example ``["RMSE", "CC"]``. See
    :func:`tigre.utilities.Measure_Quality.Measure_Quality`.
computel2 : bool, optional
    Also record the data residual norm :math:`\\|Ax - b\\|_2` at every
    iteration. Default False.
gpuids : tigre.utilities.gpu.GpuIds, optional
    GPUs to use. Default: all available."""

_SART = """\
lmbda : float, optional
    Relaxation parameter of the update. Default 1.
lmbda_red : float, optional
    Factor by which ``lmbda`` is multiplied after every iteration.
    Default {lmbda_red}.
OrderStrategy : {{"ordered", "random", "angularDistance"}}, optional
    Order in which the subsets of projections are used. Default: the order
    of ``angles``.
noneg : bool, optional
    Set negative values to zero after every update. Default True."""

_BLOCKSIZE = """\
blocksize : int, optional
    Number of projections used in each update. Default 20."""

_TV_STEP = """\
tviter : int, optional
    Iterations of TV denoising after each data step. Default {tviter}.
tvlambda : float, optional
    Data-fidelity weight of the TV denoising: larger values keep the image
    closer to the data step's result, so they denoise less. Default {tvlambda}."""

_POCS = """\
lmbda : float, optional
    Step size of the data (SART-type) step. Default 1.
lmbda_red : float, optional
    Factor by which ``lmbda`` is multiplied after every iteration. Default 0.99.
alpha : float, optional
    Step size of the TV gradient descent. Default 0.002.
alpha_red : float, optional
    Factor by which ``alpha`` is reduced when needed. Default 0.95.
rmax : float, optional
    Maximum ratio between the TV step and the data step. Default 0.95.
maxl2err : float, optional
    Data error the image is allowed to keep, :math:`\\|Ax - b\\|_2`
    (epsilon in the paper). Default: 20% of the data error of an FDK
    reconstruction.
tviter : int, optional
    TV gradient-descent steps per iteration. Default 20."""

_PCSD = """\
lmbda : float, optional
    Step size of the data (SART-type) step. Default 1.
lmbda_red : float, optional
    Factor by which ``lmbda`` is multiplied after every iteration. Default 0.99.
maxl2err : float, optional
    Data error the image is allowed to keep, :math:`\\|Ax - b\\|_2`.
    Default: 20% of the data error of an FDK reconstruction.
tviter : int, optional
    TV gradient-descent steps per iteration. Default 20."""

_AW = """\
delta : float, optional
    Edge sensitivity of the adaptive weights: closer to 0 preserves more
    edges. Default -0.005."""

_PRIOR = """\
prior : numpy.ndarray
    Prior image of the same object, shape ``geo.nVoxel``. Required.
prior_ratio : float, optional
    Weight of the prior image term against the TV term. Default 0.2."""

_RESTART = """\
restart : bool, optional
    When the residual grows (loss of orthogonality), restart from the current
    image instead of stopping. Default True."""

_HYBRID = """\
lmbda : float, optional
    Tikhonov regularisation parameter. If not given, it is chosen at every
    iteration: by the discrepancy principle if ``noiselevel`` is given,
    otherwise by generalised cross validation (GCV).
noiselevel : float, optional
    Relative noise level of the data, :math:`\\|e\\| / \\|b\\|` (for example
    0.02 for 2%), for the discrepancy principle."""

_KRYLOV_NOTE = """\
Krylov methods use the ``matched`` backprojector. They converge quickly but
start to fit the noise after a while, so use few iterations (often 10 to 30)
or a regularised variant. See :doc:`/user-guide/projectors` for why the
projector pair is not an exact transpose."""

_POCS_NOTE = """\
The algorithm alternates a SART-type data step with TV gradient-descent
steps, and stops early when the image no longer changes or the step size
becomes too small."""


def _section(title, *bodies):
    bodies = [b for b in bodies if b]
    if not bodies:
        return ""
    return f"{title}\n{'-' * len(title)}\n" + "\n".join(bodies)


def _build(summary, *, extra=(), notes="", see_also="", references=""):
    parts = [summary.strip(), _PARAMETERS, _RETURNS, _section("Other Parameters", *extra, _COMMON)]
    if see_also:
        parts.append(_section("See Also", see_also.strip()))
    if notes:
        parts.append(_section("Notes", notes.strip()))
    if references:
        parts.append(_section("References", references.strip()))
    return "\n\n".join(p for p in parts if p) + "\n"


_SART_ONE = _SART.format(lmbda_red=1)
_LSQR_REF = ".. [1] C. C. Paige and M. A. Saunders, \"LSQR: An algorithm for sparse linear equations and sparse least squares\", ACM Trans. Math. Softw. 8 (1982). https://doi.org/10.1145/355984.355989"
_LSMR_REF = ".. [1] D. C.-L. Fong and M. A. Saunders, \"LSMR: An iterative algorithm for sparse least-squares problems\", SIAM J. Sci. Comput. 33 (2011). https://doi.org/10.1137/10079687X"
_ASD_REF = ".. [1] E. Y. Sidky and X. Pan, \"Image reconstruction in circular cone-beam computed tomography by constrained, total-variation minimization\", Phys. Med. Biol. 53 (2008). https://doi.org/10.1088/0031-9155/53/17/021"
_FISTA_REF = ".. [1] A. Beck and M. Teboulle, \"A fast iterative shrinkage-thresholding algorithm for linear inverse problems\", SIAM J. Imaging Sci. 2 (2009). https://doi.org/10.1137/080716542"

DOCS = {
    # SART family
    "sart": _build(
        "SART: simultaneous algebraic reconstruction technique.\n\n"
        "Updates the image after each single projection. Converges in the fewest "
        "iterations of the SART family, but each iteration is the slowest.",
        extra=(_SART_ONE,),
        see_also="ossart, sirt",
    ),
    "sirt": _build(
        "SIRT: simultaneous iterative reconstruction technique.\n\n"
        "Updates the image once per iteration, using all projections. Each "
        "iteration is fast, but many are needed.",
        extra=(_SART_ONE,),
        see_also="ossart, sart",
    ),
    "ossart": _build(
        "OS-SART: ordered-subsets SART.\n\n"
        "Updates the image after each subset of ``blocksize`` projections, "
        "between SART (one projection) and SIRT (all). Usually the best choice "
        "of the three.",
        extra=(_BLOCKSIZE, _SART_ONE),
        see_also="sart, sirt, ossart_tv",
    ),
    "sart_tv": _build(
        "SART with a total variation denoising step after every iteration.",
        extra=(_TV_STEP.format(tviter=50, tvlambda=50), _SART.format(lmbda_red=0.99)),
        see_also="ossart_tv, asd_pocs",
    ),
    "ossart_tv": _build(
        "OS-SART with a total variation denoising step after every iteration.",
        extra=(_BLOCKSIZE, _TV_STEP.format(tviter=50, tvlambda=50), _SART.format(lmbda_red=0.99)),
        see_also="sart_tv, os_asd_pocs",
    ),
    # Krylov
    "cgls": _build(
        "CGLS: conjugate gradient for least squares.\n\n"
        "Solves :math:`\\min_x \\|Ax - b\\|_2`.",
        extra=(_RESTART,),
        notes=_KRYLOV_NOTE,
        see_also="lsqr, lsmr, hybrid_lsqr",
    ),
    "lsqr": _build(
        "LSQR: least squares by Golub-Kahan bidiagonalisation [1]_.\n\n"
        "Solves :math:`\\min_x \\|Ax - b\\|_2`. Mathematically the same as CGLS, "
        "numerically more stable.",
        extra=(_RESTART,),
        notes=_KRYLOV_NOTE,
        see_also="cgls, lsmr, hybrid_lsqr",
        references=_LSQR_REF,
    ),
    "lsmr": _build(
        "LSMR: least squares by minimising the normal-equation residual [1]_.\n\n"
        "Solves :math:`\\min_x \\|Ax - b\\|_2^2 + \\lambda^2 \\|x\\|_2^2`. Its "
        "residuals decrease more smoothly than LSQR's.",
        extra=(
            "lmbda : float, optional\n    Damping (Tikhonov) parameter :math:`\\lambda`. Default 0: plain least squares.",
            _RESTART,
        ),
        notes=_KRYLOV_NOTE,
        see_also="lsqr, hybrid_lsqr",
        references=_LSMR_REF,
    ),
    "hybrid_lsqr": _build(
        "Hybrid LSQR: LSQR with Tikhonov regularisation of the projected problem.\n\n"
        "Solves :math:`\\min_x \\|Ax - b\\|_2^2 + \\lambda^2 \\|x\\|_2^2` in the "
        "Krylov subspace, choosing :math:`\\lambda` as it goes, which makes it much "
        "less sensitive to the number of iterations than LSQR. Keeps one image "
        "and one set of projections per iteration in memory.",
        extra=(_HYBRID,),
        notes=_KRYLOV_NOTE + "\n\nIt stops early if the residual grows, which signals a loss of orthogonality.",
        see_also="lsqr, hybrid_flsqr_tv",
    ),
    "ab_gmres": _build(
        "AB-GMRES: GMRES with the backprojector as a right preconditioner.\n\n"
        "Designed for projector/backprojector pairs that are not exact "
        "transposes of each other.",
        extra=(
            "backprojector : {\"matched\", \"FDK\"}, optional\n    Backprojection used as the preconditioner: the matched "
            "backprojector (default) or the full FDK algorithm, filtering included.",
        ),
        see_also="ba_gmres, cgls",
    ),
    "ba_gmres": _build(
        "BA-GMRES: GMRES with the backprojector as a left preconditioner.\n\n"
        "Designed for projector/backprojector pairs that are not exact "
        "transposes of each other.",
        extra=(
            "backprojector : {\"matched\", \"FDK\"}, optional\n    Backprojection used as the preconditioner: the matched "
            "backprojector (default) or the full FDK algorithm, filtering included.",
        ),
        see_also="ab_gmres, cgls",
    ),
    "irn_tv_cgls": _build(
        "IRN-TV-CGLS: total variation regularisation by iteratively reweighted norms.\n\n"
        "Approximates :math:`\\min_x \\|Ax - b\\|_2^2 + \\lambda\\, \\mathrm{TV}(x)` by a "
        "sequence of weighted least-squares problems, each solved with CGLS. The "
        "total number of CGLS iterations is ``niter * niter_outer``.",
        extra=(
            "lmbda : float, optional\n    Weight of the TV term :math:`\\lambda`. Default 10.",
            "niter_outer : int, optional\n    Number of reweighting steps. Default 4.",
            _RESTART,
        ),
        see_also="hybrid_flsqr_tv, asd_pocs",
    ),
    "hybrid_flsqr_tv": _build(
        "Hybrid flexible LSQR with total variation regularisation.\n\n"
        "A flexible Krylov method in which the TV regularisation acts as a "
        "preconditioner, with the regularisation parameter chosen as it goes.",
        extra=(_HYBRID,),
        notes="It stops early if the residual grows.",
        see_also="hybrid_lsqr, irn_tv_cgls",
    ),
    # TV / POCS
    "asd_pocs": _build(
        "ASD-POCS: adaptive steepest descent, projection onto convex sets [1]_.\n\n"
        "Minimises the total variation of the image subject to "
        ":math:`\\|Ax - b\\|_2 \\le \\epsilon`. Uses one projection per data step; "
        "``os_asd_pocs`` uses subsets and is faster.",
        extra=(_POCS,),
        notes=_POCS_NOTE,
        see_also="os_asd_pocs, awasd_pocs, pcsd",
        references=_ASD_REF,
    ),
    "os_asd_pocs": _build(
        "OS-ASD-POCS: ASD-POCS with ordered subsets of projections in the data step.",
        extra=(_BLOCKSIZE, _POCS),
        notes=_POCS_NOTE,
        see_also="asd_pocs, os_awasd_pocs",
    ),
    "awasd_pocs": _build(
        "AwASD-POCS: ASD-POCS with adaptive-weighted (edge-preserving) total variation.",
        extra=(_AW, _POCS),
        notes=_POCS_NOTE,
        see_also="asd_pocs, os_awasd_pocs",
    ),
    "os_awasd_pocs": _build(
        "OS-AwASD-POCS: adaptive-weighted ASD-POCS with ordered subsets.",
        extra=(_BLOCKSIZE, _AW, _POCS),
        notes=_POCS_NOTE,
        see_also="awasd_pocs, os_asd_pocs",
    ),
    "pcsd": _build(
        "PCSD: projection-controlled steepest descent.\n\n"
        "A variant of ASD-POCS without the ``alpha``, ``alpha_red`` and ``rmax`` "
        "parameters: the size of the TV step follows the current data error.",
        extra=(_PCSD,),
        notes=_POCS_NOTE,
        see_also="os_pcsd, aw_pcsd, asd_pocs",
    ),
    "os_pcsd": _build(
        "OS-PCSD: PCSD with ordered subsets of projections in the data step.",
        extra=(_BLOCKSIZE, _PCSD),
        notes=_POCS_NOTE,
        see_also="pcsd, os_aw_pcsd",
    ),
    "aw_pcsd": _build(
        "AwPCSD: PCSD with adaptive-weighted (edge-preserving) total variation.",
        extra=(_AW, _PCSD),
        notes=_POCS_NOTE,
        see_also="pcsd, os_aw_pcsd",
    ),
    "os_aw_pcsd": _build(
        "OS-AwPCSD: adaptive-weighted PCSD with ordered subsets.",
        extra=(_BLOCKSIZE, _AW, _PCSD),
        notes=_POCS_NOTE,
        see_also="aw_pcsd, os_pcsd",
    ),
    "piccs": _build(
        "PICCS: prior image constrained compressed sensing.\n\n"
        "ASD-POCS in which the TV term also includes the difference to a prior "
        "image of the same object, such as a planning CT.",
        extra=(_PRIOR, _POCS),
        notes=_POCS_NOTE,
        see_also="os_piccs, asd_pocs",
    ),
    "os_piccs": _build(
        "OS-PICCS: PICCS with ordered subsets of projections in the data step.",
        extra=(_BLOCKSIZE, _PRIOR, _POCS),
        notes=_POCS_NOTE,
        see_also="piccs",
    ),
    # Proximal
    "FISTA": _build(
        "FISTA: fast iterative shrinkage-thresholding with a TV proximal step [1]_.\n\n"
        "Solves :math:`\\min_x \\tfrac{1}{2}\\|Ax - b\\|_2^2 + \\lambda\\, \\mathrm{TV}(x)` "
        "with accelerated proximal gradient steps.",
        extra=(
            "hyper : float or \"auto\", optional\n    Lipschitz constant of the data term's gradient (the largest "
            "eigenvalue of :math:`A^T A`); the gradient step is ``1/hyper``. Default \"auto\": estimated for "
            "this problem by power iteration.",
            "tviter : int, optional\n    Iterations of the TV proximal step. Default 20.",
            "tvlambda : float, optional\n    TV regularisation strength :math:`\\lambda`. Default 0.1.",
            "fista_p, fista_q : float, optional\n    Parameters of the momentum update. Default 1 and 1.",
        ),
        see_also="ista, asd_pocs",
        references=_FISTA_REF,
    ),
    "ISTA": _build(
        "ISTA: iterative shrinkage-thresholding with a TV proximal step.\n\n"
        "FISTA without the acceleration. Takes the same options.",
        extra=(
            "hyper : float or \"auto\", optional\n    Lipschitz constant of the data term's gradient. Default \"auto\".",
            "tviter : int, optional\n    Iterations of the TV proximal step. Default 20.",
            "tvlambda : float, optional\n    TV regularisation strength. Default 0.1.",
        ),
        see_also="fista",
    ),
    # Statistical
    "mlem": _build(
        "MLEM: maximum-likelihood expectation maximisation.\n\n"
        "Models the measurements as Poisson-distributed, which suits low-count "
        "data. The image stays non-negative. Starts from an image of ones unless "
        "``init`` is given.",
        see_also="ossart",
    ),
}
