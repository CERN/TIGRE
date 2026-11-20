# Choosing an algorithm

TIGRE has around 30 reconstruction algorithms. They all take the same four
inputs (projections, geometry, angles and, for iterative ones, the number of
iterations), so trying several on the same data is cheap.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
import tigre.algorithms as algs

img_fdk = algs.fdk(proj, geo, angles)
img_sart = algs.ossart(proj, geo, angles, 20)
img_cgls = algs.cgls(proj, geo, angles, 15)
img_tv = algs.asd_pocs(proj, geo, angles, 20)
```
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
img_fdk = FDK(proj, geo, angles);
img_sart = OS_SART(proj, geo, angles, 20);
img_cgls = CGLS(proj, geo, angles, 15);
img_tv = ASD_POCS(proj, geo, angles, 20);
```
:::
::::

## Where to start

**Plenty of good-quality projections, circular scan:** start with **FDK**
(cone beam) or **FBP** (parallel beam). They take one pass and are usually
the right answer for well-sampled data.

**Fewer projections or noisier data:** try **OS-SART** or **CGLS**. OS-SART is
robust and easy to use. CGLS reaches a low data error in very few iterations
(often 10 to 30), but it starts to fit the noise if you run it for too long,
so stop early.

**Very few projections, limited angle or very noisy data:** use a **total
variation (TV)** method such as **OS-ASD-POCS**, **FISTA** or **OS-SART-TV**.
They assume the object is made of regions of nearly constant value, which
suits most real objects, and can produce good images from poor data. They are
slower and their parameters need tuning.

**Scans that are not a circle** (arbitrary rotation axes, some helical or
moving-source geometries): FDK assumes a circular trajectory, so use an
iterative algorithm.

**Low-count data:** **MLEM** (and **OSEM** in MATLAB) model Poisson noise in
the measurements.

**A previous image of the same object is available** (for example a planning
CT): **PICCS** uses it as a prior.

Then compare: the algorithm demos (6 to 9 and 20) reconstruct the same
phantom with each family, so you can see the trade-offs.

## Algorithm families

### Filtered backprojection

| Algorithm | Python | MATLAB |
|---|---|---|
| FDK, cone beam | `fdk` | `FDK` |
| FBP, parallel beam | `fbp` | `FBP` |

Both accept a `filter` option: `ram_lak` (default), `shepp_logan`,
`cosine`, `hamming` or `hann`. Smoother filters reduce noise at the cost of
resolution. On short scans (less than a full circle), FDK applies Parker
weights so that rays measured twice are not counted twice; turn them off with
`parker=False` (`'parker', false` in MATLAB). When the detector is shifted
(`offDetector`, as in half-fan scans), FDK applies Wang weights on
full-circle scans and skips them on short scans, where they would leave
shading (demo 26).

### SART family

Algebraic methods that update the image using the error in the projections.
They differ in how many projections they use per update:

| Algorithm | Python | MATLAB | Projections per update |
|---|---|---|---|
| SART | `sart` | `SART` | 1 |
| OS-SART | `ossart` | `OS_SART` | `blocksize` (default 20) |
| SIRT | `sirt` | `SIRT` | all |

SART converges in the fewest iterations but each iteration is slowest. SIRT
is the opposite. OS-SART sits in between and is usually the best choice.

### Krylov subspace methods

Fast least-squares solvers. They use the `matched` backprojector, see
[](projectors.md).

| Algorithm | Python | MATLAB |
|---|---|---|
| CGLS | `cgls` | `CGLS` |
| LSQR | `lsqr` | `LSQR` |
| LSMR | `lsmr` | `LSMR` |
| Hybrid LSQR | `hybrid_lsqr` | `hybrid_LSQR` |
| AB-GMRES, BA-GMRES | `ab_gmres`, `ba_gmres` | `AB_GMRES`, `BA_GMRES` |
| IRN-TV-CGLS | `irn_tv_cgls` | `IRN_TV_CGLS` |
| Hybrid flexible LSQR with TV | `hybrid_flsqr_tv` | `hybrid_fLSQR_TV` |

CGLS, LSQR and LSMR need early stopping. Hybrid LSQR adds Tikhonov
regularisation and chooses its strength at every iteration (by generalised
cross validation, or by the discrepancy principle if you give the noise level
with `noiselevel` / `'Noiselevel'`), which makes it much less sensitive to the
number of iterations; set `lmbda` / `'lambda'` to fix it yourself. AB-GMRES
and BA-GMRES are designed for projector pairs that are not exact transposes.
The last two add TV regularisation.

### Total variation

These alternate a data step (SART-like) with a step that lowers the image's
total variation.

| Algorithm | Python | MATLAB |
|---|---|---|
| ASD-POCS | `asd_pocs`, `os_asd_pocs` | `ASD_POCS`, `OS_ASD_POCS` |
| Adaptive-weighted ASD-POCS (edge preserving) | `awasd_pocs`, `os_awasd_pocs` | `AwASD_POCS`, `OS_AwASD_POCS` |
| PCSD (chooses some ASD-POCS parameters for you) | `pcsd`, `os_pcsd` | `PCSD`, `OS_PCSD` |
| Adaptive-weighted PCSD | `aw_pcsd`, `os_aw_pcsd` | `AwPCSD`, `OS_AwPCSD` |
| SART-TV | `sart_tv`, `ossart_tv` | `SART_TV` |
| FISTA, ISTA with a TV proximal step | `fista`, `ista` | `FISTA` |

The `os_` versions use subsets of projections, like OS-SART, and are faster.
The options that matter most:

- ASD-POCS: `maxl2err` (called epsilon in the paper; the data error the
  image is allowed to keep), `alpha` (TV step size) and `tviter` (TV steps
  per iteration). See demo 9.
- FISTA: `hyper` (an estimate of the largest eigenvalue of $A^T A$),
  `tviter` and `tvlambda`. See demo 20.

### Prior image

When an earlier image of the same object is available (a planning CT, or a
scan from a previous visit), prior image constrained compressed sensing
(PICCS) uses it to reconstruct from far fewer or noisier projections. It works
like ASD-POCS, but the regularisation is a mix of the TV of the image and the
TV of its difference from the prior, so the reconstruction keeps the prior's
structure where the data agree with it and departs from it where they do not.

| Algorithm | Python | MATLAB |
|---|---|---|
| PICCS | `piccs`, `os_piccs` | `PICCS` |

Pass the prior image with `prior=` in Python, or as the fifth argument in
MATLAB (`PICCS(proj, geo, angles, niter, prior)`); `prior_ratio` (default
0.2) weights the prior term against the plain TV term. The prior must be
aligned with the current scan. See MATLAB demo 24.

### Statistical

| Algorithm | Python | MATLAB |
|---|---|---|
| MLEM | `mlem` | `MLEM` |
| OSEM | | `OSEM` |

## Options shared by the iterative algorithms

Most iterative algorithms accept these options, as keyword arguments in
Python and name-value pairs in MATLAB. Not every algorithm uses every
option; each algorithm's help text (`help(algs.ossart)` or `help OS_SART`)
lists what it accepts.

| Purpose | Python | MATLAB |
|---|---|---|
| Projections per subset | `blocksize` | `'BlockSize'` |
| Step size, and its reduction per iteration | `lmbda`, `lmbda_red` | `'lambda'`, `'lambda_red'` |
| Starting image | `init`: `None`, `"FDK"`, `"multigrid"` or an array | `'init'`: `'none'`, `'FDK'`, `'multigrid'` or `'image'` with `'initimg'` |
| Order of the subsets | `OrderStrategy`: `"ordered"`, `"random"`, `"angularDistance"` | `'OrderStrategy'` |
| Keep values non-negative | `noneg` (default `True`) | `'nonneg'` |
| Quality measures per iteration | `Quameasopts`, e.g. `["RMSE", "CC"]` | `'QualMeas'` |
| Choose GPUs | `gpuids` | `'gpuids'` |
| Progress messages | `verbose` | `'verbose'` |

In Python, an algorithm returns the image. With `computel2=True` or with
`Quameasopts` it returns a tuple `(image, errors)`. In MATLAB, the extra
outputs (for example `[img, resL2] = CGLS(...)`) return the residual norm
and quality measures.

Unknown keyword arguments in Python print a warning rather than raising an
error, so check the output for typos in option names.
