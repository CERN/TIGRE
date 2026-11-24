# Setting up the geometry for your scanner

A step-by-step way to go from a scanner's numbers to a correct
reconstruction, checking each step before the next. It assumes you have read
[](geometry.md), which explains every field. If TIGRE has a
[loader for your scanner](loading-data.md), it does steps 1 to 3 for you: start
at step 4.

## 1. Collect the scanner's numbers

You need, in millimetres and radians:

| Quantity | Field | Where to find it |
|---|---|---|
| Source to rotation axis | `DSO` | Scanner settings file, or the "source to object" distance |
| Source to detector | `DSD` | Scanner settings file, or the "source to detector" distance |
| Detector pixels, and their size | `nDetector`, `dDetector` | Detector specification, after any binning |
| Detector shift, if any | `offDetector` | Settings file; often 0 |
| One angle per projection | `angles` | Settings file, or evenly spaced over the scan |

Distances in some scanners' files are measured to the object or to the
detector housing rather than to the rotation axis and the detector surface.
When in doubt, check the result against an object of known size (step 7).

## 2. Size the image

Choose the voxel size from the detector pixel size at the rotation axis,
`dDetector * DSO / DSD`. Smaller voxels add no detail; larger ones lose some.
Then choose how many voxels you need to cover the object.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```python
import numpy as np
import tigre

geo = tigre.geometry()                          # empty cone-beam geometry
geo.DSD = 1100.0                                # mm
geo.DSO = 700.0                                 # mm
geo.nDetector = np.array([1000, 1200])          # (v, u) pixels
geo.dDetector = np.array([0.2, 0.2])            # (v, u) mm
geo.sDetector = geo.nDetector * geo.dDetector

voxel = geo.dDetector[1] * geo.DSO / geo.DSD    # about 0.127 mm
geo.nVoxel = np.array([900, 1000, 1000])        # (z, y, x)
geo.dVoxel = np.array([voxel, voxel, voxel])
geo.sVoxel = geo.nVoxel * geo.dVoxel

geo.offOrigin = np.array([0, 0, 0])
geo.offDetector = np.array([0, 0])
geo.mode = "cone"
geo.accuracy = 0.5
angles = np.linspace(0, 2 * np.pi, 1440, endpoint=False)
```
:::

:::{tab-item} MATLAB
:sync: matlab

```matlab
geo.DSD = 1100;                                 % mm
geo.DSO = 700;                                  % mm
geo.nDetector = [1200; 1000];                   % [u; v] pixels
geo.dDetector = [0.2; 0.2];                     % [u; v] mm
geo.sDetector = geo.nDetector .* geo.dDetector;

voxel = geo.dDetector(1) * geo.DSO / geo.DSD;   % about 0.127 mm
geo.nVoxel = [1000; 1000; 900];                 % [x; y; z]
geo.dVoxel = [voxel; voxel; voxel];
geo.sVoxel = geo.nVoxel .* geo.dVoxel;

geo.offOrigin = [0; 0; 0];
geo.offDetector = [0; 0];
geo.accuracy = 0.5;
angles = linspace(0, 2*pi, 1441); angles(end) = [];
```
:::
::::

The numbers are an example; use your scanner's.

## 3. Turn the projections into line integrals

TIGRE needs $-\log(I/I_0)$, not raw intensities. Divide by a flat-field image
(after subtracting a dark field, if you have one), take $-\log$, and store
`float32` (`single` in MATLAB). In Python the array must have shape
`(len(angles), *geo.nDetector)`; in MATLAB, `[v, u, angles]`.

## 4. Look at the data and the geometry

- Plot the geometry with `tigre.plot_geometry(geo)` (MATLAB:
  `plotgeometry(geo, 0)`). The image should sit inside the cone of rays, and
  the detector should be large enough to see all of it.
- Step through the projections with `tigre.plotproj(proj)` (MATLAB:
  `plotProj(proj, angles)`). Dense material should be bright and air dark, and
  the object should turn from left to right. If air is bright, step 3 is
  missing; if the object moves up and down, the detector axes are swapped.

## 5. Reconstruct a few slices with FDK

A full reconstruction can take minutes. While checking the geometry,
reconstruct only a few slices around the centre:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```python
import copy
import tigre.algorithms as algs

test = copy.deepcopy(geo)
test.nVoxel = np.array([8, geo.nVoxel[1], geo.nVoxel[2]])
test.sVoxel = test.nVoxel * test.dVoxel
img = algs.fdk(proj, test, angles)
tigre.plotimg(img, dim="z")
```
:::

:::{tab-item} MATLAB
:sync: matlab

```matlab
test = geo;
test.nVoxel(3) = 8;
test.sVoxel = test.nVoxel .* test.dVoxel;
img = FDK(proj, test, angles);
plotImg(img, 'Dim', 'Z');
```
:::
::::

Compare what you see with [](artefacts.md). The usual problems at this stage
are doubled edges (centre of rotation or detector offset), arcs everywhere
(angles reversed), and a bright ring at the edge (truncation).

## 6. Correct the centre of rotation

Most real scanners need a small correction. For a scan over a full circle,
estimate it from the data:

```python
geo.COR = tigre.computeCOR(proj, geo, angles)
# MATLAB: geo.COR = computeCOR(proj, geo, angles);
```

Then repeat step 5. Edges should now be sharp. If they are not, try a few
values of `geo.COR` around the estimate by hand, or check `offDetector`, the
direction of the angles and the detector roll (`rotDetector`); see demos 14
and 17.

## 7. Check the result

- **Size:** measure something of known size in the reconstruction. A wrong
  size means `DSO` or `DSD` is wrong.
- **Consistency:** forward-project your reconstruction with `tigre.Ax(img,
  geo, angles)` and compare it with the measured projections. Apart from
  noise, they should line up; a shift or a scale between them points to the
  field that is wrong.

## 8. Reconstruct

With the geometry checked, reconstruct the whole image with FDK, then try an
iterative algorithm if you need better quality from noisy data or from fewer
projections; see [](algorithms.md). Keep the geometry in a script or a
function, so that every scan from the same scanner uses the same, checked
values.
