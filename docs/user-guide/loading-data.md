# Loading scanner data

TIGRE reconstructs from **line integrals**: each projection value is
$-\log(I / I_0)$, where $I$ is the measured intensity and $I_0$ the intensity
with nothing in the beam. The scanner loaders below read a dataset, convert
it to line integrals and build the matching geometry and angles, so you can
reconstruct straight away.

## Supported scanners

| Scanner | Python (`tigre.utilities.io`) | MATLAB |
|---|---|---|
| Nikon (XTek) micro-CT | `NikonDataLoader` | `NikonDataLoader` |
| Bruker SkyScan | `BrukerDataLoader` | `BrukerDataLoader` |
| YXLON | `YXLONDataLoader` | `YXLONDataLoader` |
| Diondo | `DiondoDataLoader` | `DiondoDataLoader` |
| Varian TrueBeam on-board CBCT | `VarianDataLoader` | `VarianDataLoader` |
| DXchange (HDF5, used at synchrotrons) | `DXChangeDataLoader` | `DXChangeDataLoader` |
| DICOM (tested on Philips Allura only) | | `dicomDataLoader` |

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
import tigre.utilities.io as tigreio
import tigre.algorithms as algs

proj, geo, angles = tigreio.NikonDataLoader("path/to/Nikon/sample/")
img = algs.fdk(proj, geo, angles)
```
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
[proj, geo, angles] = NikonDataLoader('path/to/Nikon/sample/');
img = FDK(proj, geo, angles);
```
:::
::::

### Loading part of a dataset

Micro-CT datasets can have thousands of projections. The Nikon, Bruker,
YXLON and Diondo loaders can load a subset:

| Option | Effect |
|---|---|
| `sampling="equidistant"` with `num_angles=N` | N projections spread evenly over the scan (the default sampling). |
| `sampling="step"` with `sampling_step=k` | Every k-th projection. |
| `sampling="continuous"` with `num_angles=N` | The first N projections. |

In MATLAB these are name-value pairs, for example
`NikonDataLoader(folder, 'sampling', 'step', 'sampling_step', 10)`.

Bruker folders sometimes hold more than one dataset; choose one with
`dataset_number` (a number from 0, or `"all"` in Python).

### Varian on-board CBCT

`VarianDataLoader` reads TrueBeam `.xim` projections and calibration files.
By default it also applies an acceleration-deceleration correction (`acdc`),
detector point-scatter correction (`dps`) and kernel-based scatter
correction (`sc`); each can be switched off with the option of the same name.
It has been tested on clinical data from TrueBeam 2.7. The MATLAB version
has further options, including beam-hardening correction.

## Other scanners

If your scanner is not listed, or it does not save its geometry (some only
store the projections), load the data yourself and describe the geometry by
hand ([](geometry.md)). Then check the data before reconstructing:

1. **Convert to line integrals.** Plot the projections with `plotproj`
   (Python) or `plotProj` (MATLAB). Dense material should be bright and air
   dark. If it is the other way round, the data are still intensities:
   convert them with $-\log(I / I_0)$, using a flat-field image or the
   brightest air value as $I_0$.
2. **Check the orientation.** Python projections have shape
   `(angles, v, u)`, where `v` runs along the rotation axis. The object
   should rotate left to right as you step through the angles; if it moves
   up and down, swap the detector axes.
3. **Reconstruct with FDK first.** If the image looks wrong, check whether
   the angles need to go the other way, whether the geometry matches the
   scanner and, for micro-CT, whether halos around features show that the
   centre of rotation needs correcting (`COR`, or `computeCOR`).

If you write a loader for a scanner TIGRE does not support yet and can share
it, please [open a pull request](https://github.com/CERN/TIGRE/pulls).
