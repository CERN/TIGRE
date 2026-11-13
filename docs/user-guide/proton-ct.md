# Proton CT

TIGRE includes a MATLAB extension, written in CUDA, that turns proton CT
measurements into optimised proton radiographies, following the method of
[Collins-Fekete et al. (2016)](https://doi.org/10.1088/0031-9155/61/23/8232).
Protons scatter as they cross the object, so their paths are not straight
lines. The extension estimates each proton's path with a cubic spline and
bins its water-equivalent path length (WEPL) along that path.

The radiographies can then be reconstructed with TIGRE's usual algorithms to
give a full proton CT image, or used on their own. A Python version is not
available yet.

## Building the extension

From TIGRE's `MATLAB` folder, run:

```matlab
CompilePCT
```

It asks for a pixel length and an intercepts vector size; press Enter to keep
the defaults. If the pixel lengths in x and y differ, enter the smaller one.
When the code runs, it reports how many protons it rejected because of the
intercepts vector size.

## Input data

Load the measured or simulated data in any format MATLAB reads (`.mat` or raw
binary files work well). All lengths must be in millimetres.

| Parameter | Name in the demo | Unit |
|---|---|---|
| Water-equivalent path length | `WEPL` | mm |
| Entry position, x and y | `posIn` | mm |
| Exit position, x and y | `posOut` | mm |
| Entry direction, x and y | `dirIn` | mm, with the z component equal to 1 |
| Exit direction, x and y | `dirOut` | mm, with the z component equal to 1 |

Before computing the radiography, remove protons with large angular or
energy deviations from the data set (3-sigma cuts), or the result will not be
reliable.

## Computing a radiography

```matlab
proj = pCTCubicSpline_mex(posIn, posOut, dirIn, dirOut, Wepl, eIn, geo);
```

- `Wepl` is a vector of length N with the WEPL of each proton.
- `posIn`, `posOut`, `dirIn` and `dirOut` are vectors of length 2N: first the
  x values of all protons, then the y values.
- `eIn` is the initial beam energy, in MeV.

`geo` uses these fields:

| Field | Meaning | Unit |
|---|---|---|
| `geo.dDetector` | Pixel size of the radiography | mm |
| `geo.sDetector` | Size of the radiography in x and y (two values) | mm |
| `geo.DSO` | Distance from the source to the origin (positive) | mm |
| `geo.DSID` | Distance from the source to the entry detector (positive) | mm |
| `geo.DSD` | Distance from the source to the exit detector (positive) | mm |
| `geo.mode` | `'parallel'` or `'cone'` | |
| `geo.hull` | Convex hull parameters `[a, b, alpha, h]` | mm |

The origin is at the centre of the object. Set `h = 0` to compute without a
convex hull. Otherwise the hull is the elliptical cylinder

$$
\frac{(x\cos\alpha - z\sin\alpha)^2}{a^2} + \frac{(x\sin\alpha + z\cos\alpha)^2}{b^2} = 1,
\qquad -\frac{h}{2} \le y \le \frac{h}{2},
$$

that is, an ellipse rotated by $\alpha$ around the y axis.

## Testing

`MATLAB/Demos/pRad_demo.m` computes one optimised radiography in parallel
geometry from a simulated test data set (loaded with `pCTdata`, in
`MATLAB/Test_data/pCT_test_data`). Run it to check
that the extension works on your system.

## Notes

- If `CompilePCT` fails with missing symbols, check that the `libstdc++`
  shipped with MATLAB is compatible with your CUDA version and your system's
  `gcc`.
- The extension was tested on Debian 10 with CUDA 11.2, an NVIDIA RTX 4000,
  and MATLAB R2020b and R2021b.
