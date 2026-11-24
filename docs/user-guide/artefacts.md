# Why does my reconstruction look wrong?

When a reconstruction looks wrong, the cause is usually the data or the
geometry, not the algorithm. Each picture below shows a common mistake made on
purpose: the head phantom is projected with one thing wrong and reconstructed
with FDK. Find the one that looks like your image, then apply the fix.

| What you see | Likely cause |
|---|---|
| Edges doubled or smeared into arcs, worse further from the centre | [Centre of rotation](#centre-of-rotation) or [detector offset](#detector-offset) |
| The whole image blurred into arcs, or rotated | [Angles](#angles) in the wrong direction or starting elsewhere |
| Inverted contrast, a dark object in a bright disc, or strong cupping | [Missing log or flat field](#log-and-flat-field) |
| Rings and swirls, no recognisable object | [Detector axes swapped](#detector-axes-swapped) |
| A bright ring at the edge of the field of view | [Truncation](#truncation): the object is wider than the detector sees |
| Bands of shading across the image on a scan of about 180° | [Short scan](#short-scan) without Parker weights |
| The image looks fine, but sizes are wrong | [Distances](#distances) `DSO` or `DSD` |
| Blur that grows towards the top and bottom slices | [Detector roll](#detector-roll) |
| Straight streaks radiating across the image | [Too few projections](#too-few-projections) |
| Blurred, or very noisy, from an iterative algorithm | [Number of iterations](#number-of-iterations) |

The images are made by `docs/tools/make_artefacts.py`, so you can reproduce
any of them and try the fixes yourself.

## Centre of rotation

```{image} images/artefacts/cor.png
:alt: Head slice reconstructed correctly, and with the centre of rotation 1.5 mm and 4 mm off
```

**What you see:** edges that are doubled, or smeared into arcs, more so the
further they are from the centre of the image. Even 1.5 mm, about one voxel
here, softens the edges.

**Cause:** the rotation axis does not project onto the centre of the detector,
for example because the object stage is not perfectly aligned. This is the most
common problem with micro-CT data.

**Fix:** estimate the shift from the data and set `geo.COR`:

```python
geo.COR = tigre.computeCOR(proj, geo, angles)
```

In MATLAB, `geo.COR = computeCOR(proj, geo, angles);`. `computeCOR` needs a scan
over a full circle. Otherwise, try a few values of `geo.COR` on a few slices
(see [](geometry-walkthrough.md)) and keep the sharpest.

## Detector offset

```{image} images/artefacts/offset.png
:alt: Head slice reconstructed correctly, with a 6 mm detector shift ignored, and with the shift set in the geometry
```

**What you see:** the same doubled edges as a wrong centre of rotation.

**Cause:** the detector is shifted sideways (for example, on purpose to widen
the field of view), and the geometry does not say so.

**Fix:** set `geo.offDetector` to the shift in mm, `(v, u)` in Python and
`[u; v]` in MATLAB. The right-hand image is the same data with the offset set.
For a detector shifted by more than a few centimetres, see the Wang weights in
[](algorithms.md).

## Angles

```{image} images/artefacts/angles.png
:alt: Head slice reconstructed correctly, with the angles reversed, and with the angles starting 30 degrees off
```

**What you see:** with the angles going the wrong way, everything is blurred
into arcs; with the right direction but a different starting angle, a sharp
image that is rotated.

**Cause:** the scanner's angles increase in the opposite direction to
TIGRE's (see [the coordinate system](geometry.md#coordinate-system)), or its
first projection is not at TIGRE's angle 0, where the source is on the +x
axis.

**Fix:** reverse the angles (`angles = -angles`), or add the starting angle.
A rotated image is often harmless, but a rotation that differs from the one
your scanner's software produces means the angles need checking.

## Log and flat field

```{image} images/artefacts/log.png
:alt: Head slice reconstructed correctly, from raw intensities, and from the log of intensities that were not normalised
```

**What you see:** from raw intensities, a dark object inside a bright disc,
nothing like the real image. With the log taken but no flat-field
normalisation, the object is recognisable but sits on a strong bowl of
brightness (cupping) with a dark ring around it.

**Cause:** TIGRE reconstructs from line integrals, $-\log(I/I_0)$, where $I_0$
is the intensity with no object (the flat field). Intensities, or the log
without dividing by $I_0$, are not line integrals.

**Fix:** divide by a flat-field image (subtracting a dark-field image first,
if your scanner gives one), then take $-\log$. Air should be close to 0 and
dense material bright. The scanner loaders in [](loading-data.md) do this for
you.

## Detector axes swapped

```{image} images/artefacts/axes.png
:alt: Head slice reconstructed correctly, and from projections whose two detector axes were swapped
```

**What you see:** rings and swirls with no recognisable object.

**Cause:** the projections were loaded with their two detector axes swapped,
or transposed. In Python, TIGRE expects projections of shape
`(angles, v, u)`, where `v` runs along the rotation axis; in MATLAB,
`(v, u, angles)`. Images saved by other software are often in the other
order.

**Fix:** step through the projections with `plotproj` (Python) or `plotProj`
(MATLAB): the object should turn from left to right. If it moves up and down,
swap the detector axes, for example with `np.transpose(proj, (0, 2, 1))`.

## Truncation

```{image} images/artefacts/truncation.png
:alt: Head slice reconstructed correctly, and from a detector too narrow to see the whole object
```

**What you see:** a bright ring at the edge of the field of view, and wrong
values near it.

**Cause:** the object is wider than the part of it the detector sees, so some
rays through the object are missing.

**Fix:** if you can, scan again with the whole object in view, or with an
offset detector. Otherwise, iterative algorithms and a reconstruction region
limited to what the detector sees reduce the ring. Check with `plot_geometry`
(Python) or `plotgeometry` (MATLAB) that the image fits inside the cone of
rays.

## Short scan

```{image} images/artefacts/shortscan.png
:alt: Head slice from a short scan without Parker weights, its difference from a full scan, and with Parker weights
```

**What you see:** bands of shading across the image, along the directions
where the scan starts and ends. The middle panel shows the difference from a
full scan.

**Cause:** in a short scan (180° plus the fan angle) some rays are measured
twice and others once. FDK must weight them (Parker weights) to count each
ray once.

**Fix:** TIGRE applies Parker weights when it detects a short scan. If your
angles are not exactly a short scan, or you turned the weights off, pass
`parker=True` to `fdk` (in MATLAB, `'parker', true`). A short scan must cover
at least 180° plus the fan angle, $2\arctan(s_u / (2\,\mathrm{DSD}))$, where
$s_u$ is the width of the detector.

## Distances

```{image} images/artefacts/distances.png
:alt: Head slice reconstructed correctly, with DSO 10 percent too large, and their difference
```

**What you see:** almost nothing: the image looks right. But the object is the
wrong size, and edges are slightly softer. The difference image shows it.

**Cause:** `DSO` or `DSD` does not match the scanner. They set the
magnification, so a wrong value scales the reconstruction.

**Fix:** take the distances from the scanner's files (the loaders read them)
rather than from a drawing, and check the result against something of known
size. Errors of a few percent are invisible by eye but matter for
measurements.

## Detector roll

```{image} images/artefacts/roll.png
:alt: Coronal head slice reconstructed correctly, with the detector rolled 3 degrees and not corrected, and corrected
```

**What you see:** in a vertical (coronal) slice, blur that grows towards the
top and bottom of the image, while the centre stays sharp.

**Cause:** the detector is rotated slightly within its own plane.

**Fix:** set the roll in `geo.rotDetector` (roll, pitch, yaw, in radians, in
both languages). The right-hand image is the same data with the roll set.
See demo 17.

## Too few projections

```{image} images/artefacts/fewangles.png
:alt: Head slice reconstructed with FDK from 360 projections and from 40 projections
```

**What you see:** straight streaks radiating across the image.

**Cause:** FDK needs many projections, roughly as many as the detector has
pixels across. With fewer, each one leaves a streak.

**Fix:** take more projections if you can. Otherwise use an iterative
algorithm, preferably one with total-variation regularisation, which is
designed for few projections; see [](algorithms.md).

## Number of iterations

```{image} images/artefacts/iterations.png
:alt: Head slice reconstructed with CGLS from noisy data after 3, 20 and 150 iterations
```

**What you see:** after a few iterations, a blurred image; after many, a sharp
but very noisy one.

**Cause:** iterative algorithms recover the coarse image first and the fine
detail later, and with noisy data the "fine detail" they add at the end is
noise. This is called semi-convergence.

**Fix:** stop earlier, or use a regularised algorithm (total variation, or the
hybrid Krylov methods, which choose their regularisation automatically). Look
at a few iteration counts rather than trusting one; see [](algorithms.md).

## Still wrong?

Check the geometry step by step with [](geometry-walkthrough.md), and see
[](troubleshooting.md). If none of this explains it, ask in the
[discussions](https://github.com/CERN/TIGRE/discussions), with a picture of
the reconstruction and your geometry.
