---
og:description: "Measurements of the speed of TIGRE's projectors and backprojectors in version 3.3.1, on one and several GPUs."
myst:
  html_meta:
    description: "Measurements of the speed of TIGRE's projectors and backprojectors in version 3.3.1, on one and several GPUs."
---

# Projector speed in TIGRE 3.3.1

*Summary of the technical note "TIGRE 3.3.1 and ASTRA 2.5: a speed comparison
of the projection operators on one and several GPUs", Ander Biguri, November
2026.*
**{download}`Download the full note (PDF) <tigre331_vs_astra.pdf>`**

## What was measured

The four projection operators of TIGRE 3.3.1 (the Siddon and interpolated
forward projectors, and the matched and FDK-weighted backprojectors) were
timed against the ASTRA Toolbox 2.5, the other open GPU toolbox most used
from Python:

- on a desktop GPU (RTX 4070, Windows) and on a server with four Quadro RTX
  6000 (Linux);
- for cone-beam, parallel-beam and fan-beam (2D) geometries, from 32³ to
  1536³;
- with data in host memory (NumPy) and already on the GPU (PyTorch).

Both toolboxes compute the same projections of the same image (they agree to
5×10⁻⁴), and every result is a ratio of times with a confidence interval,
measured with the two toolboxes alternating.

## Main results

- **One GPU (RTX 4070):** TIGRE was faster in 129 of 205 comparisons and
  ASTRA in 64. Leaving out TIGRE's interpolated projector, which this release
  did not optimise, TIGRE was faster in 127 of 151 and ASTRA in 13.
- **3D:** the Siddon projector takes about half of ASTRA's time, and the
  backprojectors 0.54–0.67 of it.
- **2D (fan beam):** the backprojectors take about a fifth of ASTRA's time;
  the Siddon projector 1.0–1.3 times, because ASTRA's interpolating projector
  does less work per ray than exact path lengths.
- **Large problems:** forward projecting a 1536³ image with 1536 angles takes
  20.5 s on one GPU and 10.4 s on four, against 83.8 s and 24.2 s for ASTRA.
- **Several GPUs, problems that fit in one:** TIGRE gains less than ASTRA from
  extra GPUs, because each multi-GPU call copies the whole input to every
  GPU. On the server measured, one or two GPUs were the fastest choice for
  such problems. Reducing this is the next step.

The note explains every case where ASTRA is faster.

## What changed in 3.3.1

Most of the work removed costs that do not grow with the problem: allocating
memory, waiting for the GPU, page-locking host memory and copying the
geometry from Python. Those costs decide the speed of small and 2D problems,
of iterative algorithms, which call the operators hundreds of times, and of
operators inside a neural network. The largest single gain was in the
backprojectors: reading the projections from a layered 2D texture instead of
a 3D one halves the cost of each texture read, with identical results. The
Siddon projector now takes several steps of its ray walk per pass (see
[](../siddon/index.md)). Results are bit-identical to 3.3.0, except where a
sum is now ordered differently (relative differences of 10⁻⁶ or less).

The full list of changes, by pull request, is in the note and in the
[release notes](../../changelog.md).

## Run the comparison on your machine

The scripts are in
[`docs/technical-notes/speed-3.3.1`](https://github.com/CERN/TIGRE/tree/master/docs/technical-notes/speed-3.3.1).
You need TIGRE 3.3.1 or later, ASTRA 2.5, PyTorch with CUDA and NumPy in one
environment. From outside TIGRE's source tree:

```bash
python tigre_vs_astra_suite.py --check
python tigre_vs_astra_suite.py --ops fp bp --gpus 1 --label my-gpu --out suite.json --csv suite.csv
```

It takes about an hour on one GPU; add `--gpus 1 2 4` to measure several.
The thresholds TIGRE uses to choose how to run come from two GPU models, so
**results from other GPUs and multi-GPU machines are very welcome**: please
open an [issue](https://github.com/CERN/TIGRE/issues/new) with your
`suite.json` and a line about the machine. The raw results of the note are
attached to the
[3.3.1 release](https://github.com/CERN/TIGRE/releases/tag/v3.3.1).
