# Software that uses TIGRE

## Libraries

[CIL, the Core Imaging Library](https://github.com/TomographicImaging/CIL)
: A Python framework for tomographic imaging, with many optimisation
  algorithms and regularisers. Its TIGRE plugin
  ([`cil.plugins.tigre`](https://tomographicimaging.github.io/CIL/nightly/plugins))
  uses TIGRE's projectors for CIL's reconstructions. See
  [](user-guide/interoperability.md) for how to use the two together.

[LION, Learned Iterative Optimization Networks](https://github.com/CambridgeCIA/LION)
: A PyTorch library for data-driven (learned) CT reconstruction from the
  University of Cambridge: post-processing networks, unrolled methods, learned
  regularisers and plug-and-play methods, with datasets and benchmark
  experiments. It uses TIGRE for its CT operators.

## Research using TIGRE

TIGRE is used in many research projects. [Papers using TIGRE](papers.md)
lists the published ones.

## Tools TIGRE works with

TIGRE's data are NumPy arrays and MATLAB matrices, so it works with most
scientific software. [](user-guide/interoperability.md) shows how to use it
with PyTorch, CIL and ASTRA.

## Add your software

If your software uses TIGRE, tell us in a
[Discussion](https://github.com/CERN/TIGRE/discussions) or open a pull request
that adds it to this page.
