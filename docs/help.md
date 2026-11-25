---
og:description: "Where to get help with TIGRE: the guides to check first, how to ask a question, how to report a bug and what to include."
myst:
  html_meta:
    description: "Where to get help with TIGRE: the guides to check first, how to ask a question, how to report a bug and what to include."
---

# Getting help

## Check the guides first

Most questions are answered in one of these:

[](user-guide/troubleshooting.md)
: Installation errors, GPU and memory errors, and common questions.

[](user-guide/artefacts.md)
: Pictures of the usual mistakes (wrong geometry, centre of rotation, data
  not log-transformed, too few projections) and how to fix each one.

[](user-guide/geometry-walkthrough.md)
: Setting up `geo` for your scanner, step by step.

[](user-guide/real-data-tutorial.md)
: A real scan from download to final image.

Search the [issues](https://github.com/CERN/TIGRE/issues?q=is%3Aissue) and
[discussions](https://github.com/CERN/TIGRE/discussions) too: someone may have
asked already.

## Ask a question

For questions about *using* TIGRE (setting up the geometry, choosing an
algorithm, loading your scanner's data, interpreting a result), open a
[Discussion](https://github.com/CERN/TIGRE/discussions/new/choose). Other users
can answer there too, and the answer helps the next person.

## Report a bug

If TIGRE crashes, gives an error it should not, or gives a wrong result, open
an [issue](https://github.com/CERN/TIGRE/issues/new/choose) with the bug
report template. Please include:

- **A short script that shows the problem**, which someone else can run. Use
  `tigre.utilities.sample_loader.load_head_phantom` (Python) or
  `headPhantom` (MATLAB) instead of your own data where you can; if the
  problem needs your data, say how to get a small part of it.
- **The full error message**, as text rather than a screenshot.
- **Your versions.** In Python, run

  ```python
  import tigre
  tigre.show_versions()
  ```

  and in MATLAB, `tigreVersion`. These print TIGRE's version, Python or
  MATLAB, the operating system and your GPUs. (TIGRE 3.3.1 and earlier do not
  have them yet: give these by hand, plus `nvidia-smi`'s first lines.)
- If you compiled TIGRE yourself, the CUDA version (`nvcc --version`) and
  compiler.

Most TIGRE bugs depend on the GPU, the CUDA version or the geometry, so these
details usually decide whether a bug can be fixed.

## Suggest a feature

Open an issue with the feature request template, saying what you want to do
and why. For larger features, describe the design before writing much code:
see [](contributing.md).

## Contact the developers privately

For things that cannot be public (for example a security problem, or data
you cannot share), email
[tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com).
