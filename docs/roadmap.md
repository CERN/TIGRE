---
og:description: "What is planned for TIGRE: support for AMD, Intel and Apple Silicon GPUs, and beam-hardening and scatter correction."
myst:
  html_meta:
    description: "What is planned for TIGRE: support for AMD, Intel and Apple Silicon GPUs, and beam-hardening and scatter correction."
---

# Roadmap

## Other GPUs

TIGRE runs on NVIDIA GPUs (CUDA), and from Python also on AMD GPUs (HIP).
Support for the other main GPU makers is being added, so that TIGRE can run on
most computers with a GPU:

| GPUs | Backend | Status |
|---|---|---|
| AMD | HIP (ROCm) | in progress: [#747](https://github.com/CERN/TIGRE/pull/747) |
| Intel | SYCL | planned: [#885](https://github.com/CERN/TIGRE/issues/885) |
| Apple Silicon (M1 and later) | Metal | planned: [#876](https://github.com/CERN/TIGRE/issues/876) |

The build already has a `TIGRE_BACKEND` setting that chooses between them
(see [Installing TIGRE for Python](install/python.md)); at the moment `cuda`
and `hip` can be built. Help with any of these is welcome, especially testing
on hardware the developers do not have: comment on the issue.

## Corrections for real data

Real X-ray data are not exactly the line integrals that reconstruction
assumes, and two effects cause most of the difference:

- **beam hardening**: the X-ray beam has a range of energies, and the
  low-energy X-rays are absorbed first, which gives cupping and dark streaks
  between dense objects;
- **scattered X-rays** reaching the detector, which reduce contrast and add
  streaks.

We want TIGRE to include correction methods for both. Which methods is still
open: if you have one you would like to see in TIGRE, or would like to
contribute one, open a
[Discussion](https://github.com/CERN/TIGRE/discussions).

## Everything else

Other work follows what users need: bug reports, feature requests and
contributions. To ask for something, open a
[feature request](https://github.com/CERN/TIGRE/issues/new/choose); to see
what is being worked on, look at the
[open pull requests](https://github.com/CERN/TIGRE/pulls).
