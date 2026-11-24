# Roadmap

## Other GPUs

TIGRE runs on NVIDIA GPUs (CUDA). Support for the other main GPU makers is
being added, so that TIGRE can run on most computers with a GPU:

| GPUs | Backend | Status |
|---|---|---|
| AMD | HIP (ROCm) | in progress: [#747](https://github.com/CERN/TIGRE/pull/747) |
| Intel | SYCL | planned: [#885](https://github.com/CERN/TIGRE/issues/885) |
| Apple Silicon (M1 and later) | Metal | planned: [#876](https://github.com/CERN/TIGRE/issues/876) |

The build already has a `TIGRE_BACKEND` setting that chooses between them
(see [Installing TIGRE for Python](install/python.md)); at the moment only
`cuda` can be built. Help with any of these is welcome, especially testing on
hardware the developers do not have: comment on the issue.

## Everything else

Other work follows what users need: bug reports, feature requests and
contributions. To ask for something, open a
[feature request](https://github.com/CERN/TIGRE/issues/new/choose); to see
what is being worked on, look at the
[open pull requests](https://github.com/CERN/TIGRE/pulls).
