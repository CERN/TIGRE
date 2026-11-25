---
og:description: "Technical notes on how parts of TIGRE work inside, such as its Siddon ray tracing and the speed of its projectors."
myst:
  html_meta:
    description: "Technical notes on how parts of TIGRE work inside, such as its Siddon ray tracing and the speed of its projectors."
---

# Technical notes

Short write-ups of how parts of TIGRE work inside, and of the measurements
behind changes to them. Each page here summarises a note; the full note, with
the derivations, tables and references, is a PDF you can download from its
page.

[](siddon/index.md)
: How TIGRE's Siddon forward projector traces rays on the GPU, two flaws of
  the usual implementation and how the projector avoids them, and why it is
  now faster. November 2026, revised for TIGRE 3.3.1.

[](speed-3.3.1/index.md)
: The speed of TIGRE 3.3.1's projectors and backprojectors on one and
  several GPUs, from 2D to 1536³ problems, measured against the ASTRA Toolbox,
  and the changes in 3.3.1 that made them faster. November 2026.

```{toctree}
:hidden:

siddon/index
speed-3.3.1/index
```
