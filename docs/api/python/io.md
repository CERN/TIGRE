---
og:description: "Python functions that load scans from Nikon, Bruker, YXLON, Diondo and Varian scanners and DXchange files into TIGRE."
myst:
  html_meta:
    description: "Python functions that load scans from Nikon, Bruker, YXLON, Diondo and Varian scanners and DXchange files into TIGRE."
---

# Python scanner data loaders

In `tigre.utilities.io`. Each returns `proj, geo, angles`, ready to
reconstruct. See [](../../user-guide/loading-data.md).

```{eval-rst}
.. currentmodule:: tigre.utilities.io

.. autofunction:: NikonDataLoader

.. autofunction:: BrukerDataLoader

.. autofunction:: YXLONDataLoader

.. autofunction:: DiondoDataLoader

.. autofunction:: VarianDataLoader

.. autofunction:: DXChangeDataLoader
```
