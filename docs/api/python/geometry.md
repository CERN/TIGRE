---
og:description: "Python classes and functions to create, check and plot a TIGRE scanner geometry, and to find the centre of rotation."
myst:
  html_meta:
    description: "Python classes and functions to create, check and plot a TIGRE scanner geometry, and to find the centre of rotation."
---

# Python geometry

How to describe a scanner: see [](../../user-guide/geometry.md).

```{eval-rst}
.. currentmodule:: tigre

.. autofunction:: geometry

.. autoclass:: tigre.utilities.geometry.Geometry
   :no-members:

.. autoclass:: geometry_default
   :no-members:

.. autoclass:: fan_geometry_default
   :no-members:

.. autofunction:: computeCOR

.. autofunction:: tiltedAxisGeo

.. autofunction:: computeImageGeo

.. autofunction:: changeResolution

Geometries with moving sources or detectors
-------------------------------------------

.. autofunction:: staticDetectorGeo

.. autofunction:: staticDetLinearSourceGeo

.. autofunction:: ArbitrarySourceDetMoveGeo

Curved detectors
----------------

.. autofunction:: tigre.utilities.curved_detector.flatten_detector
```
