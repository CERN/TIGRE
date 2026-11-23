# Geometry

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

Geometries with moving sources or detectors
-------------------------------------------

.. autofunction:: staticDetectorGeo

.. autofunction:: staticDetLinearSourceGeo

.. autofunction:: ArbitrarySourceDetMoveGeo

Per-view vectors (ASTRA's ``cone_vec`` layout)
----------------------------------------------

.. autofunction:: from_vector_to_TIGRE

.. autofunction:: from_TIGRE_to_vector

Curved detectors
----------------

.. autofunction:: tigre.utilities.curved_detector.flatten_detector
```
