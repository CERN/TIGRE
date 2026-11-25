---
og:description: "Python utilities in TIGRE: noise simulation, TV denoising, image quality measures, phantoms, GPU selection and version information."
myst:
  html_meta:
    description: "Python utilities in TIGRE: noise simulation, TV denoising, image quality measures, phantoms, GPU selection and version information."
---

# Python utilities

```{eval-rst}
Simulation and test data
------------------------

.. autofunction:: tigre.utilities.CTnoise.add

.. autofunction:: tigre.utilities.sample_loader.load_head_phantom

Denoising
---------

.. autofunction:: tigre.utilities.im_3d_denoise.im3ddenoise

Image quality
-------------

.. autofunction:: tigre.utilities.Measure_Quality.Measure_Quality

GPUs
----

See also :doc:`/user-guide/gpus`.

.. autoclass:: tigre.utilities.gpu.GpuIds
   :no-members:

.. autofunction:: tigre.utilities.gpu.getGpuNames

.. autofunction:: tigre.utilities.gpu.getGpuIds

Version
-------

``tigre.__version__`` is TIGRE's version as text, for example ``"3.3.1"``.

.. autofunction:: tigre.show_versions
```
