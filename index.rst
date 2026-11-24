:orphan:

Demos
=====

Worked examples, generated from TIGRE's demos. Each page shows the Python demo and its output; open the drop-down at the top of a page for the MATLAB version. The demos themselves are in `Python/demos <https://github.com/CERN/TIGRE/tree/master/Python/demos>`__ and `MATLAB/Demos <https://github.com/CERN/TIGRE/tree/master/MATLAB/Demos>`__.


.. raw:: html

  <div id='sg-tag-list' class='sphx-glr-tag-list'></div>


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. thumbnail-parent-div-close

.. raw:: html

    </div>

Getting started
---------------

Describe a scanner, simulate data, reconstruct and look at the result.


.. raw:: html

  <div id='sg-tag-list' class='sphx-glr-tag-list'></div>


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Create a cone-beam geometry field by field, and the default geometries that save typing. See /user-guide/geometry.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d01_CreateGeometry_thumb.png
    :alt:

  :doc:`/demos/01_basics/d01_CreateGeometry`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Describe the scanner geometry</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Load the head phantom and the other images that ship with TIGRE.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d02_SampleData_thumb.gif
    :alt:

  :doc:`/demos/01_basics/d02_SampleData`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Load the sample data</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Forward-project an image and add realistic Poisson and Gaussian noise.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d03_generateData_thumb.gif
    :alt:

  :doc:`/demos/01_basics/d03_generateData`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Simulate projections and noise</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Reconstruct simulated data with FDK and OS-SART and compare the results.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d04_SimpleReconstruction_thumb.gif
    :alt:

  :doc:`/demos/01_basics/d04_SimpleReconstruction`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">A first reconstruction</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="All the options of plotimg and plotproj.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d05_plottingFunctions_thumb.gif
    :alt:

  :doc:`/demos/01_basics/d05_plottingFunctions`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Plot images and projections</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Compute error and similarity measures against a reference image, also during the iterations of an algorithm.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d10_MeasureQuality_thumb.png
    :alt:

  :doc:`/demos/01_basics/d10_MeasureQuality`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Measure image quality</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Denoise a reconstruction with total variation after it is computed.">

.. only:: html

  .. image:: /demos/01_basics/images/thumb/sphx_glr_d11_PostProcessing_thumb.gif
    :alt:

  :doc:`/demos/01_basics/d11_PostProcessing`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Post-process a reconstruction</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Algorithms
----------

Each family of reconstruction algorithms on the same data. See :doc:`/user-guide/algorithms` for how to choose.


.. raw:: html

  <div id='sg-tag-list' class='sphx-glr-tag-list'></div>


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Filtered backprojection for cone beam, with the available filters.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d06_Algorithms01_thumb.gif
    :alt:

  :doc:`/demos/02_algorithms/d06_Algorithms01`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">FDK and its filters</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="The SART family, the options they share, and how they compare.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d07_Algorithms02_thumb.png
    :alt:

  :doc:`/demos/02_algorithms/d07_Algorithms02`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">SART, OS-SART and SIRT</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="CGLS, LSQR, LSMR, hybrid LSQR and the GMRES variants.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d08_Algorithms03_thumb.gif
    :alt:

  :doc:`/demos/02_algorithms/d08_Algorithms03`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Krylov subspace methods</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="ASD-POCS and its variants, and their parameters.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d09_Algorithms04_thumb.gif
    :alt:

  :doc:`/demos/02_algorithms/d09_Algorithms04`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Total variation algorithms</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="FISTA, which reaches a TV-regularised solution faster than gradient methods.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d20_Algorithms05_thumb.gif
    :alt:

  :doc:`/demos/02_algorithms/d20_Algorithms05`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">FISTA with a TV proximal step</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Use a prior image of the same object, such as a planning CT.">

.. only:: html

  .. image:: /demos/02_algorithms/images/thumb/sphx_glr_d24_Algorithms06_thumb.png
    :alt:

  :doc:`/demos/02_algorithms/d24_Algorithms06`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Prior images with PICCS</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Geometries and operators
------------------------

Projection operators and the scanner geometries TIGRE can describe.


.. raw:: html

  <div id='sg-tag-list' class='sphx-glr-tag-list'></div>


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Siddon and interpolated projection, FDK-weighted and matched backprojection. See /user-guide/projectors.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d12_ProjectionOperations_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d12_ProjectionOperations`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Projectors and backprojectors</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Move the image along the rotation axis at each angle.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d13_HelicalGeometry_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d13_HelicalGeometry`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Helical scans</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Shift the image and the detector, per scan or per angle.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d14_Offsets_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d14_Offsets`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Image and detector offsets</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Reconstruct 3D parallel-beam data, as from a synchrotron.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d15_3Dparallel_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d15_3Dparallel`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Parallel-beam geometry</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Fan-beam and parallel-beam 2D problems, as one-slice 3D problems.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d16_2Dtomography_thumb.png
    :alt:

  :doc:`/demos/03_geometry/d16_2Dtomography`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">2D tomography</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Describe a detector that is rotated or tilted.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d17_DetectorRotation_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d17_DetectorRotation`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Detector rotation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Trajectories that are not a circle around z, given as Euler angles.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d18_ArbitraryAxisOfRotation_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d18_ArbitraryAxisOfRotation`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Arbitrary axis of rotation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="A limited-angle cone-beam scan, reconstructed with FDK and OS-SART.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d19_DBTexample_thumb.gif
    :alt:

  :doc:`/demos/03_geometry/d19_DBTexample`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Digital breast tomosynthesis</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Resample curved-detector data onto a flat detector, then reconstruct.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d24_CurvedDetector_thumb.png
    :alt:

  :doc:`/demos/03_geometry/d24_CurvedDetector`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Curved detectors</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Why Wang weights are skipped on short scans, shown on a uniform cylinder.">

.. only:: html

  .. image:: /demos/03_geometry/images/thumb/sphx_glr_d26_ShortScanWangWeights_thumb.png
    :alt:

  :doc:`/demos/03_geometry/d26_ShortScanWangWeights`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Shifted detector on a short scan</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Real data, GPUs and other tools
-------------------------------

Loading scanner data, choosing GPUs and using TIGRE from other software. These demos need data or hardware that the website build does not have, so they are shown without output.


.. raw:: html

  <div id='sg-tag-list' class='sphx-glr-tag-list'></div>


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Read Nikon, Bruker, YXLON, Diondo and Varian datasets. See /user-guide/loading-data.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_d21_ScannerDataLoader_thumb.png
    :alt:

  :doc:`/demos/04_data/d21_ScannerDataLoader`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Load scanner data</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Find the installed GPUs and choose which ones TIGRE uses. See /user-guide/gpus.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_d22_ListGPUs_thumb.png
    :alt:

  :doc:`/demos/04_data/d22_ListGPUs`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">List and choose GPUs</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Differentiable TIGRE operators inside a learned reconstruction network.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_d25_Pytorch_thumb.png
    :alt:

  :doc:`/demos/04_data/d25_Pytorch`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Use TIGRE from PyTorch</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Material decomposition from two scans at different tube voltages.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_d23_DualEnergy_thumb.png
    :alt:

  :doc:`/demos/04_data/d23_DualEnergy`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Dual-energy CBCT</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Reconstruct the SophiaBeads dataset, a tube of glass beads.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_SophiaBeads_demo_thumb.png
    :alt:

  :doc:`/demos/04_data/SophiaBeads_demo`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Real micro-CT data: SophiaBeads</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Compute an optimised proton radiography. See /user-guide/proton-ct.">

.. only:: html

  .. image:: /demos/04_data/images/thumb/sphx_glr_pRad_demo_thumb.png
    :alt:

  :doc:`/demos/04_data/pRad_demo`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Proton radiography</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>


.. toctree::
   :hidden:
   :includehidden:


   /demos/01_basics/index.rst
   /demos/02_algorithms/index.rst
   /demos/03_geometry/index.rst
   /demos/04_data/index.rst



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
