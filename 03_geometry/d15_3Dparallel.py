"""
Parallel-beam geometry
======================

Reconstruct 3D parallel-beam data, as from a synchrotron.

Source: `Python/demos/d15_3Dparallel.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d15_3Dparallel.py>`__ · `MATLAB/Demos/d15_3Dparallel.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d15_3Dparallel.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d15_3Dparallel.m
      :language: matlab
"""

# %%
# Initialize
# ----------

import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs

# %%
# Geometry
# --------

# While you can defien you rown geometry, the easiest way is:

# just give as nVoxel the same number as your projection size.
geo = tigre.geometry(
    mode="parallel", nVoxel=np.array([512, 512, 512])
)  # Parallel beam geometry does not require anything other than the image size.

# %%
# Load data and generate projections
# ----------------------------------

# define angles
angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
# Load thorax phantom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

# recon
imgFBP = algs.fbp(projections, geo, angles)
imgOSSART = algs.ossart(projections, geo, angles, 40)

# plot
tigre.plotImg(np.concatenate([head, imgFBP, imgOSSART], axis=1), dim="z", step=3)
