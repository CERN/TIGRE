"""
Load the sample data
====================

Load the head phantom and the other images that ship with TIGRE.

Source: `Python/demos/d02_SampleData.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d02_SampleData.py>`__ · `MATLAB/Demos/d02_SampleData.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d02_SampleData.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d02_SampleData.m
      :language: matlab
"""

import tigre
from tigre.utilities import sample_loader
from tigre.utilities import sl3d

# %%
# Define geometry
# ---------------

geo = tigre.geometry_default(high_resolution=False)

# %%

# --------------------------------------------------------------------------
# 3D Shepp-Logan
#
# Type of shepp logan phantom. The shape will be the same, but Hounsfield
# values will be different
# Options:
#  phantom_type: str
#      One of {"kak-slaney", "yu-ye-wang", "toft-schabel"}, optional
#      Default is "yu-ye-wang"
#      The type of phantom.
#  size_out : int or list whose length is three, optional
#      Default is [128, 128, 128]
#      The number of voxels of phantom. [nVoxelZ, nVoxelY, nVoxelX]
# --------------------------------------------------------------------------
# phantom_type="kak-slaney"    # Air is 0. Water is 1. Proportional to Hounsfield value.
phantom_type = (
    "yu-ye-wang"  # Default of Python TIGRE Shepp-Logan phantom. Improved visual perception
)
# phantom_type="toft-schabel"  # Default of MATLAB TIGRE Shepp-Logan phantom.
shepp = sl3d.shepp_logan_3d(
    geo.nVoxel, phantom_type=phantom_type
)  # Default are 128^3 and "yu-ye-wang"
tigre.plotImg(shepp, dim="z")

# %%
# load head phantom
# -----------------

head = sample_loader.load_head_phantom(geo.nVoxel)
# check the shape
print(head.shape)
# show it
tigre.plotImg(head, dim="z")
