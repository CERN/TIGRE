# CUDACBCT

Code to compute projection and backprojection of a Cone Beam Computed Tomography geometry from Matlab, using GPU acceleration.

Before use, compile with:

mex -v -largeArrayDims ./Mex_files/minTV.cpp ./Mex_files/POCS_TV.cu
mex -v -largeArrayDims ./Mex_files/tvDenoise.cpp ./Mex_files/tvdenoising.cu
mex -v -largeArrayDims ./Mex_files/Ax.cpp ./Mex_files/ray_interpolated_projection.cu ./Mex_files/Siddon_projection.cu
mex -v -largeArrayDims ./Mex_files/Atb.cpp ./Mex_files/voxel_backprojection.cu ./Mex_files/voxel_backprojection2.cu

** FEX submissions used:**

https://www.mathworks.com/matlabcentral/fileexchange/50974-3d-shepp-logan-phantom
