TIGRE further tuning
======

You may want to accelerate TIGRE more, and are wondering how to do so. 
This File gives advice on what to change in the source files of TIGRE in order to tune it to your specific computer.

Instructions are given for MATLAB, but are easily adaptable to Python.

We highly suggest you start with the first suggestion. It means 0.3s per operation, considering that a backprojection of 512^3 with projections of size 1000^2 x 500 take approximately 1s, 0.3s from that is a significant acceleration, it can be over 50% in some cases. 

The other suggestions may not improve your performance, or improve it very slightly. 

Current values in the source code are for GTX 1080 Ti.

## All CUDA calls, ~0.3s per GPU saved

All .cu files have a function called `checkFreeMemory` defined in the bottom of the file. 
This function checks each available GPU and queries how much free memory there is available.
The function is fundamental for the correct splitting without errors of the problem. 

**However**, if you know how much memory (in bytes) your GPUs have, you can replace this function by just a constant number, in bytes, of the amount of free memory your GPU has. Consider giving it only 95% of the actual free memory, and consider that other processes (such as your desktop, internet browser) may be using the GPU too. Adding a number greater than the available memory will cause the code to error. 

Just comment all the lines inside the function and add
`*mem_GPU_global= your number in bytes;` 

e.g. for a GTX 1080 Ti, with 11 GB of DRAM, this value is around 11721500000. We recommend using less than 95% of this value, something like 11000000000. Lower if other processes are also using your GPUs.

Compile again, and you save up to 0.3s per GPU, per CUDA call. This can be significant (512^3 image, 512^2 projections 512 angles take around this time to backproject). 

## TV-minimization (SART-TV, im3Ddenoise)

**Expected acceleration is not more than 5%**, likely 0%. So change it if time is really critical to you. 


 1. Define your expected geometry. 
 
 2. Create a sample image using `head=headPhantom(geo.nVoxel)+rand(geo.nVoxel.')*0.2` 
 
 3. Time it using `tic; tv=im3DDenoise(head,'TV',150,5); toc` 
 
 4. Open `Source/tvDenoise.cu` , lines 54,55 should define `MAX_BUFFER` and `BLOCK_SIZE`. 
 
 5. Modify `MAX_BUFFER` by steps of size 10, modify `BLOCK_SIZE` by steps of size 2. `MAX_BUFFER` is expected to have a bigger influence.
 
 6. Recompile and time. Find the appropriate `MAX_BUFFER` and `BLOCK_SIZE` for your GPU, System and image size
 
## All algorithms (Atb)

**Expected acceleration is not more than 10%**, likely 0%. So change it if time is really critical to you. The steps show how to improve FDK backprojection, similar steps are needed for matched backprojector.


 1. Define your expected geometry. 
 
 2. Create a sample projection set, geometry and angles.` 
 
 3. Time it using `tic; bp=Atb(projection,geo,angles); toc` 
 
 4. Open `Source/voxel_backprojection.cu` , lines 107,108 should define `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD`. 
 
 5. Modify `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD` by steps of size 1. `PROJ_PER_KERNEL*PROJ_PER_KERNEL*VOXEL_PER_THREAD` must not exceed your maximum number of threads, likely 1024
 
 6. Recompile and time again. Find the appropriate `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD` for your GPU, System and image size
 
## All algoritms (Ax)

**Expected acceleration is not more than 10%**, likely 0%. So change it if time is really critical to you. The steps show how to improve FDK backprojection, similar steps are needed for matched backprojector.

 1. Define your expected geometry. 
 
 2. Create a sample image, geometry and angles.` 
 
 3. Time it using `tic; fp=Ax(image,geo,angles); toc` 
 
 4. Open `Source/siddon_projection.cu` , lines 107,108 should define `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD`. 
 
 5. Modify `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD` by steps of size 1. `PROJ_PER_KERNEL*PROJ_PER_KERNEL*VOXEL_PER_THREAD` must not exceed your maximum number of threads, likely 1024
 
 6. Recompile and time again. Find the appropriate `PROJ_PER_KERNEL` and `VOXEL_PER_THREAD` for your GPU, System and image size
 
