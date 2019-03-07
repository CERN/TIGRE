TIGRE further tunning
======

You may want to accelerate TIGRE more, and are wondering how to do so. 
This File gives advice on what to change in the source files of TIGRE in order to tune it to your specific computer.

Instructions are give for MATLAB, but are easily adaptable to Python.

We highly suggest changing the following optimization. It means 0.3s per operation, considering that 
a backprojection of 512^3 with projections of size 1000^2 x 500 take a pproximately 1s, 0.3s from that is a significant
acceleration, it can be over 50% in some cases. 

## All CUDA calls, ~0.3s per GPU saved

All .cu files have a function called `checkFreeMemory` defined in the bottom of the file. 
This function checks each available GPU and queries how much free memory there is available.
The fucntion is fundamental for the correct splitting without errors of the problem. 

*However*, if you know how much memory (in bytes) your GPUs have, you can replace this function by
just a constant number, in bytes, of the amount of free memory your GPU has. Consider giving it only 95% of the actuall free memory,
and consider that other processes (such as your desktop, internet browser) may be using the GPU too. Adding a number greated than the available memory 
will cause the code to error. 

Just comment all the lines inside the function and add

`*mem_GPU_global= your number in bytes;`

Compile again, and you saved 0.3s per GPU, per cuda call. This can be significant. 

## TV-minimization (SART-TV, im3Ddenoise)

*Expected acceleration is not more than 5%*, so change it if time is really critical to you. 


 1. Define your expected geometry. 
 
 2. Create a sample image using `head=headPhantom(geo.nVoxel)+rand(geo.nVoxel.')*0.2` 
 
 3. time it using `tic; tv=im3DDenoise(head,'TV',150,5); toc` 
 
 4. Open `Source/tvDenoise.cu` , lines 54,55 should define `MAX_BUFFER` and `BLOCK_SIZE`. 
 
 5. Modify `MAX_BUFFER` by steps of size 10, modify `BLOCK_SIZE` by steps of size 2. `MAX_BUFFER` is expected to have a bigger influence.
 
 6. Recompile and time. Find the appropiate `MAX_BUFFER` and `BLOCK_SIZE` for your GPU, System and image size
 
 






#
