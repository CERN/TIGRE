/*-------------------------------------------------------------------------
 *
 * MATLAB MEX  functions for TV image denoising. Check inputs and parses
 * MATLAB data to C++ data.
 *
 *
 * CODE by   Imanol Luengo
 *           PhD student University of Nottingham
 *           imaluengo@gmail.com
 *           2015
 *           Modified by Ander Biguri for multi-GPU 
 * ---------------------------------------------------------------------------
 * ---------------------------------------------------------------------------
 * Copyright (c) 2015, University of Bath and CERN- European Organization for
 * Nuclear Research
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * ---------------------------------------------------------------------------
 *
 * Contact: tigre.toolbox@gmail.com
 * Codes  : https://github.com/CERN/TIGRE
 * ---------------------------------------------------------------------------
 */



// http://gpu4vision.icg.tugraz.at/papers/2010/knoll.pdf#pub47
#define MAXTREADS 1024
#define MAX_BUFFER 60
#define BLOCK_SIZE 10  // BLOCK_SIZE^3 must be smaller than MAXTREADS

#include "tvdenoising.hpp"
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                cudaDeviceReset();\
                mexPrintf("%s \n",msg);\
                        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
    
    __device__ __inline__
            float divergence(const float* pz, const float* py, const float* px,
            long z, long y, long x, long depth, long rows, long cols,
            float dz, float dy, float dx)
    {
        long size2d = rows*cols;
        long idx = z * size2d + y * cols + x;
        float _div = 0.0f;
        
        if ( z - 1 >= 0 ) {
            _div += (pz[idx] - pz[(z-1)*size2d + y*cols + x]) / dz;
        } else {
            _div += pz[idx];
        }
        
        if ( y - 1 >= 0 ) {
            _div += (py[idx] - py[z*size2d + (y-1)*cols + x]) / dy;
        } else {
            _div += py[idx];
        }
        
        if ( x - 1 >= 0 ) {
            _div += (px[idx] - px[z*size2d + y*cols + (x-1)]) / dx;
        } else {
            _div += px[idx];
        }
        
        return _div;
    }
    
    __device__ __inline__
            void gradient(const float* u, float* grad,
            long z, long y, long x,
            long depth, long rows, long cols,
            float dz, float dy, float dx)
    {
        long size2d = rows*cols;
        long idx = z * size2d + y * cols + x;
        
        float uidx = u[idx];
        
        if ( z + 1 < depth ) {
            grad[0] = (u[(z+1)*size2d + y*cols + x] - uidx) / dz;
        }
        
        if ( y + 1 < rows ) {
            grad[1] = (u[z*size2d + (y+1)*cols + x] - uidx) / dy;
        }
        
        if ( x + 1 < cols ) {
            grad[2] = (u[z*size2d + y*cols + (x+1)] - uidx) / dx;
        }
    }
    
    
    __global__
            void update_u(const float* f, const float* pz, const float* py, const float* px, float* u,
            float tau, float lambda,
            long depth, long rows, long cols,
            float dz, float dy, float dx)
    {
        long x = threadIdx.x + blockIdx.x * blockDim.x;
        long y = threadIdx.y + blockIdx.y * blockDim.y;
        long z = threadIdx.z + blockIdx.z * blockDim.z;
        long idx = z * rows * cols + y * cols + x;
        
        if ( x >= cols || y >= rows || z >= depth )
            return;
        
        float _div = divergence(pz, py, px, z, y, x, depth, rows, cols, dz, dy, dx);
        
        u[idx] = u[idx] * (1.0f - tau) + tau * (f[idx] + (1.0f/lambda) * _div);
    }
    
    
    __global__
            void update_p(const float* u, float* pz, float* py, float* px,
            float tau, long depth, long rows, long cols,
            float dz, float dy, float dx)
    {
        long x = threadIdx.x + blockIdx.x * blockDim.x;
        long y = threadIdx.y + blockIdx.y * blockDim.y;
        long z = threadIdx.z + blockIdx.z * blockDim.z;
        long idx = z * rows * cols + y * cols + x;
        
        if ( x >= cols || y >= rows || z >= depth )
            return;
        
        float grad[3] = {0,0,0}, q[3];
        gradient(u, grad, z, y, x, depth, rows, cols, dz, dy, dx);
        
        q[0] = pz[idx] + tau * grad[0];
        q[1] = py[idx] + tau * grad[1];
        q[2] = px[idx] + tau * grad[2];
        
        float norm = fmaxf(1.0f, sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]));
        
        pz[idx] = q[0] / norm;
        py[idx] = q[1] / norm;
        px[idx] = q[2] / norm;
    }
    
    
// Main function
    void tvdenoising(float* src, float* dst, float lambda,
            const float* spacing, const long* image_size, int maxIter, const GpuIds& gpuids) {
        
        // Prepare for MultiGPU
        int deviceCount = gpuids.GetLength();
        cudaCheckErrors("Device query fail");
        if (deviceCount == 0) {
            mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPUselect","There are no available device(s) that support CUDA\n");
        }
        //
        // CODE assumes
        // 1.-All available devices are usable by this code
        // 2.-All available devices are equal, they are the same machine (warning thrown)
        // Check the available devices, and if they are the same
        if (!gpuids.AreEqualDevices()) {
            mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
        }
        int dev;

        // We don't know if the devices are being used. lets check that. and only use the amount of memory we need.
        
        size_t mem_GPU_global;
        checkFreeMemory(gpuids, &mem_GPU_global);
        
        
        // %5 of free memory should be enough, we have almost no variables in these kernels
       size_t total_pixels           = image_size[0] * image_size[1] * image_size[2] ;
       const size_t pixels_per_slice = image_size[0] * image_size[1] ;
       const size_t mem_slice_image  = sizeof(float)* pixels_per_slice  ;
       const size_t mem_size_image   = sizeof(float)* total_pixels;
        
        // Decide how are we handling the distribution of computation
        size_t mem_img_each_GPU;
        
        unsigned int buffer_length=1;
        //Does everything fit in the GPU?
        unsigned int slices_per_split;
        unsigned int splits=1; // if the number does not fit in an uint, you have more serious trouble than this.
        if(mem_GPU_global> 5*mem_size_image+5*mem_slice_image*buffer_length*2){
            // We only need to split if we have extra GPUs
            slices_per_split=(image_size[2]+deviceCount-1)/deviceCount;
            mem_img_each_GPU=mem_slice_image*(  (image_size[2]+deviceCount-1)/deviceCount  + buffer_length*2);
        }else{
            // As mem_auxiliary is not expected to be a large value (for a 2000^3 image is around 28Mbytes), lets for now assume we need it all
            size_t mem_free=mem_GPU_global;
            
            splits=(unsigned int)(ceil(((float)(5*mem_size_image)/(float)(deviceCount))/mem_free));
            // Now, there is an overhead here, as each splits should have 2 slices more, to accoutn for overlap of images.
            // lets make sure these 2 slices fit, if they do not, add 1 to splits.
            slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits);
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // if the new stuff does not fit in the GPU, it measn we are in the edge case where adding that extra slice will overflow memory
            if (mem_GPU_global< 5*mem_img_each_GPU){
                // one more split should do the job, as its an edge case.
                splits++;
                //recompute for later
                slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits); // amount of slices that fit on a GPU. Later we add 2 to these, as we need them for overlap
                mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            }
            
            // How many EXTRA buffer slices should be able to fit in here??!?!
            mem_free=mem_GPU_global-(5*mem_img_each_GPU);
            unsigned int extra_buff=(mem_free/mem_slice_image);
            buffer_length=(extra_buff/2)/5; // we need double whatever this results in, rounded down.
            
            buffer_length=min(MAX_BUFFER,buffer_length);
            
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // Assert
            if (mem_GPU_global< 5*mem_img_each_GPU){
                mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","Bad assert. Logic behind spliting flawed! Please tell: ander.biguri@gmail.com\n");
            }
        }
        
        
        // Lets try to make the host memory pinned:
        // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
        int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
        cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
        if (isHostRegisterSupported & splits>1){
            cudaHostRegister(src ,image_size[2]*image_size[1]*image_size[0]*sizeof(float),cudaHostRegisterPortable);
            cudaHostRegister(dst ,image_size[2]*image_size[1]*image_size[0]*sizeof(float),cudaHostRegisterPortable);
        }
        cudaCheckErrors("Error pinning memory");
        
        
        
        // Lets allocate auxiliary  variables.
        float* buffer_u, *buffer_px, *buffer_py, *buffer_pz;
        float* h_px, *h_py, *h_pz, *h_u;
        if(splits>1){
            
            //These take A LOT of memory and A LOT of time to use. If we can avoid using them, better.
            if (buffer_length<maxIter){ // if we do only 1 big iter, they are not needed.
                mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:Memory","TV dneoising requires 5 times the image memory. Your GPU(s) do not have the required memory.\n This memory will be attempted to allocate on the CPU, Whic may fail or slow the computation by a very significant amount.\n If you want to kill the execution: CTRL+C");
                
                cudaMallocHost((void**)&h_px,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
                
                cudaMallocHost((void**)&h_py,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
                
                cudaMallocHost((void**)&h_pz,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
            }
            h_u=dst;
        }else{
            cudaMallocHost((void**)&buffer_u,  pixels_per_slice*sizeof(float));
            cudaMallocHost((void**)&buffer_px, pixels_per_slice*sizeof(float));
            cudaMallocHost((void**)&buffer_py, pixels_per_slice*sizeof(float));
            cudaMallocHost((void**)&buffer_pz, pixels_per_slice*sizeof(float));
            
        }
        // We should be good to go memory wise.
        
        
        float** d_src   =(float**)malloc(deviceCount*sizeof(float*));
        float** d_u     =(float**)malloc(deviceCount*sizeof(float*));
        float** d_px    =(float**)malloc(deviceCount*sizeof(float*));
        float** d_py    =(float**)malloc(deviceCount*sizeof(float*));
        float** d_pz    =(float**)malloc(deviceCount*sizeof(float*));
        
        //Malloc
        for(dev=0;dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            // F
            cudaMalloc((void**)&d_src[dev], mem_img_each_GPU);
            // U
            cudaMalloc((void**)&d_u [dev],  mem_img_each_GPU);
            // PX
            cudaMalloc((void**)&d_px[dev],  mem_img_each_GPU);
            // PY
            cudaMalloc((void**)&d_py[dev],  mem_img_each_GPU);
            // PZ
            cudaMalloc((void**)&d_pz[dev],  mem_img_each_GPU);
        }
        cudaDeviceSynchronize();
        cudaCheckErrors("Malloc  error");
        
        
        // Create streams
        int nStream_device=5;
        int nStreams=deviceCount*nStream_device;
        cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
        
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            for (int i = 0; i < nStream_device; ++i){
                cudaStreamCreate(&stream[i+dev*nStream_device]);
            }
        }
        cudaCheckErrors("Stream creation fail");
        
        
        
        
        // Allocate CPU buffer if needed, warn user if not.
        
        
        
        unsigned int curr_slices;
        unsigned long long curr_pixels;
        size_t linear_idx_start;
        
        unsigned long long buffer_pixels=buffer_length*pixels_per_slice;
        
        unsigned long long* offset_device=(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* offset_host  =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* bytes_device =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        bool is_first_chunk;
        bool is_last_chunk;
        
        float tau2, tau1;

        for(unsigned int i=0;i<maxIter;i+=(buffer_length)){
           
            for(unsigned int sp=0;sp<splits;sp++){
                
                // For each iteration we need to compute all the image. The ordering of these loops
                // need to be like this due to the bounding layers between splits. If more than 1 split is needed
                // for each GPU then there is no other way that taking the entire memory out of GPU and putting it back.
                // If the memory can be shared between GPUs fully without extra splits, then there is an easy way of synchronizing the memory
                
                // Copy image to memory
                for (dev = 0; dev < deviceCount; dev++){
                    // Precompute indices and needed bytes
                    curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                    curr_pixels=curr_slices*pixels_per_slice;
                    linear_idx_start=pixels_per_slice*slices_per_split*(sp*deviceCount+dev);
                    
                    // Check if its the first or last chunck
                    is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                    is_first_chunk=!(sp*deviceCount+dev);
                    
                    // lets compute where we start copyes and how much. This avoids 3 calls to Memcpy
                    offset_device[dev]=buffer_pixels*is_first_chunk;
                    offset_host[dev]=linear_idx_start-buffer_pixels*!is_first_chunk;
                    bytes_device[dev]=curr_pixels+buffer_pixels*!is_first_chunk+buffer_pixels*!is_last_chunk;
                }
                if(i==0){
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaMemcpyAsync(d_src[dev]+offset_device[dev], src+offset_host[dev]  , bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        // All these are async
                        cudaMemcpyAsync(d_u[dev]  +offset_device[dev], d_src[dev]+offset_device[dev], bytes_device[dev]*sizeof(float), cudaMemcpyDeviceToDevice,stream[dev*nStream_device+1]);
                        cudaMemsetAsync(d_px[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                        cudaMemsetAsync(d_py[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                        cudaMemsetAsync(d_pz[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                    }
                    // we need all the stream to finish
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                    }
                    cudaCheckErrors("Memcpy failure");
                    
                }
                // if we need to split and its not the first iteration, then we need to copy from Host memory.
                // d_src is the original image, with no change.
                if (splits>1 & i>0){

                    for (dev = 0; dev < deviceCount; dev++){   
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+1]);
                        cudaMemcpyAsync(d_u [dev] +offset_device[dev], h_u +offset_host[dev],  bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                       
                    }

                    for (dev = 0; dev < deviceCount; dev++){   
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+2]);
                        cudaMemcpyAsync(d_px[dev]+offset_device[dev], h_px+offset_host[dev],  bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+2]);
                       
                    }
                    for (dev = 0; dev < deviceCount; dev++){   
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+3]);
                        cudaMemcpyAsync(d_py[dev] +offset_device[dev], h_py+offset_host[dev],  bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+3]);
                        
                    }
                    for (dev = 0; dev < deviceCount; dev++){   
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+4]);
                        cudaMemcpyAsync(d_pz[dev] +offset_device[dev], h_pz+offset_host[dev],  bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+4]);
                        
                    } 
                    for (dev = 0; dev < deviceCount; dev++){   

                        
                        cudaStreamSynchronize(stream[dev*nStream_device+1]);
                        cudaMemcpyAsync(d_src[dev]+offset_device[dev], src +offset_host[dev],  bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                        

                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                        cudaCheckErrors("Memcpy failure on multi split");
                    }
                }
                
                // Inter interations.
                for(unsigned int ib=0;  (ib<(buffer_length)) && ((i+ib)<maxIter);  ib++){

                    tau2 = 0.3f + 0.02f * (i+ib);
                    tau1 = (1.f/tau2) * ((1.f/6.f) - (5.f/(15.f+(i+ib))));
                    // bdim and gdim
                    
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (curr_slices+buffer_length*2+block.z-1)/block.z);
                        
                        update_u<<<grid, block,0,stream[dev*nStream_device]>>>(d_src[dev], d_pz[dev], d_py[dev], d_px[dev], d_u[dev], tau1, lambda,
                                (long)(curr_slices+buffer_length*2), image_size[1],image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (curr_slices+buffer_length*2+block.z-1)/block.z);
                        
                        update_p<<<grid, block,0,stream[dev*nStream_device]>>>(d_u[dev], d_pz[dev], d_py[dev], d_px[dev], tau2,
                                (long)(curr_slices+buffer_length*2), image_size[1], image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                }// END internal iter
                
                // Synchronize mathematics, make sure bounding pixels are correct
                for(dev=0; dev<deviceCount;dev++){
                    cudaSetDevice(gpuids[dev]);
                    cudaDeviceSynchronize();
                }
                if(splits==1){
                    for(dev=0; dev<deviceCount;dev++){
                        if (dev<deviceCount-1){
                            // U
                            cudaSetDevice(gpuids[dev+1]);
                            cudaMemcpyAsync(buffer_u , d_u[dev+1] , buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device+1]);
                            cudaMemcpyAsync(buffer_px, d_px[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device+2]);
                            cudaMemcpyAsync(buffer_py, d_py[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device+3]);
                            cudaMemcpyAsync(buffer_pz, d_pz[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device+4]);

                            
                            cudaSetDevice(gpuids[dev]);
                            cudaStreamSynchronize(stream[(dev+1)*nStream_device+1]);
                            cudaMemcpyAsync(d_u[dev] +slices_per_split*pixels_per_slice+buffer_pixels, buffer_u , buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+1]);
                            cudaStreamSynchronize(stream[(dev+1)*nStream_device+2]);
                            cudaMemcpyAsync(d_px[dev]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_px, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+2]);
                            cudaStreamSynchronize(stream[(dev+1)*nStream_device+3]);
                            cudaMemcpyAsync(d_py[dev]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_py, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+3]);
                            cudaStreamSynchronize(stream[(dev+1)*nStream_device+4]);
                            cudaMemcpyAsync(d_pz[dev]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_pz, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+4]);
                            
                            
                        }
                        cudaDeviceSynchronize();
                        if (dev>0){
                            // U
                            cudaSetDevice(gpuids[dev-1]);
                            cudaMemcpyAsync(buffer_u,  d_u[dev-1] +slices_per_split*pixels_per_slice+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev-1)*nStream_device+1]);
                            cudaMemcpyAsync(buffer_px, d_px[dev-1]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev-1)*nStream_device+2]);
                            cudaMemcpyAsync(buffer_py, d_py[dev-1]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev-1)*nStream_device+3]);
                            cudaMemcpyAsync(buffer_pz, d_pz[dev-1]+slices_per_split*pixels_per_slice+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev-1)*nStream_device+4]);
                            
                            
                            cudaSetDevice(gpuids[dev]);
                            cudaStreamSynchronize(stream[(dev-1)*nStream_device+1]);
                            cudaMemcpyAsync(d_u[dev] ,buffer_u , buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+1]);
                            cudaStreamSynchronize(stream[(dev-1)*nStream_device+2]);
                            cudaMemcpyAsync(d_px[dev],buffer_px, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+2]);
                            cudaStreamSynchronize(stream[(dev-1)*nStream_device+3]);
                            cudaMemcpyAsync(d_py[dev],buffer_py, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+3]);
                            cudaStreamSynchronize(stream[(dev-1)*nStream_device+4]);
                            cudaMemcpyAsync(d_pz[dev],buffer_pz, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev)*nStream_device+4]);
                            
                            
                        }
                    }
                }else{
                    // We need to take it out :(
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices      = ((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        linear_idx_start = pixels_per_slice*slices_per_split*(sp*deviceCount+dev);
                        total_pixels     = curr_slices*pixels_per_slice;
                        cudaMemcpyAsync(&h_u[linear_idx_start],  d_u [dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                    }
                    if ((i+buffer_length)<maxIter){ // If its the last iteration, we don't need to get these out.
                        for(dev=0; dev<deviceCount;dev++){
                            cudaSetDevice(gpuids[dev]);
                            curr_slices      = ((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                            linear_idx_start = pixels_per_slice*slices_per_split*(sp*deviceCount+dev);
                            total_pixels     = curr_slices*pixels_per_slice;
                            cudaMemcpyAsync(&h_px[linear_idx_start], d_px[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+2]);
                            cudaMemcpyAsync(&h_py[linear_idx_start], d_py[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+3]);
                            cudaMemcpyAsync(&h_pz[linear_idx_start], d_pz[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+4]);
                            
                        }
                    }

                }
            }//END splits
        }//END main iter
        
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("TV minimization");
        
        if(splits==1){
            for(dev=0; dev<deviceCount;dev++){
                cudaSetDevice(gpuids[dev]);
                curr_slices  = ((dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*dev;
                total_pixels = curr_slices*pixels_per_slice;
                cudaMemcpyAsync(dst+slices_per_split*pixels_per_slice*dev, d_u[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
            }
        }
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("Copy result back");
        for(dev=0; dev<deviceCount;dev++){
            
            cudaFree(d_src[dev]);
            cudaFree(d_u [dev]);
            cudaFree(d_pz[dev]);
            cudaFree(d_py[dev]);
            cudaFree(d_px[dev]);
        }
        if(splits>1 && buffer_length<maxIter){
            cudaFreeHost(h_px);
            cudaFreeHost(h_py);
            cudaFreeHost(h_pz);
        }else if(splits==1){
            cudaFreeHost(buffer_u);
            cudaFreeHost(buffer_px);
            cudaFreeHost(buffer_py);
            cudaFreeHost(buffer_pz);
        }
        
        for (int i = 0; i < nStreams; ++i)
           cudaStreamDestroy(stream[i]) ;

        if (isHostRegisterSupported & splits>1){
            cudaHostUnregister(src);
            cudaHostUnregister(dst);
        }
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("Copy free ");
        
    }
    
    
void checkFreeMemory(const GpuIds& gpuids,size_t *mem_GPU_global){
        size_t memfree;
        size_t memtotal;
        const int deviceCount = gpuids.GetLength();

        for (int dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaMemGetInfo(&memfree,&memtotal);
            if(dev==0) *mem_GPU_global=memfree;
            if(memfree<memtotal/2){
                mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
            }
            cudaCheckErrors("Check mem error");
            
            *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
        }
        *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
        
        //*mem_GPU_global= insert your known number here, in bytes.
}
