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
 *           Sligtly modified by Ander Biguri
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
    void tvdenoising(const float* src, float* dst, float lambda,
            const float* spacing, const long* image_size, int maxIter){
        
        // Prepare for MultiGPU
        int deviceCount = 0;
        cudaGetDeviceCount(&deviceCount);
        cudaCheckErrors("Device query fail");
        if (deviceCount == 0) {
            mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPUselect","There are no available device(s) that support CUDA\n");
        }
        //
        // CODE assumes
        // 1.-All available devices are usable by this code
        // 2.-All available devices are equal, they are the same machine (warning trhown)
        int dev;
        char * devicenames;
        cudaDeviceProp deviceProp;
        
        for (dev = 0; dev < deviceCount; dev++) {
            cudaSetDevice(dev);
            cudaGetDeviceProperties(&deviceProp, dev);
            if (dev>0){
                if (strcmp(devicenames,deviceProp.name)!=0){
                    mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed. \n POCS_TV.cu line 277.");
                    break;
                }
            }
            devicenames=deviceProp.name;
        }
        
        // We don't know if the devices are being used. lets check that. and only use the amount of memory we need.
        size_t memfree;
        size_t memtotal;
        size_t mem_GPU_global;
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(dev);
            cudaMemGetInfo(&memfree,&memtotal);
            if(dev==0) mem_GPU_global=memfree;
            if(memfree<memtotal/2){
                mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
            }
            cudaCheckErrors("Check mem error");
            
            mem_GPU_global=(memfree<mem_GPU_global)?memfree:mem_GPU_global;
        }
        mem_GPU_global=(size_t)((double)mem_GPU_global*0.95);
        
        
        // %5 of free memory shoudl be enough, we have almsot no variables in these kernels
        size_t total_pixels              = image_size[0] * image_size[1]  * image_size[2] ;
        size_t mem_slice_image           = sizeof(float)* image_size[0] * image_size[1]  ;
        size_t mem_size_image            = sizeof(float)* total_pixels;
        
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
                // one more splot shoudl do the job, as its an edge case.
                splits++;
                //recompute for later
                slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits); // amountf of slices that fit on a GPU. Later we add 2 to these, as we need them for overlap
                mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            }
            
            // How many EXTRA buffer slices shoudl be able to fit in here??!?!
            mem_free=mem_GPU_global-(5*mem_img_each_GPU);
            unsigned int extra_buff=(mem_free/mem_slice_image);
            buffer_length=(extra_buff/2)/5; // we need double whatever this results in, rounded down.
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // Assert
            if (mem_GPU_global< 5*mem_img_each_GPU){
                mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","Bad assert. Logic behind spliting flawed! Please tell: ander.biguri@gmail.com\n");
            }
        }
        
        
        float* buffer;
        float* h_px, *h_py, *h_pz, *h_u;
        if(splits>1){
            mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:Memory","TV dneoising requires 5 times the image memory. Your GPU(s) do not have the required memory.\n This memory will be attempted to allocate on the CPU, Whic may fail or slow the computation by a very significant amount.\n If you want to kill the execution: CTRL+C");
                        
            cudaMallocHost((void**)&h_px,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
            cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
            
            cudaMallocHost((void**)&h_py,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
            cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
            
            cudaMallocHost((void**)&h_pz,image_size[0]*image_size[1]*image_size[2]*sizeof(float));
            cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
            h_u=dst;
//             if (h_u==NULL){
//                 mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising","Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");\
//             }

            
        }else{
            cudaMallocHost((void**)&buffer,image_size[0]*image_size[1]*sizeof(float));

        }
        // We shoudl be good to go memory wise.
        
        
        float** d_src   =(float**)malloc(deviceCount*sizeof(float*));
        float** d_u     =(float**)malloc(deviceCount*sizeof(float*));
        float** d_px    =(float**)malloc(deviceCount*sizeof(float*));
        float** d_py    =(float**)malloc(deviceCount*sizeof(float*));
        float** d_pz    =(float**)malloc(deviceCount*sizeof(float*));
        
        //Malloc
        for(dev=0;dev<deviceCount;dev++){
            cudaSetDevice(dev);
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
        

        
        // Allocate CPU buffer if needed, warn user if not.
        
        
        
        unsigned int curr_slices;
        unsigned long long curr_pixels;
        unsigned long long buffer_pixels=buffer_length*image_size[0]*image_size[1];
        
        
        float tau2, tau1;
        
        for(unsigned int i=0;i<maxIter;i+=(buffer_length)){
            
            for(unsigned int sp=0;sp<splits;sp++){
                
                // For each iteration we need to comptue all the image. The ordering of these loops
                // need to be like this due to the boudnign layers between slpits. If more than 1 split is needed
                // for each GPU then there is no other way that taking the entire memory out of GPU and putting it back.
                // If the memory can be shared ebtween GPUs fully without extra splits, then there is an easy way of syncronizing the memory
                
                // Copy image to memory
                size_t linear_idx_start;
                if(i==0){
                    for (dev = 0; dev < deviceCount; dev++){
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        curr_pixels=curr_slices*image_size[0]*image_size[1];
                        linear_idx_start=image_size[0]*image_size[1]*slices_per_split*(sp*deviceCount+dev);
                        
                        cudaSetDevice(dev);
                        cudaMemcpyAsync(d_src[dev]+buffer_pixels, &src[linear_idx_start]  , curr_pixels*sizeof(float), cudaMemcpyHostToDevice  );
                        cudaMemcpyAsync(d_u[dev]  +buffer_pixels, d_src[dev]+buffer_pixels, curr_pixels*sizeof(float), cudaMemcpyDeviceToDevice);
                        
                        
                        cudaMemset(d_px[dev], 0, mem_img_each_GPU);
                        cudaMemset(d_py[dev], 0, mem_img_each_GPU);
                        cudaMemset(d_pz[dev], 0, mem_img_each_GPU);

                        // if its not the last, copy also the intersection buffer.
                        if((sp*deviceCount+dev)<deviceCount*splits-1){
                            cudaMemcpyAsync(d_src[dev]+curr_pixels+buffer_pixels, src+linear_idx_start+curr_pixels  , buffer_pixels*sizeof(float), cudaMemcpyHostToDevice  );       
                            cudaMemcpyAsync(d_u[dev]  +curr_pixels+buffer_pixels, d_src[dev]+curr_pixels+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToDevice);

                        }
                        // if its not the first, copy also the intersection buffer.
                        if((sp*deviceCount+dev)){
                            cudaMemcpyAsync(d_src[dev], &src[linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice  );
                            cudaMemcpyAsync(d_u[dev]  , d_src[dev]                          , buffer_pixels*sizeof(float), cudaMemcpyDeviceToDevice);

                        }
                        
                    }
                }
                cudaDeviceSynchronize();
                cudaCheckErrors("Memcpy failure");
                // if we need to split and its not the first iteration, then we need to copy from Host memory.
                // d_src is the original image, with no change.
                if (splits>1 & i>0){
                    for (dev = 0; dev < deviceCount; dev++){
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        linear_idx_start=image_size[0]*image_size[1]*slices_per_split*(sp*deviceCount+dev);
                        curr_pixels=curr_slices*image_size[0]*image_size[1];
                        
                        cudaSetDevice(dev);
                        cudaMemcpyAsync(d_src[dev]+buffer_pixels, &src [linear_idx_start], curr_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        cudaMemcpyAsync(d_u [dev] +buffer_pixels, &h_u [linear_idx_start], curr_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        cudaMemcpyAsync(d_px[dev] +buffer_pixels, &h_px[linear_idx_start], curr_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        cudaMemcpyAsync(d_py[dev] +buffer_pixels, &h_py[linear_idx_start], curr_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        cudaMemcpyAsync(d_pz[dev] +buffer_pixels, &h_pz[linear_idx_start], curr_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        
                        // if its not the last, copy also the intersection buffer.
                        if((sp*deviceCount+dev)<deviceCount*splits-1){
                            
                            cudaMemcpyAsync(d_src[dev]+curr_pixels+buffer_pixels, &src [linear_idx_start]+curr_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_u [dev] +curr_pixels+buffer_pixels, &h_u [linear_idx_start]+curr_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_px[dev] +curr_pixels+buffer_pixels, &h_px[linear_idx_start]+curr_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_py[dev] +curr_pixels+buffer_pixels, &h_py[linear_idx_start]+curr_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_pz[dev] +curr_pixels+buffer_pixels, &h_pz[linear_idx_start]+curr_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            
                        }
                        // if its not the first, copy also the intersection buffer.
                        if((sp*deviceCount+dev)){
                            cudaMemcpyAsync(d_src[dev], &src [linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_u [dev] , &h_u [linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_px[dev] , &h_px[linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_py[dev] , &h_py[linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            cudaMemcpyAsync(d_pz[dev] , &h_pz[linear_idx_start]-buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                        }
                        
                    }
                    
                }
                cudaDeviceSynchronize();
                cudaCheckErrors("Memcpy failure on multi split");
                
                for(unsigned int ib=0;  (ib<(buffer_length)) && ((i+ib)<maxIter);  ib++){
                    
                    
                    
                    
                    tau2 = 0.3f + 0.02f * (i+ib);
                    tau1 = (1.f/tau2) * ((1.f/6.f) - (5.f/(15.f+(i+ib))));
                    // bdim and gdim
                    
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(dev);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        dim3 block(10, 10, 10);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (curr_slices+buffer_length*2+block.z-1)/block.z);
                        
                        update_u<<<grid, block>>>(d_src[dev], d_pz[dev], d_py[dev], d_px[dev], d_u[dev], tau1, lambda,
                                (long)(curr_slices+buffer_length*2), image_size[1],image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(dev);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        dim3 block(10, 10, 10);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (curr_slices+buffer_length*2+block.z-1)/block.z);
                        
                        update_p<<<grid, block>>>(d_u[dev], d_pz[dev], d_py[dev], d_px[dev], tau2,
                                (long)(curr_slices+buffer_length*2), image_size[1], image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                    
                }// END internal iter
                
                // Syncronize mathematics, make sure bounding pixels are correct
                cudaDeviceSynchronize();
                
                if(splits==1){
                    for(dev=0; dev<deviceCount;dev++){
                        if (dev<deviceCount-1){
                            // U
                            cudaSetDevice(dev+1);
                            cudaMemcpy(buffer, d_u[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_u[dev]+slices_per_split+buffer_pixels, buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            // px
                            cudaSetDevice(dev+1);
                            cudaMemcpy(buffer, d_px[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_px[dev]+slices_per_split+buffer_pixels, buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            // py
                            cudaSetDevice(dev+1);
                            cudaMemcpy(buffer, d_py[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_py[dev]+slices_per_split+buffer_pixels, buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            // pz
                            cudaSetDevice(dev+1);
                            cudaMemcpy(buffer, d_pz[dev+1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_pz[dev]+slices_per_split+buffer_pixels, buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            
                            
                        }
                        if (dev>0){
                            // U
                            cudaSetDevice(dev-1);
                            cudaMemcpy(buffer, d_u[dev-1]+slices_per_split+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_u[dev],buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            //px
                            cudaSetDevice(dev-1);
                            cudaMemcpy(buffer, d_py[dev-1]+slices_per_split+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_px[dev],buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            //py
                            cudaSetDevice(dev-1);
                            cudaMemcpy(buffer, d_py[dev-1]+slices_per_split+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_py[dev],buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            //pz
                            cudaSetDevice(dev-1);
                            cudaMemcpy(buffer, d_pz[dev-1]+slices_per_split+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                            cudaSetDevice(dev);
                            cudaMemcpy(d_pz[dev],buffer, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice);
                            
                            
                        }
                    }
                }else{
                    // We need to take it out :(
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(dev);
                        
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        linear_idx_start=image_size[0]*image_size[1]*slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*image_size[0]*image_size[1];
                        cudaMemcpy(&h_u[linear_idx_start],  d_u [dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                        cudaMemcpy(&h_px[linear_idx_start], d_px[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                        cudaMemcpy(&h_py[linear_idx_start], d_py[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                        cudaMemcpy(&h_pz[linear_idx_start], d_pz[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                        
                    }
                }
                cudaDeviceSynchronize();
                cudaCheckErrors("Memory gather error");
                
            }//END splits
        }//END main iter
        
        
        cudaDeviceSynchronize();
        cudaCheckErrors("TV minimization");
        
        if(splits==1){
            for(dev=0; dev<deviceCount;dev++){
                cudaSetDevice(dev);
                curr_slices=((dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*dev;
                total_pixels=curr_slices*image_size[0]*image_size[1];
                cudaMemcpy(dst+slices_per_split*image_size[0]*image_size[1]*dev, d_u[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
            }
        }
        cudaDeviceSynchronize();
        cudaCheckErrors("Copy result back");
        
        for(dev=0; dev<deviceCount;dev++){
            
            cudaFree(d_src[dev]);
            cudaFree(d_u [dev]);
            cudaFree(d_pz[dev]);
            cudaFree(d_py[dev]);
            cudaFree(d_px[dev]);
        }
        if(splits>1){
           cudaFreeHost(h_px); 
           cudaFreeHost(h_py); 
           cudaFreeHost(h_pz); 
        }else{
            cudaFreeHost(buffer);
        }
            
        cudaDeviceReset();
    }