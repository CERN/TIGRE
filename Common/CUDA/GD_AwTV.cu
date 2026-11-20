/*-------------------------------------------------------------------------
 *
 * CUDA functions for Steepest descend in POCS-type algorithms.
 *
 * This file will iteratively minimize by steepest descend the total variation
 * of the input image, with the parameters given, using GPUs.
 *
 * CODE by       Ander Biguri
 *
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







#define MAXTHREADS 1024
#define MAX_BUFFER 60

#include "GD_AwTV.hpp"
#include "gpuUtils.hpp"
#include "tigre_release_on_error.hpp"




#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                TIGRE_MATLAB_ONLY(cudaDeviceReset());\
                tigreCudaError("CBCT:CUDA:GD_TV", msg, cudaGetErrorString(__err));\
        } \
} while (0)
    
// CUDA kernels
//https://stackoverflow.com/questions/21332040/simple-cuda-kernel-optimization/21340927#21340927
    __global__ void divideArrayScalar(float* vec,float scalar,const size_t n){
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]/=scalar;
        }
    }
    __global__ void multiplyArrayScalar(float* vec,float scalar,const size_t n)
    {
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]*=scalar;
        }
    }
    __global__ void substractArrays(float* vec,float* vec2,const size_t n)
    {
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]-=vec2[i];
        }
    }
    
    // Backward differences at voxel (z,y,x) of the buffer; 0 outside the image. In z the image is the slots
    // [zlo,zhi) of the buffer: its halos hold the neighbouring GPUs' or pieces' slices, but there is nothing
    // beyond the first and the last slice of the whole image, as in x and y (see GD_TV.cu).
    __device__ __inline__
            void gradient(const float* u, float* grad,
            long z, long y, long x,
            long zlo, long zhi, long rows, long cols){
        if (z < zlo || z >= zhi || y >= rows || x >= cols)
            return;
        unsigned long size2d = rows*cols;
        size_t idx = (size_t)z * size2d + y * cols + x;
        
        float uidx = u[idx];
        
        if ( z - 1 >= zlo) {
            grad[0] = (uidx-u[(size_t)(z-1)*size2d + y*cols + x]) ;
        }
        
        if ( y - 1 >= 0){
            grad[1] = (uidx-u[(size_t)z*size2d + (y-1)*cols + x]) ;
        }
        
        if ( x - 1 >= 0) {
            grad[2] = (uidx-u[(size_t)z*size2d + y*cols + (x-1)]);
        }
    }
    
    __global__ void gradientTV(const float* f, float* dftv,
            long depth, long zlo, long zhi, long rows, long cols,const float delta){
        unsigned long x = threadIdx.x + blockIdx.x * blockDim.x;
        unsigned long y = threadIdx.y + blockIdx.y * blockDim.y;
        unsigned long z = threadIdx.z + blockIdx.z * blockDim.z;
        size_t idx = (size_t)z * (rows * cols) + y * cols + x;
        if ( x >= cols || y >= rows || z >= depth )
            return;
        
        
        float df[3] ={0.f,0.f,0.f};
        float dfi[3]={0.f,0.f,0.f}; // dfi== \partial f_{i+1,j,k}
        float dfj[3]={0.f,0.f,0.f};
        float dfk[3]={0.f,0.f,0.f};
        gradient(f,df  ,z  ,y  ,x  , zlo,zhi,rows,cols);
        gradient(f,dfi ,z  ,y  ,x+1, zlo,zhi,rows,cols);
        gradient(f,dfj ,z  ,y+1,x  , zlo,zhi,rows,cols);
        gradient(f,dfk ,z+1,y  ,x  , zlo,zhi,rows,cols);
        float eps=0.00000001; //% avoid division by zero
        
        float wx=__expf(-(df[0]/delta)*(df[0]/delta));
        float wy=__expf(-(df[1]/delta)*(df[1]/delta));
        float wz=__expf(-(df[2]/delta)*(df[2]/delta));
        
        float wxi=__expf(-(dfi[0]/delta)*(dfi[0]/delta));
        float wyi=__expf(-(dfi[1]/delta)*(dfi[1]/delta));
        float wzi=__expf(-(dfi[2]/delta)*(dfi[2]/delta));
        
        float wxj=__expf(-(dfj[0]/delta)*(dfj[0]/delta));
        float wyj=__expf(-(dfj[1]/delta)*(dfj[1]/delta));
        float wzj=__expf(-(dfj[2]/delta)*(dfj[2]/delta));
        
        float wxk=__expf(-(dfk[0]/delta)*(dfk[0]/delta));
        float wyk=__expf(-(dfk[1]/delta)*(dfk[1]/delta));
        float wzk=__expf(-(dfk[2]/delta)*(dfk[2]/delta));

        
        // this hsould do the trick I think
        
        dftv[idx]=(wx*df[0]+wy*df[1]+wz*df[2])/(sqrt(wx*df[0] *df[0] +wy*df[1] *df[1] +wz*df[2] *df[2])+eps)
        -wzi*dfi[2]/(sqrt(wxi*dfi[0]*dfi[0]+wyi*dfi[1]*dfi[1]+wzi*dfi[2]*dfi[2]) +eps)     // I wish I coudl precompute this, but if I do then Id need to recompute the gradient.
        -wyj*dfj[1]/(sqrt(wxj*dfj[0]*dfj[0]+wyj*dfj[1]*dfj[1]+wzj*dfj[2]*dfj[2]) +eps)
        -wxk*dfk[0]/(sqrt(wxk*dfk[0]*dfk[0]+wyk*dfk[1]*dfk[1]+wzk*dfk[2]*dfk[2]) +eps);
        
    
        return;
        
    }
    
    template<typename T>
    __device__ void warpReduce(volatile T *sdata, size_t tid) {
        sdata[tid] += sdata[tid + 32];
        sdata[tid] += sdata[tid + 16];
        sdata[tid] += sdata[tid + 8];
        sdata[tid] += sdata[tid + 4];
        sdata[tid] += sdata[tid + 2];
        sdata[tid] += sdata[tid + 1];
    }
    
    __global__ void  reduceNorm2(float *g_idata, float *g_odata, size_t n){
        extern __shared__ volatile float sdata[];
        //http://stackoverflow.com/a/35133396/1485872
        size_t tid = threadIdx.x;
        size_t i = (size_t)blockIdx.x*blockDim.x + tid;
        size_t gridSize = blockDim.x*gridDim.x;
        float mySum = 0;
        float value=0;
        while (i < n) {
            value=g_idata[i]; //avoid reading twice
            mySum += value*value;
            i += gridSize;
        }
        sdata[tid] = mySum;
        __syncthreads();
        
        if (tid < 512)
            sdata[tid] += sdata[tid + 512];
        __syncthreads();
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
        
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
        
        if (tid <  64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
        
        
#if (CUDART_VERSION >= 9000)
        if ( tid < 32 )
        {
            mySum = sdata[tid] + sdata[tid + 32];
            for (int offset = warpSize/2; offset > 0; offset /= 2) {
                mySum += __shfl_down_sync(0xFFFFFFFF, mySum, offset,32);
            }
        }
#else
        if (tid < 32) {
            warpReduce(sdata, tid);
            mySum = sdata[0];
        }
#endif
        if (tid == 0) g_odata[blockIdx.x] = mySum;
    }
    
    // Adds the per-block partial sums of reduceNorm2 in double (#856): each GPU returns its part of the norm
    // in double and the host adds the parts in double, so the total barely depends on how the image is split.
    __global__ void  reduceSum(float *g_idata, double *g_odata, size_t n){
        extern __shared__ volatile double sdata_d[];
        volatile double* sdata = sdata_d;
        //http://stackoverflow.com/a/35133396/1485872
        size_t tid = threadIdx.x;
        size_t i = (size_t)blockIdx.x*blockDim.x + tid;
        size_t gridSize = blockDim.x*gridDim.x;
        double mySum = 0;
        // float value=0;
        while (i < n) {
            mySum += g_idata[i];
            i += gridSize;
        }
        sdata[tid] = mySum;
        __syncthreads();
        
        if (tid < 512)
            sdata[tid] += sdata[tid + 512];
        __syncthreads();
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
        
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
        
        if (tid <  64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
        
        
#if (CUDART_VERSION >= 9000)
        if ( tid < 32 )
        {
            mySum = sdata[tid] + sdata[tid + 32];
            for (int offset = warpSize/2; offset > 0; offset /= 2) {
                mySum += __shfl_down_sync(0xFFFFFFFF, mySum, offset,32);
            }
        }
#else
        if (tid < 32) {
            warpReduce(sdata, tid);
            mySum = sdata[0];
        }
#endif
        if (tid == 0) g_odata[blockIdx.x] = mySum;
    }
    
    
    
    
// main function
void aw_pocs_tv(float* img,float* dst,float alpha,const long* image_size, int maxIter,const float delta, const GpuIds& gpuids){
        // Prepare for MultiGPU
        int deviceCount = gpuids.GetLength();
        cudaCheckErrors("Device query fail");
        if (deviceCount == 0) {
            mexErrMsgIdAndTxt("minimizeAwTV:GD_AwTV:GPUselect","There are no available device(s) that support CUDA\n");
        }
        //
        // CODE assumes
        // 1.-All available devices are usable by this code
        // 2.-All available devices are equal, they are the same machine (warning thrown)
        // Check the available devices, and if they are the same
        if (!gpuids.AreEqualDevices()) {
            mexWarnMsgIdAndTxt("minimizeAwTV:GD_AwTV:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
        }
        // onError releases what the call holds if it is left by an error (see
        // tigre_release_on_error.hpp). Without it, an error such as out of memory
        // left the caller's arrays registered with cudaHostRegister; a later array
        // at the same address then used the old memory pages, and calls returned
        // zeros without any error.
        tigre::ReleaseOnError onError;
        std::vector<int> devIds;
        for (int d = 0; d < deviceCount; d++){ devIds.push_back(gpuids[d]); onError.addDevice(gpuids[d]); }

        int dev;
        
        // We don't know if the devices are being used. lets check that. and only use the amount of memory we need.
        // check free memory
        size_t mem_GPU_global;
        checkFreeMemory(gpuids, &mem_GPU_global);

        
        
        // %5 of free memory should be enough, we have almost no variables in these kernels
        size_t total_pixels              = (size_t)image_size[0]*image_size[1]  * image_size[2] ;
        size_t mem_slice_image           = sizeof(float)* (size_t)image_size[0]*image_size[1]  ;
        size_t mem_size_image            = sizeof(float)* total_pixels;
        size_t mem_auxiliary             = sizeof(float)* (total_pixels + MAXTHREADS - 1) / MAXTHREADS;
        
        // Decide how are we handling the distribution of computation
        size_t mem_img_each_GPU;
        
        unsigned int buffer_length=2;
        //Does everything fit in the GPU?
        unsigned int slices_per_split;
        
        // if it is a thin problem (no need to split), just use one GPU
        if (image_size[2]<4){deviceCount=1;}
        
        unsigned int splits=1; // if the number does not fit in an uint, you have more serious trouble than this.
        if(mem_GPU_global> 3*mem_size_image+3*(deviceCount-1)*mem_slice_image*buffer_length+mem_auxiliary) {
            // We only need to split if we have extra GPUs
            slices_per_split=(image_size[2]+deviceCount-1)/deviceCount;
            // No empty pieces (see tvPiecesWithoutEmpty); this may use fewer GPUs or smaller pieces.
            tvPiecesWithoutEmpty((unsigned long)image_size[2], &deviceCount, &splits, &slices_per_split);
            mem_img_each_GPU=mem_slice_image*((slices_per_split+buffer_length*2));
        }else{
            // The smallest working set: on each GPU, 3 images of one slice plus the
            // minimum overlap (2 slices on each side), and the auxiliary memory. If
            // not even that fits, say so now: the sizes below are unsigned and would
            // wrap around past zero.
            const size_t mem_min=3*mem_slice_image*(1+2*2)+mem_auxiliary;
            if (mem_GPU_global<=mem_min)
                tigreNotEnoughGpuMemory("minimizeAwTV:GD_AwTV:GPUmemory","Adaptive-weighted total variation minimization (AwminTV)",mem_min,mem_GPU_global);
            // As mem_auxiliary is not expected to be a large value (for a 2000^3 image is around 28Mbytes), lets for now assume we need it all
            size_t mem_free=mem_GPU_global-mem_auxiliary;
            
            splits=(unsigned int)(ceil(((float)(3*mem_size_image)/(float)(deviceCount))/mem_free));
            // Now, there is an overhead here, as each splits should have 2 slices more, to account for overlap of images.
            // lets make sure these 2 slices fit, if they do not, add 1 to splits.
            slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits);
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // if the new stuff does not fit in the GPU, it means we are in the edge case where adding that extra slice will overflow memory
            // One more split usually does it; with very little memory it may take more
            // (one slice per piece always fits, checked above).
            while (mem_GPU_global< 3*mem_img_each_GPU+mem_auxiliary && slices_per_split>1){
                splits++;
                //recompute for later
                slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits); // amount of slices that fit on a GPU. Later we add 2 to these, as we need them for overlap
                mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            }
            // No empty pieces (see tvPiecesWithoutEmpty); this may use fewer GPUs or smaller pieces.
            tvPiecesWithoutEmpty((unsigned long)image_size[2], &deviceCount, &splits, &slices_per_split);
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));


            // How many EXTRA buffer slices should be able to fit in here??!?!
            // Only do it if there are splits needed. 
            if(splits>1){
                mem_free=mem_GPU_global-(3*mem_img_each_GPU+mem_auxiliary);
                unsigned int extra_buff=(mem_free/mem_slice_image); 
                buffer_length=(extra_buff/2)/3; // we need double whatever this results in, rounded down.
                buffer_length=max(buffer_length,2);// minimum 2
                buffer_length=min(MAX_BUFFER,buffer_length);

                mem_img_each_GPU=mem_slice_image*(slices_per_split+buffer_length*2);
                
            }else{
                buffer_length=2;
            }

            // Assert
            if (mem_GPU_global< 3*mem_img_each_GPU+mem_auxiliary){
                mexErrMsgIdAndTxt("minimizeAwTV:GD_AwTV:GPU","Assertion Failed. Logic behind splitting flawed! Please tell: ander.biguri@gmail.com\n");
            }
        }
        
        
         // Assert
       
        if ((slices_per_split+buffer_length*2)*(size_t)image_size[0]*image_size[1]* sizeof(float)!= mem_img_each_GPU){
            mexErrMsgIdAndTxt("minimizeAwTV:GD_AwTV:GPU","Assertion Failed. Memory needed calculation broken! Please tell: ander.biguri@gmail.com\n");
        }
        
        
        
        
        
        
        float** d_image=    (float**)calloc(deviceCount,sizeof(float*));
        float** d_dimgTV=   (float**)calloc(deviceCount,sizeof(float*));
        float** d_norm2aux= (float**)calloc(deviceCount,sizeof(float*));
        float** d_norm2=    (float**)calloc(deviceCount,sizeof(float*));
        onError.add([d_image, d_dimgTV, d_norm2aux, d_norm2, deviceCount, devIds]{
            for (int d = 0; d < deviceCount; d++){
                cudaSetDevice(devIds[d]);
                cudaFree(d_image[d]); cudaFree(d_dimgTV[d]); cudaFree(d_norm2aux[d]); cudaFree(d_norm2[d]);
            }
            free(d_image); free(d_dimgTV); free(d_norm2aux); free(d_norm2);
        });
         
        // allocate memory in each GPU
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            
            tigre::checkMalloc(cudaMalloc((void**)&d_image[dev], mem_img_each_GPU), "minimizeAwTV:GD_AwTV:GPU", "Adaptive-weighted total variation minimization (AwminTV)", mem_img_each_GPU);
            cudaMemset(         d_image[dev],0  , mem_img_each_GPU); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            tigre::checkMalloc(cudaMalloc((void**)&d_dimgTV[dev], mem_img_each_GPU), "minimizeAwTV:GD_AwTV:GPU", "Adaptive-weighted total variation minimization (AwminTV)", mem_img_each_GPU);
            cudaMemset(         d_dimgTV[dev],0 , mem_img_each_GPU); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            tigre::checkMalloc(cudaMalloc((void**)&d_norm2[dev], slices_per_split*mem_slice_image), "minimizeAwTV:GD_AwTV:GPU", "Adaptive-weighted total variation minimization (AwminTV)", slices_per_split*mem_slice_image);
            cudaMemset(         d_norm2[dev],0  , slices_per_split*mem_slice_image); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            tigre::checkMalloc(cudaMalloc((void**)&d_norm2aux[dev], mem_auxiliary), "minimizeAwTV:GD_AwTV:GPU", "Adaptive-weighted total variation minimization (AwminTV)", mem_auxiliary);
            cudaMemset(         d_norm2aux[dev],0 , mem_auxiliary); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            cudaCheckErrors("Malloc  error");
            
            
        }
       unsigned long long buffer_pixels=buffer_length*(size_t)image_size[0]*image_size[1];
        // Halo exchange e=2*dev copies the last slices of dev to dev+1, e=2*dev+1 the first slices of dev+1 to dev.
        // Neighbours with peer access (peer_access[dev]) copy directly; the others go through their own part of
        // the pinned buffer, so all of them can run at once.
        int nExchanges=2*(deviceCount-1);
        bool* peer_access=(bool*)malloc((deviceCount-1)*sizeof(bool));
        float* buffer=NULL;
        cudaEvent_t* exchange_event=(cudaEvent_t*)calloc(nExchanges>0? nExchanges : 1,sizeof(cudaEvent_t));
        if(splits>1){
            mexWarnMsgIdAndTxt("minimizeAwTV:GD_AwTV:Image_split","Your image can not be fully split between the available GPUs. The computation of minTV will be significantly slowed due to the image size.\nApproximated mathematics turned on for computational speed.");
        }else if(deviceCount>1){
            cudaMallocHost((void**)&buffer,nExchanges*buffer_pixels*sizeof(float));
            onError.add([buffer]{ cudaFreeHost(buffer); });
            for (dev = 0; dev < deviceCount-1; dev++){
                peer_access[dev]=UsePeerHaloExchange(gpuids[dev],gpuids[dev+1]);
                cudaSetDevice(gpuids[dev]);
                cudaEventCreateWithFlags(&exchange_event[2*dev],cudaEventDisableTiming);
                cudaSetDevice(gpuids[dev+1]);
                cudaEventCreateWithFlags(&exchange_event[2*dev+1],cudaEventDisableTiming);
            }
        }
        
        
        
        // Lets try to make the host memory pinned:
        // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
        onError.add([exchange_event, nExchanges, peer_access]{
            for (int e = 0; e < nExchanges; e++) if (exchange_event[e]) cudaEventDestroy(exchange_event[e]);
            free(exchange_event); free(peer_access);
        });
        int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
        cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
        // splits>2 is completely empirical observation
        // Only what this call registers is unregistered, on error or at the end.
        bool img_registered=false, dst_registered=false;
        if (isHostRegisterSupported & splits>2){
            img_registered=tigre::hostRegister(img ,(size_t)image_size[2]*image_size[1]*image_size[0]*sizeof(float));
            dst_registered=tigre::hostRegister(dst ,(size_t)image_size[2]*image_size[1]*image_size[0]*sizeof(float));
        }
        if (img_registered) onError.add([img]{ cudaHostUnregister(img); });
        if (dst_registered) onError.add([dst]{ cudaHostUnregister(dst); });
        cudaCheckErrors("Error pinning memory");

        
        
                // Create streams
        int nStream_device=2;
        int nStreams=deviceCount*nStream_device;
        cudaStream_t* stream=(cudaStream_t*)calloc(nStreams,sizeof(cudaStream_t));
        onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
        
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            for (int i = 0; i < nStream_device; ++i){
                cudaStreamCreate(&stream[i+dev*nStream_device]);
            }
        }
        cudaCheckErrors("Stream creation fail");

        
        // For the reduction

        double totalsum_prev;
        double totalsum;
        double sum_curr_spl;
        double * sumnorm2;
        cudaMallocHost((void**)&sumnorm2,deviceCount*sizeof(double));
        onError.add([sumnorm2]{ cudaFreeHost(sumnorm2); });
        
        unsigned int curr_slices;
        unsigned long long curr_pixels;
        size_t linear_idx_start;
        unsigned long long* offset_device=(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* offset_host  =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* bytes_device =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        onError.add([offset_device, offset_host, bytes_device]{ free(offset_device); free(offset_host); free(bytes_device); });
        bool is_first_chunk;
        bool is_last_chunk;
        for(unsigned int i=0;i<maxIter;i+=(buffer_length-1)){
            if(splits>1){
                totalsum_prev=0;
            }
            for(unsigned int sp=0;sp<splits;sp++){
                
                // For each iteration we need to compute all the image. The ordering of these loops
                // need to be like this due to the bounding layers between splits. If more than 1 split is needed
                // for each GPU then there is no other way that taking the entire memory out of GPU and putting it back.
                // If the memory can be shared between GPUs fully without extra splits, then there is an easy way of synchronizing the memory
                
                // Copy image to memory
                for (dev = 0; dev < deviceCount; dev++){
                    curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                    curr_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                    linear_idx_start=(size_t)image_size[0]*image_size[1]*slices_per_split*(sp*deviceCount+dev);
                    
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
                        
                        cudaMemcpyAsync(d_image[dev]+offset_device[dev], img+offset_host[dev]  , bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                        
                        
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                    }
                }
                // if we need to split and its not the first iteration, then we need to copy from Host memory the previosu result.
                if (splits>1 & i>0){
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaMemcpyAsync(d_image[dev]+offset_device[dev], dst+offset_host[dev]  , bytes_device[dev]*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                        
                        
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                    }
                }
                cudaCheckErrors("Memcpy failure on multi split");
                
                for(unsigned int ib=0;  (ib<(buffer_length-1)) && ((i+ib)<maxIter);  ib++){
                    
                    // For the gradient
                    
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        // Compute the gradient of the TV norm
                        
                        // The image's slices in this buffer are [zlo,zhi): the halos of the first and last piece of the
                        // whole image lie beyond it.
                        const long piece=(long)sp*deviceCount+dev;
                        const long depth=(long)(curr_slices+buffer_length*2-1);
                        const long zlo=(piece==0)? (long)buffer_length : 0;
                        const long zhi=(piece==(long)deviceCount*splits-1)? (long)(buffer_length+curr_slices) : depth;
                        // Grid for this piece's own depth (it was sized from whichever piece was computed last, so a
                        // longer piece could miss its last slices).
                        dim3 blockGrad(10, 10, 10);
                        dim3 gridGrad((image_size[0]+blockGrad.x-1)/blockGrad.x, (image_size[1]+blockGrad.y-1)/blockGrad.y, (depth+blockGrad.z-1)/blockGrad.z);
                        gradientTV<<<gridGrad, blockGrad,0,stream[dev*nStream_device]>>>(d_image[dev],d_dimgTV[dev],depth,zlo,zhi, image_size[1],image_size[0],delta);
                        
                    }
                    
                    
                    
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        // no need to copy the 2 aux slices here
                        cudaStreamSynchronize(stream[dev*nStream_device]);
                        cudaMemcpyAsync(d_norm2[dev], d_dimgTV[dev]+buffer_pixels, (size_t)image_size[0]*image_size[1]*curr_slices*sizeof(float), cudaMemcpyDeviceToDevice,stream[dev*nStream_device+1]);
                    }
                    
                    
                    // Compute the L2 norm of the gradient. For that, reduction is used.
                    //REDUCE
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        
                        size_t dimblockRed = MAXTHREADS;
                        size_t dimgridRed = (total_pixels + MAXTHREADS - 1) / MAXTHREADS;
                        
                        cudaStreamSynchronize(stream[dev*nStream_device+1]);
                        reduceNorm2 << <dimgridRed, dimblockRed, MAXTHREADS*sizeof(float),stream[dev*nStream_device]>> >(d_norm2[dev], d_norm2aux[dev], total_pixels);
                        
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        size_t dimblockRed = MAXTHREADS;
                        size_t dimgridRed = (total_pixels + MAXTHREADS - 1) / MAXTHREADS;

                        if (dimgridRed > 1) {
                            // d_norm2 is free here (reduceNorm2 has read it) and holds more than 1024 floats: room for the double
                            reduceSum << <1, dimblockRed, MAXTHREADS*sizeof(double),stream[dev*nStream_device] >> >(d_norm2aux[dev], (double*)d_norm2[dev], dimgridRed);
                            cudaStreamSynchronize(stream[dev*nStream_device]);
                            cudaMemcpyAsync(&sumnorm2[dev], d_norm2[dev], sizeof(double), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                        }
                        else {
                            cudaStreamSynchronize(stream[dev*nStream_device]);
                            float onesum;
                            cudaMemcpy(&onesum, d_norm2aux[dev], sizeof(float), cudaMemcpyDeviceToHost);
                            sumnorm2[dev]=onesum;
                        }
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                     }
                    cudaCheckErrors("Reduction error");
                    
                    
                    // Accumulate the norm accross devices
                    sum_curr_spl=0;
                    // this is CPU code
                    for (dev = 0; dev < deviceCount; dev++){
                        sum_curr_spl+=sumnorm2[dev];
                    }
                    sum_curr_spl+=0.0000001f; // avoid division by zero
                    
                    // If we have more than one splits, lets use the result from prior calls
                    if(i>0 && splits>1){
                        // this is already stored:
                        //totalsum=totalsum_prev; 
                    }else{
                        totalsum=sum_curr_spl;
                    }
                    
                    
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        //NORMALIZE
                        //in a Tesla, maximum blocks =15 SM * 4 blocks/SM
                        divideArrayScalar  <<<60,MAXTHREADS,0,stream[dev*nStream_device]>>>(d_dimgTV[dev]+buffer_pixels,(float)sqrt(totalsum),total_pixels);
                        //MULTIPLY HYPERPARAMETER
                        multiplyArrayScalar<<<60,MAXTHREADS,0,stream[dev*nStream_device]>>>(d_dimgTV[dev]+buffer_pixels,alpha,   total_pixels);
                    }
                     for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                     }
                    cudaCheckErrors("Scalar operations error");
                    
                    //SUBSTRACT GRADIENT
                    //////////////////////////////////////////////
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        
                        substractArrays<<<60,MAXTHREADS,0,stream[dev*nStream_device]>>>(d_image[dev]+buffer_pixels,d_dimgTV[dev]+buffer_pixels, total_pixels);
                    }
                }

                // Synchronize mathematics, make sure bounding pixels are correct
                 for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                     }
                
                if(splits==1){
                    // Each device moves its last slices on stream 1 and its first slices on stream 0. Queue all downloads,
                    // then the uploads, each waiting for its own download; the synchronization below waits for all of them.
                    for(dev=0; dev<deviceCount-1;dev++){
                        if (peer_access[dev]) continue;
                        int to_next=2*dev, to_prev=2*dev+1;
                        cudaSetDevice(gpuids[dev]);
                        cudaMemcpyAsync(buffer+to_next*buffer_pixels, d_image[dev]+(size_t)slices_per_split*(size_t)image_size[0]*image_size[1], buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                        cudaEventRecord(exchange_event[to_next],stream[dev*nStream_device+1]);
                        cudaSetDevice(gpuids[dev+1]);
                        cudaMemcpyAsync(buffer+to_prev*buffer_pixels, d_image[dev+1]+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device]);
                        cudaEventRecord(exchange_event[to_prev],stream[(dev+1)*nStream_device]);
                    }
                    for(dev=0; dev<deviceCount-1;dev++){
                        int to_next=2*dev, to_prev=2*dev+1;
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        if (peer_access[dev]){
                            cudaSetDevice(gpuids[dev]);
                            cudaMemcpyPeerAsync(d_image[dev+1], gpuids[dev+1], d_image[dev]+(size_t)slices_per_split*(size_t)image_size[0]*image_size[1], gpuids[dev], buffer_pixels*sizeof(float), stream[dev*nStream_device+1]);
                            cudaSetDevice(gpuids[dev+1]);
                            cudaMemcpyPeerAsync(d_image[dev]+total_pixels+buffer_pixels, gpuids[dev], d_image[dev+1]+buffer_pixels, gpuids[dev+1], buffer_pixels*sizeof(float), stream[(dev+1)*nStream_device]);
                            continue;
                        }
                        cudaSetDevice(gpuids[dev+1]);
                        cudaStreamWaitEvent(stream[(dev+1)*nStream_device],exchange_event[to_next],0);
                        cudaMemcpyAsync(d_image[dev+1],buffer+to_next*buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev+1)*nStream_device]);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamWaitEvent(stream[dev*nStream_device+1],exchange_event[to_prev],0);
                        cudaMemcpyAsync(d_image[dev]+total_pixels+buffer_pixels,buffer+to_prev*buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);
                    }
                }else{
                    
                    // We need to take it out :(
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(gpuids[dev]);
                        
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        linear_idx_start=(size_t)image_size[0]*image_size[1]*slices_per_split*(sp*deviceCount+dev);
                        total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                        cudaMemcpyAsync(&dst[linear_idx_start], d_image[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                    }
                }
                
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(gpuids[dev]);
                    cudaDeviceSynchronize();
                }
                cudaCheckErrors("Memory gather error");

                totalsum_prev+=sum_curr_spl;
            }
            totalsum=totalsum_prev;
        }
        // If there has not been splits, we still have data in memory
        if(splits==1){
            for(dev=0; dev<deviceCount;dev++){
                cudaSetDevice(gpuids[dev]);
                
                curr_slices=((dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*dev;
                total_pixels=curr_slices*(size_t)image_size[0]*image_size[1];
                cudaMemcpy(dst+slices_per_split*(size_t)image_size[0]*image_size[1]*dev, d_image[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost);
            }
        }
        cudaCheckErrors("Copy result back");
        onError.dismiss();  // from here on the function releases everything itself
        
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaFree(d_image[dev]);
            cudaFree(d_norm2aux[dev]);
            cudaFree(d_dimgTV[dev]);
            cudaFree(d_norm2[dev]);
        }
        if (splits==1 && deviceCount>1){
            cudaFreeHost(buffer);
            for (int e = 0; e < nExchanges; e++){
                cudaEventDestroy(exchange_event[e]);
            }
        }
        free(exchange_event);
        free(peer_access);
        
        if (img_registered) cudaHostUnregister(img);
        if (dst_registered) cudaHostUnregister(dst);
        for (int i = 0; i < nStreams; ++i)
           cudaStreamDestroy(stream[i]) ;
        cudaCheckErrors("Memory free");
//         cudaDeviceReset();
    }
        
void checkFreeMemory(const GpuIds& gpuids, size_t *mem_GPU_global){
        size_t memfree;
        size_t memtotal;
        const int deviceCount = gpuids.GetLength();
        for (int dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaMemGetInfo(&memfree,&memtotal);
            if(dev==0) *mem_GPU_global=memfree;
            if(memfree<memtotal/2){
                mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
            }
            cudaCheckErrors("Check mem error");
            
            *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
        }
        *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
        
        //*mem_GPU_global= insert your known number here, in bytes.
}
