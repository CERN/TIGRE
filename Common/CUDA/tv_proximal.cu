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

#include "tv_proximal.hpp"
#include "gpuUtils.hpp"
#include "tigre_release_on_error.hpp"
#include "page_locking.hpp"
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                TIGRE_MATLAB_ONLY(cudaDeviceReset());\
                tigreCudaError("CBCT:CUDA:TVdenoising", msg, cudaGetErrorString(__err));\
        } \
} while (0)
void cpy_from_host(float* device_array,float* host_array, 
                   unsigned long long bytes_device,unsigned long long offset_device,unsigned long long offset_host, 
                   unsigned long long pixels_per_slice, unsigned int buffer_length, 
                   cudaStream_t stream, bool is_first_chunk, bool is_last_chunk,const long* image_size);  
    
    
    __global__ void multiplyArrayScalar(float* vec,float scalar,const size_t n)
    {
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]*=scalar;
        }
    }

    __device__ __inline__
            float divergence(const float* pz, const float* py, const float* px,
            long z, long y, long x, long depth, long rows, long cols,
            float dz, float dy, float dx)
    {
        long size2d = rows*cols;
        size_t idx = (size_t)z * size2d + y * cols + x;
        float _div = 0.0f;
        
        if ( z - 1 >= 0 ) {
            _div += (pz[idx] - pz[(size_t)(z-1)*size2d + y*cols + x]) / dz;
        } else {
            _div += pz[idx];
        }
        
        if ( y - 1 >= 0 ) {
            _div += (py[idx] - py[(size_t)z*size2d + (y-1)*cols + x]) / dy;
        } else {
            _div += py[idx];
        }
        
        if ( x - 1 >= 0 ) {
            _div += (px[idx] - px[(size_t)z*size2d + y*cols + (x-1)]) / dx;
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
        size_t idx = (size_t)z * size2d + y * cols + x;
        
        float uidx = u[idx];
        
        if ( z + 1 < depth ) {
            grad[0] = (u[(size_t)(z+1)*size2d + y*cols + x] - uidx) / dz;
        }
        
        if ( y + 1 < rows ) {
            grad[1] = (u[(size_t)z*size2d + (y+1)*cols + x] - uidx) / dy;
        }
        
        if ( x + 1 < cols ) {
            grad[2] = (u[(size_t)z*size2d + y*cols + (x+1)] - uidx) / dx;
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
        size_t idx = (size_t)z * (rows * cols) + y * cols + x;
        
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
        size_t idx = (size_t)z * (rows * cols) + y * cols + x;
        
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
        
        size_t mem_GPU_global;
        checkFreeMemory(gpuids, &mem_GPU_global);
        
        
        // %5 of free memory should be enough, we have almost no variables in these kernels
       size_t total_pixels           = (size_t)image_size[0]*image_size[1] * image_size[2] ;
       const size_t pixels_per_slice = (size_t)image_size[0]*image_size[1] ;
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
            // No empty pieces (see tvPiecesWithoutEmpty); this may use fewer GPUs or smaller pieces.
            tvPiecesWithoutEmpty((unsigned long)image_size[2], &deviceCount, &splits, &slices_per_split);
            mem_img_each_GPU=mem_slice_image*(slices_per_split + buffer_length*2);
        }else{
            // The smallest working set: on each GPU, 5 images of one slice plus the
            // minimum overlap (1 slice on each side). If not even that fits, say so
            // now: the sizes below are unsigned and would wrap around past zero.
            const size_t mem_min=5*mem_slice_image*(1+2*1);
            if (mem_GPU_global<=mem_min)
                tigreNotEnoughGpuMemory("tvDenoise:tvdenoising:GPUmemory","Total variation denoising (tvdenoise)",mem_min,mem_GPU_global);
            // As mem_auxiliary is not expected to be a large value (for a 2000^3 image is around 28Mbytes), lets for now assume we need it all
            size_t mem_free=mem_GPU_global;
            
            splits=(unsigned int)(ceil(((float)(5*mem_size_image)/(float)(deviceCount))/mem_free));
            // Now, there is an overhead here, as each splits should have 2 slices more, to account for overlap of images.
            // lets make sure these 2 slices fit, if they do not, add 1 to splits.
            slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits);
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // if the new stuff does not fit in the GPU, it means we are in the edge case where adding that extra slice will overflow memory
            // One more split usually does it; with very little memory it may take more
            // (one slice per piece always fits, checked above).
            while (mem_GPU_global< 5*mem_img_each_GPU && slices_per_split>1){
                splits++;
                //recompute for later
                slices_per_split=(image_size[2]+deviceCount*splits-1)/(deviceCount*splits); // amount of slices that fit on a GPU. Later we add 2 to these, as we need them for overlap
                mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            }
            // No empty pieces (see tvPiecesWithoutEmpty); this may use fewer GPUs or smaller pieces.
            tvPiecesWithoutEmpty((unsigned long)image_size[2], &deviceCount, &splits, &slices_per_split);
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // How many EXTRA buffer slices should be able to fit in here??!?!
            mem_free=mem_GPU_global-(5*mem_img_each_GPU);
            unsigned int extra_buff=(mem_free/mem_slice_image);
            buffer_length=(extra_buff/2)/5; // we need double whatever this results in, rounded down.
            // At least 1: the main loop advances by buffer_length slices, so 0 never ended.
            if (buffer_length<1) buffer_length=1;
            buffer_length=min(MAX_BUFFER,buffer_length);
            // No wider than a piece: a piece's halo then comes from its neighbour piece only (see "carry" below).
            buffer_length=min(slices_per_split,buffer_length);
            
            mem_img_each_GPU=(mem_slice_image*(slices_per_split+buffer_length*2));
            
            // Assert
            if (mem_GPU_global< 5*mem_img_each_GPU){
                mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","Bad assert. Logic behind splitting flawed! Please tell: ander.biguri@gmail.com\n");
            }
        }
        
        
        // Lets try to make the host memory pinned:
        // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
        int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
        cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
        // Only what this call registers is unregistered, on error or at the end.
        bool src_registered=false, dst_registered=false;
        if (isHostRegisterSupported & splits>1){
            src_registered=tigre::hostRegister(src ,(size_t)image_size[2]*image_size[1]*image_size[0]*sizeof(float));
            dst_registered=tigre::hostRegister(dst ,(size_t)image_size[2]*image_size[1]*image_size[0]*sizeof(float));
        }
        if (src_registered) onError.add([src]{ cudaHostUnregister(src); });
        if (dst_registered) onError.add([dst]{ cudaHostUnregister(dst); });
        cudaCheckErrors("Error pinning memory");
        
        
        
        // Lets allocate auxiliary  variables.
        // Halo exchange e=2*dev copies the last buffer_length slices of dev to dev+1, e=2*dev+1 the first buffer_length
        // slices of dev+1 to dev. Neighbours with peer access (peer_access[dev]) copy directly; for the others each
        // variable v (u, px, py, pz) goes through its own part 4*e+v of the pinned buffer, buffer_length slices each,
        // so all of them can run at once. (Each part used to be one slice: with buffer_length>1 the parts overlapped
        // and the last one wrote past the end of the buffer, #891.)
        int nExchanges=2*(deviceCount-1);
        bool* peer_access=(bool*)malloc((deviceCount-1)*sizeof(bool));
        float* buffer=NULL;
        cudaEvent_t* exchange_event=(cudaEvent_t*)calloc(nExchanges>0? 4*nExchanges : 1,sizeof(cudaEvent_t));
        float* h_px=NULL, *h_py=NULL, *h_pz=NULL, *h_u;
        // With splits>1 the pieces run one after another and write their results back to the host. The lower halo
        // of the first piece of a split is the top of the last piece of the previous split, which has already been
        // advanced by then. "carry" keeps those slices (u, px, py, pz) as they were before that piece ran, so every
        // piece starts from neighbours at its own iteration count.
        float* carry=NULL;
        onError.add([exchange_event, nExchanges, peer_access]{
            for (int e = 0; e < 4*nExchanges; e++) if (exchange_event[e]) cudaEventDestroy(exchange_event[e]);
            free(exchange_event); free(peer_access);
        });
        if(splits>1){
            
            //These take A LOT of memory and A LOT of time to use. If we can avoid using them, better.
            if (buffer_length<maxIter){ // if we do only 1 big iter, they are not needed.
                mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:Memory","TV denoising requires 5 times the image memory. Your GPU(s) do not have the required memory.\n This memory will be attempted to allocate on the CPU, which may fail or slow the computation by a very significant amount.\n If you want to kill the execution: CTRL+C");
                
                cudaMallocHost((void**)&h_px,(size_t)image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                onError.add([h_px]{ cudaFreeHost(h_px); });
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
                
                cudaMallocHost((void**)&h_py,(size_t)image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                onError.add([h_py]{ cudaFreeHost(h_py); });
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");
                
                cudaMallocHost((void**)&h_pz,(size_t)image_size[0]*image_size[1]*image_size[2]*sizeof(float));
                onError.add([h_pz]{ cudaFreeHost(h_pz); });
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.\n Your image is too big to use SART_TV or im3Ddenoise in your current machine");

                cudaMallocHost((void**)&carry,4*(size_t)buffer_length*pixels_per_slice*sizeof(float));
                onError.add([carry]{ cudaFreeHost(carry); });
                cudaCheckErrors("Malloc error on auxiliary variables on CPU.");
            }
            h_u=dst;
        }else if(deviceCount>1){
            cudaMallocHost((void**)&buffer, 4*nExchanges*(size_t)buffer_length*pixels_per_slice*sizeof(float));
            onError.add([buffer]{ cudaFreeHost(buffer); });
            for (dev = 0; dev < deviceCount-1; dev++){
                peer_access[dev]=UsePeerHaloExchange(gpuids[dev],gpuids[dev+1]);
                for (int v = 0; v < 4; v++){
                    cudaSetDevice(gpuids[dev]);
                    cudaEventCreateWithFlags(&exchange_event[8*dev+v],cudaEventDisableTiming);
                    cudaSetDevice(gpuids[dev+1]);
                    cudaEventCreateWithFlags(&exchange_event[8*dev+4+v],cudaEventDisableTiming);
                }
            }
        }
        // We should be good to go memory wise.
        
        
        float** d_src   =(float**)calloc(deviceCount,sizeof(float*));
        float** d_u     =(float**)calloc(deviceCount,sizeof(float*));
        float** d_px    =(float**)calloc(deviceCount,sizeof(float*));
        float** d_py    =(float**)calloc(deviceCount,sizeof(float*));
        float** d_pz    =(float**)calloc(deviceCount,sizeof(float*));
        onError.add([d_src, d_u, d_px, d_py, d_pz, deviceCount, devIds]{
            for (int d = 0; d < deviceCount; d++){
                cudaSetDevice(devIds[d]);
                cudaFree(d_src[d]); cudaFree(d_u[d]); cudaFree(d_px[d]); cudaFree(d_py[d]); cudaFree(d_pz[d]);
            }
            free(d_src); free(d_u); free(d_px); free(d_py); free(d_pz);
        });
        
        //Malloc
        for(dev=0;dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            // F
            tigre::checkMalloc(cudaMalloc((void**)&d_src[dev], mem_img_each_GPU), "tvDenoise:tvdenoising:GPU", "Total variation denoising (tvdenoise)", mem_img_each_GPU);
            // U
            tigre::checkMalloc(cudaMalloc((void**)&d_u[dev], mem_img_each_GPU), "tvDenoise:tvdenoising:GPU", "Total variation denoising (tvdenoise)", mem_img_each_GPU);
            // PX
            tigre::checkMalloc(cudaMalloc((void**)&d_px[dev], mem_img_each_GPU), "tvDenoise:tvdenoising:GPU", "Total variation denoising (tvdenoise)", mem_img_each_GPU);
            // PY
            tigre::checkMalloc(cudaMalloc((void**)&d_py[dev], mem_img_each_GPU), "tvDenoise:tvdenoising:GPU", "Total variation denoising (tvdenoise)", mem_img_each_GPU);
            // PZ
            tigre::checkMalloc(cudaMalloc((void**)&d_pz[dev], mem_img_each_GPU), "tvDenoise:tvdenoising:GPU", "Total variation denoising (tvdenoise)", mem_img_each_GPU);
        }
        cudaDeviceSynchronize();
        cudaCheckErrors("Malloc  error");
        
        
        // Create streams
        int nStream_device=5;
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
        
        
        
        
        // Allocate CPU buffer if needed, warn user if not.
        
        
        
        unsigned int curr_slices;
        unsigned long long curr_pixels;
        size_t linear_idx_start;
        
        unsigned long long buffer_pixels=buffer_length*pixels_per_slice;
        
        unsigned long long* offset_device=(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* offset_host  =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        unsigned long long* bytes_device =(unsigned long long*)malloc(deviceCount*sizeof(unsigned long long));
        // Halo slices below and above each piece: buffer_length, or fewer where the piece is closer than that to the
        // image's first or last slice. There the kernels' own boundary conditions apply, in z as in x and y, so the
        // result does not depend on how the image is split. On the device, a piece's interior always starts at
        // buffer_pixels; halo_lo[dev] slices below it and halo_hi[dev] above it are used.
        unsigned int* halo_lo=(unsigned int*)malloc(deviceCount*sizeof(unsigned int));
        unsigned int* halo_hi=(unsigned int*)malloc(deviceCount*sizeof(unsigned int));
        onError.add([offset_device, offset_host, bytes_device, halo_lo, halo_hi]{ free(offset_device); free(offset_host); free(bytes_device); free(halo_lo); free(halo_hi); });
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
                    
                    // lets compute where we start copies and how much. This avoids 3 calls to Memcpy
                    size_t start_slice=(size_t)slices_per_split*(sp*deviceCount+dev);
                    halo_lo[dev]=(unsigned int)min((size_t)buffer_length, start_slice);
                    halo_hi[dev]=(unsigned int)min((size_t)buffer_length, (size_t)image_size[2]-start_slice-curr_slices);
                    offset_device[dev]=(buffer_length-halo_lo[dev])*pixels_per_slice;
                    offset_host[dev]=linear_idx_start-halo_lo[dev]*pixels_per_slice;
                    bytes_device[dev]=curr_pixels+(halo_lo[dev]+halo_hi[dev])*pixels_per_slice;
                }
                // copy data to the GPU if we are just starting
                if(i==0){
                    for (dev = 0; dev < deviceCount; dev++){
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);

                        cudaSetDevice(gpuids[dev]);
                        if (is_last_chunk) {cudaMemsetAsync(d_src[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+1]);}
                        cpy_from_host(d_src[dev],src,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+1],  is_first_chunk,  is_last_chunk, image_size);
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);

                        cudaMemcpyAsync(d_u[dev], d_src[dev], mem_img_each_GPU, cudaMemcpyDeviceToDevice,stream[dev*nStream_device+1]);
                        cudaMemsetAsync(d_px[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                        cudaMemsetAsync(d_py[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                        cudaMemsetAsync(d_pz[dev], 0, mem_img_each_GPU,stream[dev*nStream_device]);
                    }

                    // Sync
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
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+1]);
                        if (is_last_chunk) {cudaMemsetAsync(d_u[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+1]);}
                        cpy_from_host(d_u[dev],h_u,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+1],  is_first_chunk,  is_last_chunk, image_size);
                    }

                    for (dev = 0; dev < deviceCount; dev++){ 
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+2]);
                        if (is_last_chunk) {cudaMemsetAsync(d_px[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+2]);}
                        cpy_from_host(d_px[dev],h_px,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+2],  is_first_chunk,  is_last_chunk, image_size);
                    }
                    for (dev = 0; dev < deviceCount; dev++){ 
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+3]);
                        if (is_last_chunk) {cudaMemsetAsync(d_py[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+3]);}
                        cpy_from_host(d_py[dev],h_py,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+3],  is_first_chunk,  is_last_chunk, image_size);
                    }
                    for (dev = 0; dev < deviceCount; dev++){ 
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+4]);
                        if (is_last_chunk) {cudaMemsetAsync(d_pz[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+4]);}
                        cpy_from_host(d_pz[dev],h_pz,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+4],  is_first_chunk,  is_last_chunk, image_size);
                    }
                    // The lower halo of the first piece of this split: the slices the previous split had before it ran.
                    if (sp>0){
                        float** d_var[4]={d_u, d_px, d_py, d_pz};
                        cudaSetDevice(gpuids[0]);
                        for (int v = 0; v < 4; v++){
                            cudaMemcpyAsync(d_var[v][0], carry+(size_t)v*buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[1+v]);
                        }
                    }
                    for (dev = 0; dev < deviceCount; dev++){ 
                        is_last_chunk=!((sp*deviceCount+dev)<deviceCount*splits-1);
                        is_first_chunk=!(sp*deviceCount+dev);
                        cudaSetDevice(gpuids[dev]);
                        cudaStreamSynchronize(stream[dev*nStream_device+1]);
                        if (is_last_chunk) {cudaMemsetAsync(d_src[dev], 0, mem_img_each_GPU,stream[dev*nStream_device+1]);}
                        cpy_from_host(d_src[dev],src,bytes_device[dev], offset_device[dev],offset_host[dev], pixels_per_slice, buffer_length, stream[dev*nStream_device+1],  is_first_chunk,  is_last_chunk, image_size);
                    }

                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                        cudaCheckErrors("Memcpy failure on multi split");
                    }
                }
                // Before the last piece of this split runs, keep its top buffer_length slices for the next split.
                // (Its interior starts at buffer_pixels and has slices_per_split >= buffer_length slices.)
                if (splits>1 && sp<splits-1 && i>0){ // at i==0 every piece starts from src, no carry needed
                    float** d_var[4]={d_u, d_px, d_py, d_pz};
                    dev=deviceCount-1;
                    cudaSetDevice(gpuids[dev]);
                    for (int v = 0; v < 4; v++){
                        cudaMemcpy(carry+(size_t)v*buffer_pixels, d_var[v][dev]+(size_t)slices_per_split*pixels_per_slice, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost);
                    }
                    cudaCheckErrors("Memcpy failure on multi split");
                }
                
                // Inter interations.
                for(unsigned int ib=0;  (ib<(buffer_length)) && ((i+ib)<maxIter);  ib++){

                    tau2 = 0.3f + 0.02f * (i+ib);
                    tau1 = (1.f/tau2) * ((1.f/6.f) - (5.f/(15.f+(i+ib))));
                    // bdim and gdim
                    
                    // Only the piece and its halo_lo/halo_hi halo slices (see above).
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        size_t off=(size_t)(buffer_length-halo_lo[dev])*pixels_per_slice;
                        long depth=(long)curr_slices+halo_lo[dev]+halo_hi[dev];
                        dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (depth+block.z-1)/block.z);
                        
                        update_u<<<grid, block,0,stream[dev*nStream_device]>>>(d_src[dev]+off, d_pz[dev]+off, d_py[dev]+off, d_px[dev]+off, d_u[dev]+off, tau1, lambda,
                                depth, image_size[1],image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                    for (dev = 0; dev < deviceCount; dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices=((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        size_t off=(size_t)(buffer_length-halo_lo[dev])*pixels_per_slice;
                        long depth=(long)curr_slices+halo_lo[dev]+halo_hi[dev];
                        dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
                        dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (depth+block.z-1)/block.z);
                        
                        update_p<<<grid, block,0,stream[dev*nStream_device]>>>(d_u[dev]+off, d_pz[dev]+off, d_py[dev]+off, d_px[dev]+off, tau2,
                                depth, image_size[1], image_size[0],
                                spacing[2], spacing[1], spacing[0]);
                    }
                }// END internal iter
                
                // Synchronize mathematics, make sure bounding pixels are correct
                for(dev=0; dev<deviceCount;dev++){
                    cudaSetDevice(gpuids[dev]);
                    cudaDeviceSynchronize();
                }

                // We have done as many iterations as our buffer allowed. We now need to syncronize the buffers.
            
                if(splits==1){
                    // If everything fits in the GPUs, we can just share the updates between GPUs directly.
                    // We iterate for each device, and we copy the buffer pixels from each device. 
                    // "buffer" variables are just host auxiliary variables to allow the copy from GPU to GPU. 
                    // Essentially this code takes for each device (exceptions for the first and last devices included) "buffer_pixels" amount 
                    // of the beggining and end of each important variable and passes it to the next/previous GPU (the one containing the next/previous chunck)

                    // Every exchange only reads interior slices and writes halo slices, so they are independent. Queue all
                    // downloads, then the uploads, each waiting for its own download, and synchronize once.
                    float** d_var[4]={d_u, d_px, d_py, d_pz};
                    for(dev=0; dev<deviceCount-1;dev++){
                        for (int v = 0; v < 4; v++){
                            int to_next=8*dev+v, to_prev=8*dev+4+v;
                            if (peer_access[dev]){
                                cudaSetDevice(gpuids[dev]);
                                cudaMemcpyPeerAsync(d_var[v][dev+1], gpuids[dev+1], d_var[v][dev]+slices_per_split*pixels_per_slice, gpuids[dev], buffer_pixels*sizeof(float), stream[dev*nStream_device+1+v]);
                                cudaSetDevice(gpuids[dev+1]);
                                cudaMemcpyPeerAsync(d_var[v][dev]+slices_per_split*pixels_per_slice+buffer_pixels, gpuids[dev], d_var[v][dev+1]+buffer_pixels, gpuids[dev+1], buffer_pixels*sizeof(float), stream[(dev+1)*nStream_device+1+v]);
                                continue;
                            }
                            cudaSetDevice(gpuids[dev]);
                            cudaMemcpyAsync(buffer+to_next*buffer_pixels, d_var[v][dev]+slices_per_split*pixels_per_slice, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1+v]);
                            cudaEventRecord(exchange_event[to_next],stream[dev*nStream_device+1+v]);
                            cudaSetDevice(gpuids[dev+1]);
                            cudaMemcpyAsync(buffer+to_prev*buffer_pixels, d_var[v][dev+1]+buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[(dev+1)*nStream_device+1+v]);
                            cudaEventRecord(exchange_event[to_prev],stream[(dev+1)*nStream_device+1+v]);
                        }
                    }
                    for(dev=0; dev<deviceCount-1;dev++){
                        if (peer_access[dev]) continue;
                        for (int v = 0; v < 4; v++){
                            int to_next=8*dev+v, to_prev=8*dev+4+v;
                            cudaSetDevice(gpuids[dev+1]);
                            cudaStreamWaitEvent(stream[(dev+1)*nStream_device+1+v],exchange_event[to_next],0);
                            cudaMemcpyAsync(d_var[v][dev+1], buffer+to_next*buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[(dev+1)*nStream_device+1+v]);
                            cudaSetDevice(gpuids[dev]);
                            cudaStreamWaitEvent(stream[dev*nStream_device+1+v],exchange_event[to_prev],0);
                            cudaMemcpyAsync(d_var[v][dev]+slices_per_split*pixels_per_slice+buffer_pixels, buffer+to_prev*buffer_pixels, buffer_pixels*sizeof(float), cudaMemcpyHostToDevice,stream[dev*nStream_device+1+v]);
                        }
                    }
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                    }
                // This is the case when we can't solely use GPU memory, as the total size of the images+variables exceeds total amounf of memory among GPUs.
                // This situation requires partial results and full memory allocation in the host. 
                }else{
                    // Vopy all the U variable into the host.
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(gpuids[dev]);
                        curr_slices      = ((sp*deviceCount+dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*(sp*deviceCount+dev);
                        linear_idx_start = pixels_per_slice*slices_per_split*(sp*deviceCount+dev);
                        total_pixels     = curr_slices*pixels_per_slice;
                        cudaMemcpyAsync(&h_u[linear_idx_start],  d_u [dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                    }
                    if ((i+buffer_length)<maxIter){ // If its the last iteration, we don't need to get these out.
                        // if its not, copy them to host fully. 
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
                    // Wait for all devices to finish copying to the host before the next split reads from it:
                    // the first device of the next split reads its halo from what the last device of this split is writing.
                    for(dev=0; dev<deviceCount;dev++){
                        cudaSetDevice(gpuids[dev]);
                        cudaDeviceSynchronize();
                    }

                }
            }//END splits
        }//END main iter
        
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("TV minimization");
        
        // We are done. If we were solely using GPU memory, because the problem fitted fully on all GPU memory available, then the result is still inside the GPU.
        // lets get it out. 
        if(splits==1){
            for(dev=0; dev<deviceCount;dev++){
                cudaSetDevice(gpuids[dev]);
                curr_slices  = ((dev+1)*slices_per_split<image_size[2])?  slices_per_split:  image_size[2]-slices_per_split*dev;
                total_pixels = curr_slices*pixels_per_slice;
                cudaMemcpyAsync(dst+slices_per_split*pixels_per_slice*dev, d_u[dev]+buffer_pixels,total_pixels*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
            }
        } // done, everything in GPU and auxiliary variables are good to go. 

        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("Copy result back");
        onError.dismiss();  // from here on the function releases everything itself
        for(dev=0; dev<deviceCount;dev++){
            cudaSetDevice(gpuids[dev]);
            cudaFree(d_src[dev]);
            cudaFree(d_u [dev]);
            cudaFree(d_pz[dev]);
            cudaFree(d_py[dev]);
            cudaFree(d_px[dev]);
        }
        free(d_src); free(d_u); free(d_px); free(d_py); free(d_pz);  // were freed only on error
        free(offset_device); free(offset_host); free(bytes_device); free(halo_lo); free(halo_hi);
        if(splits>1 && buffer_length<maxIter){
            cudaFreeHost(h_px);
            cudaFreeHost(h_py);
            cudaFreeHost(h_pz);
            cudaFreeHost(carry);
        }else if(splits==1 && deviceCount>1){
            cudaFreeHost(buffer);
            for (int e = 0; e < 4*nExchanges; e++){
                cudaEventDestroy(exchange_event[e]);
            }
        }
        free(exchange_event);
        free(peer_access);
        
        for (int i = 0; i < nStreams; ++i)
           cudaStreamDestroy(stream[i]) ;
        free(stream);

        if (src_registered) cudaHostUnregister(src);
        if (dst_registered) cudaHostUnregister(dst);
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
                mexWarnMsgIdAndTxt("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
            }
            cudaCheckErrors("Check mem error");
            
            *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
        }
        *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
        
        //*mem_GPU_global= insert your known number here, in bytes.
}

void cpy_from_host(float* device_array,float* host_array, 
                   unsigned long long bytes_device,unsigned long long offset_device,unsigned long long offset_host, 
                   unsigned long long pixels_per_slice, unsigned int buffer_length, 
                   cudaStream_t stream, bool is_first_chunk, bool is_last_chunk,const long* image_size)
{

    // The halo slots beyond the image's first and last slice are not used (the kernels skip them and apply
    // their own boundary conditions there), so nothing is copied into them.
    cudaStreamSynchronize(stream);
    cudaMemcpyAsync(device_array +offset_device, host_array +offset_host,  bytes_device*sizeof(float), cudaMemcpyHostToDevice,stream);
}