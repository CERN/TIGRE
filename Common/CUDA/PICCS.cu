/*-------------------------------------------------------------------------
 *
 * CUDA functions for Steepest descend in POCS-type algorithms.
 *
 * This file will iteratively minimize by stepest descend the total variation 
 * of the input image, with the parameters given, using GPUs.
 *
 * CODE by       Ander Biguri
 *
---------------------------------------------------------------------------
---------------------------------------------------------------------------
Copyright (c) 2015, University of Bath and CERN- European Organization for 
Nuclear Research
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
 ---------------------------------------------------------------------------

Contact: tigre.toolbox@gmail.com
Codes  : https://github.com/CERN/TIGRE
--------------------------------------------------------------------------- 
 */







#define MAXTHREADS 1024

#include "PICCS.hpp"




#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("ERROR in: %s \n",msg);\
                mexErrMsgIdAndTxt("err",cudaGetErrorString(__err));\
        } \
} while (0)
    
// CUDA kernels
//https://stackoverflow.com/questions/21332040/simple-cuda-kernel-optimization/21340927#21340927
    __global__ void divideArrayScalar(float* vec,float scalar,const size_t n)
    {
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]/=scalar;
        }
    }
    __global__ void isnan_device(float* vec,const size_t n,bool* result)
    {
        *result = false;
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            if(isnan(vec[i]))
                *result=true;
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
    __global__ void addArrays(float* vec,float* vec2,const size_t n)
    {
        unsigned long long i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]+=vec2[i];
        }
    }
    __device__ __inline__
            void gradient(const float* u, float* grad,
            long z, long y, long x,
            long depth, long rows, long cols)
    {
        unsigned long size2d = rows*cols;
        unsigned long long idx = z * size2d + y * cols + x;
        
        float uidx = u[idx];
        
        if ( z - 1 >= 0 && z<depth) {
            grad[0] = (uidx-u[(z-1)*size2d + y*cols + x]) ;
        }
        
        if ( y - 1 >= 0 && y<rows){
            grad[1] = (uidx-u[z*size2d + (y-1)*cols + x]) ;
        }
        
        if ( x - 1 >= 0 && x<cols) {
            grad[2] = (uidx-u[z*size2d + y*cols + (x-1)]);
        }
    }
    
    __global__ void gradientTV(const float* f, float* dftv,
            long depth, long rows, long cols){
        unsigned long x = threadIdx.x + blockIdx.x * blockDim.x;
        unsigned long y = threadIdx.y + blockIdx.y * blockDim.y;
        unsigned long z = threadIdx.z + blockIdx.z * blockDim.z;
        unsigned long long idx = z * rows * cols + y * cols + x;
        if ( x >= cols || y >= rows || z >= depth )
            return;
        
        float df[3] ={0,0,0};
        float dfi[3]={0,0,0}; // dfi== \partial f_{i+1,j,k}
        float dfj[3]={0,0,0};
        float dfk[3]={0,0,0};
        gradient(f,df  ,z  ,y  ,x  , depth,rows,cols);
        gradient(f,dfi ,z  ,y  ,x+1, depth,rows,cols);
        gradient(f,dfj ,z  ,y+1,x  , depth,rows,cols);
        gradient(f,dfk ,z+1,y  ,x  , depth,rows,cols);
        float eps=0.000001; //% avoid division by zero
        dftv[idx]=(df[0]+df[1]+df[2])/(sqrt(df[0] *df[0] +df[1] *df[1] +df[2] *df[2])+eps)
        -dfi[2]/(sqrt(dfi[0]*dfi[0]+dfi[1]*dfi[1]+dfi[2]*dfi[2]) +eps)     // I wish I coudl precompute this, but if I do then Id need to recompute the gradient.
        -dfj[1]/(sqrt(dfj[0]*dfj[0]+dfj[1]*dfj[1]+dfj[2]*dfj[2]) +eps)
        -dfk[0]/(sqrt(dfk[0]*dfk[0]+dfk[1]*dfk[1]+dfk[2]*dfk[2]) +eps);
        
    }
    
    __device__ void warpReduce(volatile float *sdata, size_t tid) {
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
        size_t i = blockIdx.x*blockDim.x + tid;
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
        
        
#if (__CUDART_VERSION >= 9000)
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
    __global__ void  reduceSum(float *g_idata, float *g_odata, size_t n){
        extern __shared__ volatile float sdata[];
        //http://stackoverflow.com/a/35133396/1485872
        size_t tid = threadIdx.x;
        size_t i = blockIdx.x*blockDim.x + tid;
        size_t gridSize = blockDim.x*gridDim.x;
        float mySum = 0;
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
        
        
#if (__CUDART_VERSION >= 9000)
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
    

bool isnan_cuda(float* vec, size_t size){
    bool*d_nan;
    bool h_nan;
    cudaMalloc((void **)&d_nan, sizeof (bool));
    isnan_device<<<60,MAXTHREADS>>>(vec,size,d_nan);
    cudaMemcpy(&h_nan, d_nan, sizeof(bool), cudaMemcpyDeviceToHost);
    return h_nan;

}
    
// main function
 void piccs_tv(const float* img,const float* prior, float* dst,float alpha,float ratio, const long* image_size, int maxIter, const GpuIds& gpuids){
        
     
        
    
        size_t total_pixels = image_size[0] * image_size[1]  * image_size[2] ;
        size_t mem_size = sizeof(float) * total_pixels;
        
        float *d_image,*d_prior,*d_dpiccsTV, *d_dimgTV,*d_aux_small,*d_aux_image, *d_norm2;
        // memory for image
        cudaMalloc(&d_image, mem_size);
        cudaMalloc(&d_prior, mem_size);

        cudaCheckErrors("Malloc Image error");
        cudaMemcpy(d_image, img, mem_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_prior, prior, mem_size, cudaMemcpyHostToDevice);
        cudaCheckErrors("Memory Malloc and Memset: SRC");
        // memory for df
        cudaMalloc(&d_dimgTV, mem_size);
        cudaMalloc(&d_dpiccsTV, mem_size);
        cudaCheckErrors("Memory Malloc and Memset: TV");
        cudaMalloc(&d_norm2, mem_size);
        cudaCheckErrors("Memory Malloc and Memset: TV");
        cudaMalloc(&d_aux_image, mem_size);
        cudaCheckErrors("Memory Malloc and Memset: TV");
        
        // memory for L2norm auxiliar
        cudaMalloc(&d_aux_small, sizeof(float)*(total_pixels + MAXTHREADS - 1) / MAXTHREADS);
        cudaCheckErrors("Memory Malloc and Memset: NORMAux");
        
        
        
        // For the gradient
        dim3 blockGrad(10, 10, 10);
        dim3 gridGrad((image_size[0]+blockGrad.x-1)/blockGrad.x, (image_size[1]+blockGrad.y-1)/blockGrad.y, (image_size[2]+blockGrad.z-1)/blockGrad.z);
        
        // For the reduction
        float sumnorm2;
        size_t dimblockRed = MAXTHREADS;
        size_t dimgridRed = (total_pixels + MAXTHREADS - 1) / MAXTHREADS;


        for(unsigned int i=0;i<maxIter;i++){
            
            cudaMemcpy( d_aux_image,d_image, mem_size, cudaMemcpyDeviceToDevice);
//             mexPrintf("Iteration %d\n",(int)i);

            // Compute the gradient of the TV norm
            gradientTV<<<gridGrad, blockGrad>>>(d_image,d_dimgTV,image_size[2], image_size[1],image_size[0]);
            cudaDeviceSynchronize();
            cudaCheckErrors("Gradient");
//             mexPrintf("Gradient is nan: %s\n",isnan_cuda(d_dimgTV,total_pixels) ? "true" : "false");


            multiplyArrayScalar<<<60,MAXTHREADS>>>(d_dimgTV,(1-ratio),   total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("Multiplication error");

            substractArrays<<<60,MAXTHREADS>>>(d_aux_image,d_prior, total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("Substraction error");
            
            gradientTV<<<gridGrad, blockGrad>>>(d_aux_image,d_dpiccsTV,image_size[2], image_size[1],image_size[0]);
            cudaDeviceSynchronize();
            cudaCheckErrors("Gradient");
//             mexPrintf("Gradient piccs is nan: %s\n",isnan_cuda(d_dimgTV,total_pixels) ? "true" : "false");

            multiplyArrayScalar<<<60,MAXTHREADS>>>(d_dpiccsTV,ratio,   total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("Multiplication error");
//             mexPrintf("Multiplication is nan: %s\n",isnan_cuda(d_dimgTV,total_pixels) ? "true" : "false");
                
            
            addArrays<<<60,MAXTHREADS>>>(d_dimgTV,d_dpiccsTV,total_pixels);
            cudaDeviceSynchronize();
            //NOMRALIZE via reduction
            //mexPrintf("Pre-norm2 is nan: %s\n",isnan_cuda(d_dimgTV,total_pixels) ? "true" : "false");
            cudaMemcpy(d_norm2, d_dimgTV, mem_size, cudaMemcpyDeviceToDevice);
            cudaCheckErrors("Copy from gradient call error");
            reduceNorm2 << <dimgridRed, dimblockRed, MAXTHREADS*sizeof(float) >> >(d_norm2, d_aux_small, total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("reduce1");
            if (dimgridRed > 1) {
                reduceSum << <1, dimblockRed, MAXTHREADS*sizeof(float) >> >(d_aux_small, d_norm2, dimgridRed);
                cudaDeviceSynchronize();
                cudaCheckErrors("reduce2");
                cudaMemcpy(&sumnorm2, d_norm2, sizeof(float), cudaMemcpyDeviceToHost);
                cudaCheckErrors("cudaMemcpy");

            }
            else {
                cudaMemcpy(&sumnorm2, d_aux_small, sizeof(float), cudaMemcpyDeviceToHost);
                cudaCheckErrors("cudaMemcpy");
            }
//             mexPrintf("alpha/sqrt(sumnorm2): %f\n",alpha/sqrt(sumnorm2));
            //MULTIPLY HYPERPARAMETER sqrt(sumnorm2)
            multiplyArrayScalar<<<60,MAXTHREADS>>>(d_dimgTV,alpha/sqrt(sumnorm2),  total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("Multiplication error");
            //SUBSTRACT GRADIENT
            substractArrays    <<<60,MAXTHREADS>>>(d_image,d_dimgTV, total_pixels);
            cudaDeviceSynchronize();
            cudaCheckErrors("Substraction error");
//             mexPrintf("Final update is nan: %s\n",isnan_cuda(d_image,total_pixels) ? "true" : "false");
//             mexPrintf("\n");
            sumnorm2=0;
        }
        
        cudaCheckErrors("TV minimization");
        
        cudaMemcpy(dst, d_image, mem_size, cudaMemcpyDeviceToHost);
        cudaCheckErrors("Copy result back");
        
        cudaFree(d_image);
        cudaFree(d_dpiccsTV);
        cudaFree(d_aux_image);
        cudaFree(d_aux_small);
        cudaFree(d_prior);
        cudaFree(d_norm2);


        cudaCheckErrors("Memory free");
        cudaDeviceReset();
    }
    
