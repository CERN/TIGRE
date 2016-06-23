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

// Code by:
// Imanol Luengo
// PhD student University of Nottingham
// imaluengo@gmail.com
// 2015
// Sligtly modified by Ander Biguri

// http://gpu4vision.icg.tugraz.at/papers/2010/knoll.pdf#pub47
#define MAXTREADS 1024

#include "tvdenoising.hpp"
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
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
                 const float* spacing, const long* image_size, int maxIter)
{
    // Init params
    size_t total_pixels = image_size[0] * image_size[1]  * image_size[2] ;
    size_t mem_size = sizeof(float) * total_pixels;



    float *d_src, *d_u, *d_px, *d_py, *d_pz;

    // F
    cudaMalloc(&d_src, mem_size);
    cudaMemcpy(d_src, src, mem_size, cudaMemcpyHostToDevice);
    cudaCheckErrors("Memory Malloc and Memset: SRC");

    

    // U
    cudaMalloc(&d_u, mem_size);
    cudaMemcpy(d_u, d_src, mem_size, cudaMemcpyDeviceToDevice);
    cudaCheckErrors("Memory Malloc and Memset: U");
    // PX
    cudaMalloc(&d_px, mem_size);
    cudaMemset(d_px, 0, mem_size);
    cudaCheckErrors("Memory Malloc and Memset: PX");
    // PY
    cudaMalloc(&d_py, mem_size);
    cudaMemset(d_py, 0, mem_size);
    cudaCheckErrors("Memory Malloc and Memset: PY");
    // PZ
    cudaMalloc(&d_pz, mem_size);
    cudaMemset(d_pz, 0, mem_size);
    cudaCheckErrors("Memory Malloc and Memset: PZ");

    // bdim and gdim
    dim3 block(10, 10, 10);
    dim3 grid((image_size[0]+block.x-1)/block.x, (image_size[1]+block.y-1)/block.y, (image_size[2]+block.z-1)/block.z);

    int i = 0;

    float tau2, tau1;
    
    for ( i = 0; i < maxIter; i++ )
    {
        tau2 = 0.3f + 0.02f * i;
        tau1 = (1.f/tau2) * ((1.f/6.f) - (5.f/(15.f+i)));

        update_u<<<grid, block>>>(d_src, d_pz, d_py, d_px, d_u, tau1, lambda,
                                  image_size[2], image_size[1],image_size[0],
                                  spacing[2], spacing[1], spacing[0]);

        update_p<<<grid, block>>>(d_u, d_pz, d_py, d_px, tau2,
                                  image_size[2], image_size[1], image_size[0],
                                  spacing[2], spacing[1], spacing[0]);
        
    }

    cudaCheckErrors("TV minimization"); 

    cudaMemcpy(dst, d_u, mem_size, cudaMemcpyDeviceToHost);
    cudaCheckErrors("Copy result back");

    cudaFree(d_src);
    cudaFree(d_u);
    cudaFree(d_pz);
    cudaFree(d_py);
    cudaFree(d_px);
    //cudaDeviceReset();
}