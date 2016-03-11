#define MAXTHREADS 1024

#include "POCS_TV.hpp"




#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("ERROR in: %s \n",msg);\
                mexErrMsgIdAndTxt("CBCT:CUDA:POCS_TV",cudaGetErrorString(__err));\
        } \
} while (0)
    
// CUDA kernels
//https://stackoverflow.com/questions/21332040/simple-cuda-kernel-optimization/21340927#21340927
    __global__ void divideArrayScalar(float* vec,float scalar,const size_t n)
    {
        unsigned int i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]/=scalar;
        }
    }
    __global__ void multiplyArrayScalar(float* vec,float scalar,const size_t n)
    {
        unsigned int i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]*=scalar;
        }
    }
    __global__ void substractArrays(float* vec,float* vec2,const size_t n)
    {
        unsigned int i = (blockIdx.x * blockDim.x) + threadIdx.x;
        for(; i<n; i+=gridDim.x*blockDim.x) {
            vec[i]-=vec2[i];
        }
    }
    
    __device__ __inline__
            void gradient(const float* u, float* grad,
            long z, long y, long x,
            long depth, long rows, long cols)
    {
        long size2d = rows*cols;
        long idx = z * size2d + y * cols + x;
        
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
        long x = threadIdx.x + blockIdx.x * blockDim.x;
        long y = threadIdx.y + blockIdx.y * blockDim.y;
        long z = threadIdx.z + blockIdx.z * blockDim.z;
        long idx = z * rows * cols + y * cols + x;
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
        float eps=0.00000001; //% avoid division by zero
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
        
        
#if (__CUDA_ARCH__ >= 300)
        if ( tid < 32 )
        {
            mySum = sdata[tid] + sdata[tid + 32];
            for (int offset = warpSize/2; offset > 0; offset /= 2) {
                mySum += __shfl_down(mySum, offset);
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
        float value=0;
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
        
        
#if (__CUDA_ARCH__ >= 300)
        if ( tid < 32 )
        {
            mySum = sdata[tid] + sdata[tid + 32];
            for (int offset = warpSize/2; offset > 0; offset /= 2) {
                mySum += __shfl_down(mySum, offset);
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
    void pocs_tv(const float* img,float* dst,float alpha,const long* image_size, int maxIter){
        
        // Init cuda memory
        // BEFORE DOING ANYTHING: Use the proper CUDA enabled GPU: Tesla K40c or Gforce GT 740M
        int deviceCount = 0;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            mexErrMsgIdAndTxt("cudaGetDeviceCount","No CUDA enabled NVIDIA GPUs found");
        }
        bool found=false;
        for (int dev = 0; dev < deviceCount; ++dev)
        {
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, dev);
            
            if (strcmp(deviceProp.name, "Tesla K40c") == 0 || strcmp(deviceProp.name, "GeForce GT 740M") == 0){
                cudaSetDevice(dev);
                found=true;
                break;
            }
        }
        if (!found)
            mexErrMsgIdAndTxt("cudaDevice","No Supported GPU found");
        
        
        
        size_t total_pixels = image_size[0] * image_size[1]  * image_size[2] ;
        size_t mem_size = sizeof(float) * total_pixels;
        
        float *d_image, *d_dimgTV,*d_norm2aux,*d_norm2;
        // memory for image
        cudaMalloc(&d_image, mem_size);
        cudaMemcpy(d_image, img, mem_size, cudaMemcpyHostToDevice);
        cudaCheckErrors("Memory Malloc and Memset: SRC");
        // memory for df
        cudaMalloc(&d_dimgTV, mem_size);
        cudaCheckErrors("Memory Malloc and Memset: TV");
        
        cudaMalloc(&d_norm2, mem_size);
        cudaCheckErrors("Memory Malloc and Memset: TV");
        
        
        
        
        
        // For the gradient
        dim3 blockGrad(10, 10, 10);
        dim3 gridGrad((image_size[0]+blockGrad.x-1)/blockGrad.x, (image_size[1]+blockGrad.y-1)/blockGrad.y, (image_size[2]+blockGrad.z-1)/blockGrad.z);
        
        // For the reduction
        float sumnorm2;
        
        // memory for L2norm auxiliar
        cudaMalloc(&d_norm2aux, sizeof(float)*(total_pixels + MAXTHREADS - 1) / MAXTHREADS);
        cudaCheckErrors("Memory Malloc and Memset: NORMAux");
        
        for(int i=0;i<maxIter;i++){
            
            
            // Compute the gradient of the TV norm
            gradientTV<<<gridGrad, blockGrad>>>(d_image,d_dimgTV,image_size[2], image_size[1],image_size[0]);
            cudaCheckErrors("Gradient");
//             cudaMemcpy(dst, d_dimgTV, mem_size, cudaMemcpyDeviceToHost);
            
            
            cudaMemcpy(d_norm2, d_dimgTV, mem_size, cudaMemcpyDeviceToDevice);
            
            // Compute the L2 norm of the gradietn. For that, reduction is used.
            //REDUCE
            size_t dimblockRed = MAXTHREADS;
            size_t dimgridRed = (total_pixels + MAXTHREADS - 1) / MAXTHREADS;
            reduceNorm2 << <dimgridRed, dimblockRed, MAXTHREADS*sizeof(float) >> >(d_norm2, d_norm2aux, total_pixels);
            cudaCheckErrors("reduce1");
            if (dimgridRed > 1) {
                reduceSum << <1, dimblockRed, MAXTHREADS*sizeof(float) >> >(d_norm2aux, d_norm2, dimgridRed);
                cudaCheckErrors("reduce2");
                cudaMemcpy(&sumnorm2, d_norm2, sizeof(float), cudaMemcpyDeviceToHost);
                cudaCheckErrors("cudaMemcpy");
                
            }
            else {
                cudaMemcpy(&sumnorm2, d_norm2aux, sizeof(float), cudaMemcpyDeviceToHost);
                cudaCheckErrors("cudaMemcpy");
            }
            //mexPrintf("%f ",sqrt(sumnorm2));
            //NOMRALIZE
            //in a Tesla, maximum blocks =15 SM * 4 blocks/SM
            divideArrayScalar  <<<60,MAXTHREADS>>>(d_dimgTV,sqrt(sumnorm2),total_pixels);
            //MULTIPLY HYPERPARAMETER
            multiplyArrayScalar<<<60,MAXTHREADS>>>(d_dimgTV,alpha,   total_pixels);
            //SUBSTRACT GRADIENT
            substractArrays    <<<60,MAXTHREADS>>>(d_image,d_dimgTV, total_pixels);
            sumnorm2=0;
        }
        
        cudaCheckErrors("TV minimization");
        
        cudaMemcpy(dst, d_image, mem_size, cudaMemcpyDeviceToHost);
        cudaCheckErrors("Copy result back");
        
        cudaFree(d_image);
        cudaFree(d_norm2aux);
        cudaFree(d_dimgTV);
        
        cudaCheckErrors("Memory free");
        
    }
    
