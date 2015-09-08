/**
 * @file mxGPUArray.h
 * @brief Exported interface to GPU data.
 *
 * This header declares an interface that allows MEX function access to data on
 * the GPU. It is designed to decouple the MEX author from the details of the
 * MathWorks Parallel Computing Toolbox (PCT) GPU implementation.
 *
 * \par Notes on creating and destroying mxGPUArray objects
 * Whenever you call any functions that return an mxGPUArray object, you should
 * always call mxGPUDestroyGPUArray() on the result when you are done.  This
 * advice differs from the advice for general MATLAB mxArrays. General mxArrays
 * should only be destroyed if they are created inside the present mex function
 * and are not being returned. For mxGPUArrays, a small wrapper object is
 * created on the CPU that refers to the internal data structure. That wrapper
 * object must be destroyed before you leave the mex function. This applies even
 * if you plan to return the array to MATLAB.
 *
 * \par
 * We guarantee that all asynchronous operations on the CUDA device have
 * finished before we give you a pointer to data on the device. We expect that
 * your MEX function will similarly complete any asynchronous operations on the
 * CUDA device before returning a gpuArray to MATLAB. Otherwise, undefined
 * behavior will result.
 *
 * \par
 * Do not attempt to wrap the same mxGPUArray object in two different
 * mxArrays. This will cause undefined behavior.
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#if defined(_MSC_VER)
# pragma once
#endif
#if defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ > 3))
# pragma once
#endif

#ifndef _MX_GPUARRAY_PUBLISHED_API_H
#define _MX_GPUARRAY_PUBLISHED_API_H

#ifndef LIBMWGPU_API
#  define LIBMWGPU_API
#endif

#ifndef EXTERN_C
#  ifdef __cplusplus
#    define EXTERN_C extern "C"
#  else
#    define EXTERN_C extern
#  endif
#endif

/* Incomplete typedef for mxGPUArray */
#ifndef _MX_GPUARRAY_DEFINED
#  ifdef __cplusplus
class mxGPUArray;
#  else
typedef struct mxGPUArray_tag mxGPUArray;
#  endif
/*lint -esym(1923,_MX_GPUARRAY_DEFINED) // MACRO input cannot be converted to const variable in C*/
#define _MX_GPUARRAY_DEFINED 1.0
#endif

#include "matrix.h"

enum
{
    MX_GPU_SUCCESS = 0,
    MX_GPU_FAILURE = 1
};

typedef enum
{
    MX_GPU_DO_NOT_INITIALIZE = 0,
    MX_GPU_INITIALIZE_VALUES = 1
} mxGPUInitialize;

/**
 * @brief Initialize the MathWorks GPU library on the currently selected device.
 *
 * Before you use any CUDA code in your MEX file, you must initialize the
 * MathWorks GPU library if you intend to use any mxGPU functionality in MEX or
 * any GPU calls in MATLAB. There are many ways to initialize the MathWorks GPU
 * API, but here are three possible ways:
 * 1) Call mxInitGPU at the beginning of your MEX file.
 * 2) Call gpuDevice(deviceIndex) in MATLAB.
 * 3) Create a gpuArray in MATLAB.
 *
 * We recommend that you call mxInitGPU at the beginning of your MEX file unless
 * you have an alternate way of guaranteeing that the MathWorks GPU API has
 * been initialized at the start of your MEX file.
 *
 * If the library is already initialized, this function returns without doing
 * any work. If the library has not been initialized, the function initializes
 * the default device.
 * @note At present, MATLAB can only work with one GPU device at a time.
 *
 * @returns MX_GPU_SUCCESS if the MathWorks GPU API has been successfully
 * initialized, MX_GPU_FAILURE otherwise.
 */
EXTERN_C LIBMWGPU_API int                mxInitGPU(void);

/**
 * @brief Return non-zero if the input mxArray contains GPU data.
 */
EXTERN_C LIBMWGPU_API int                mxIsGPUArray(mxArray const * const);

/**
 * @brief Create a new read-only mxGPUArray object from an input mxArray.
 * @param in  An input mxArray from MATLAB containing either GPU or CPU data.
 *
 * Produces a read-only mxGPUArray object from an mxArray passed as an input to
 * the function. If the input mxArray contains a gpuArray, this function
 * extracts a reference to the GPU data from an mxArray passed as an input to
 * the function. If the input mxArray contains CPU data, the data will be copied
 * to the GPU, but the returned object is still read-only. If you need a
 * writable copy, use mxGPUCopyFromMxArray instead.
 * @note This function allocates a new mxGPUArray object on the CPU. Call
 * mxGPUDestroyGPUArray on the result when you are done.
 */
EXTERN_C LIBMWGPU_API mxGPUArray const * mxGPUCreateFromMxArray(mxArray const * const);

/**
 * @brief Create a new mxGPUArray copy from an input mxArray.
 * @param in  An mxArray containing either GPU or CPU data.
 *
 * Produces a new mxGPUArray object with the same characteristics as the
 * mxArray. If the input mxArray contains a gpuArray, this function makes a new
 * copy of the data on the GPU. If the input mxArray contains numeric or logical
 * CPU data, the data will be copied to the GPU. Either way, this function
 * always allocates memory on the GPU.
 * @note This function allocates a new mxGPUArray object on the CPU. Call
 * mxGPUDestroyGPUArray on the result when you are done.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCopyFromMxArray(mxArray const * const);

/**
 * @brief Create a new mxGPUArray object, allocating memory on the GPU.
 *
 * Create a new mxGPUArray object with the given size, type, and complexity,
 * optionally filling the memory with zeros. Allocates memory on the GPU.
 * @note This function allocates a new mxGPUArray object on the CPU. Call
 * mxGPUDestroyGPUArray on the result when you are done.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCreateGPUArray(mwSize const ndims,
                                                             mwSize const *dims,
                                                             mxClassID const theClassID,
                                                             mxComplexity const theComplexity,
                                                             mxGPUInitialize const initializeToZero);

/**
 * @brief Return the mxClassID associated with the data on the GPU.
 */
EXTERN_C LIBMWGPU_API mxClassID          mxGPUGetClassID(mxGPUArray const *);

/**
 * @brief Return the complexity of the data on the GPU.
 */

EXTERN_C LIBMWGPU_API mxComplexity       mxGPUGetComplexity(mxGPUArray const *);

/**
 * @brief Returns a raw pointer to the underlying data.
 *
 * Returns a raw pointer to the underlying data. Cast it to the type of data
 * that you wish to use on the device.
 * @warning It is the caller's responsibility to check that the data inside the
 * array has the appropriate type.
 */
EXTERN_C LIBMWGPU_API void *             mxGPUGetData(mxGPUArray * const);

/**
 * @brief Returns a read-only raw pointer to the underlying data.
 * Returns a read-only raw pointer to the underlying data. Cast it to the type
 * of data that you wish to use on the device.
 * @warning It is the caller's responsibility to check that the data inside the
 * array has the appropriate type.
 */
EXTERN_C LIBMWGPU_API void const *       mxGPUGetDataReadOnly(mxGPUArray const * const);

/**
 * @brief Return an array of dimensions.
 * @returns A pointer to an array of mwSize that must be deleted using mxFree.
 */
EXTERN_C LIBMWGPU_API mwSize const *     mxGPUGetDimensions(mxGPUArray const * const);

/**
 * @brief Return the size of the dimension array for this array.
 */
EXTERN_C LIBMWGPU_API mwSize             mxGPUGetNumberOfDimensions(mxGPUArray const * const);

/**
 * @brief Return the total number of elements on the GPU for this array.
 */
EXTERN_C LIBMWGPU_API mwSize             mxGPUGetNumberOfElements(mxGPUArray const * const);

/**
 * @brief Create an mxArray for returning GPU data to MATLAB.
 * @returns A new mxArray containing the GPU data.
 *
 * Put the mxGPUArray into an mxArray for return to MATLAB. The data remains on
 * the GPU and the returned class in MATLAB will be 'gpuArray'. After this call,
 * the mxGPUArray object is no longer needed and can be destroyed.
 */
EXTERN_C LIBMWGPU_API mxArray *          mxGPUCreateMxArrayOnGPU(mxGPUArray const * const);

/**
 * @brief Create an mxArray for returning CPU data to MATLAB, copying data from the GPU.
 * @returns A new mxArray containing CPU data that is a copy of the GPU data.
 *
 * Copy the GPU data from the mxGPUArray into an mxArray on the CPU for return
 * to MATLAB. After this call, the mxGPUArray object is no longer needed and can
 * be destroyed.
 */
EXTERN_C LIBMWGPU_API mxArray *          mxGPUCreateMxArrayOnCPU(mxGPUArray const * const);

/**
 * @brief Delete an mxGPUArray object.
 *
 * Deletes an mxGPUArray object on the CPU.  Call if you created a new
 * mxGPUArray with mxGPUCreateGPUArray, mxGPUCreateFromMxArray, mxGPUCopyFromMxArray,
 * mxGPUCopyReal, mxGPUCopyImag, or mxGPUCreateComplexGPUArray.
 *
 * This call will deallocate memory on the GPU, unless some other mxArray holds
 * a reference to the same data. For example, if the mxGPUArray was extracted
 * from an input mxArray, or wrapped in an mxArray for an output, then the data
 * on the GPU will remain.
 */
EXTERN_C LIBMWGPU_API void               mxGPUDestroyGPUArray(mxGPUArray const * const);

/**
 * @brief Duplicate (deep copy) an mxGPUArray object.
 *
 * Produces a new array on the GPU and copies the data. Returns a new mxGPUArray
 * that refers to the copy. Call mxGPUDestroyGPUArray on the object when you are
 * done with it.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCopyGPUArray(mxGPUArray const * const);

/**
 * @brief Copy the real part of an mxGPUArray.
 *
 * Copies the real part of GPU data. Returns a new mxGPUArray object that refers
 * to the copy. Call mxGPUDestroyGPUArray on the object when you are done with it.
 *
 * Allocates memory on the GPU and copies data. If the input is real rather than
 * complex it will return a copy of the input.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCopyReal(mxGPUArray const * const);

/**
 * @brief Copy the imaginary part of an mxGPUArray.
 *
 * Copies the imaginary part of GPU data. Returns a new mxGPUArray object that
 * refers to the copy. Call mxGPUDestroyGPUArray on the object when you are done
 * with it.
 *
 * Allocates memory on the GPU and copies data. If the input is real rather than
 * complex it will return an array of zeros.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCopyImag(mxGPUArray const * const);

/**
 * @brief Create a complex GPU array by copying two real GPU arrays.
 *
 * Creates a new complex mxGPUArray from two real mxGPUArrays. Call
 * mxGPUDestroyGPUArray on the object when you are done with it.
 *
 * Allocates memory on the GPU and copies data. Requires that each input is real
 * and that the input sizes and classes match.
 */
EXTERN_C LIBMWGPU_API mxGPUArray *       mxGPUCreateComplexGPUArray(mxGPUArray const * const realPart,
                                                                    mxGPUArray const * const imagPart);

/**
 * @brief Returns non-zero if two mxGPUArrays refer to the same GPU data.
 */
EXTERN_C LIBMWGPU_API int                mxGPUIsSame(mxGPUArray const * const in1,
                                                     mxGPUArray const * const in2);

/**
 * @brief Returns non-zero if the mxArray is a pointer to valid GPU data.
 *
 * If the device is reinitialized in MATLAB using gpuDevice(), all GPU data
 * becomes invalid, but the CPU structures that describe the GPU data still
 * exist. This function checks whether the mxArray is a container of valid GPU
 * data. It returns false for CPU data and false for invalid GPU data.
 */
EXTERN_C LIBMWGPU_API int                mxGPUIsValidGPUData(mxArray const * const);

#endif
