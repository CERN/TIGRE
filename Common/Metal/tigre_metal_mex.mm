/*
 * Apple Metal backend for MATLAB TIGRE.
 *
 * All numerical operations execute as Metal compute kernels. Buffers use
 * shared storage on Apple Silicon and the device, queue, library, and compute
 * pipelines remain cached for the lifetime of the loaded MEX file.
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "mex.h"
#include "matrix.h"
#include "../Geometry/tigre_geometry.hpp"
#include "tigre_metal_limits.hpp"

#include <algorithm>
#ifdef TIGRE_METAL_TIMING
#include <chrono>
#endif
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using Geometry = tigre_geometry::Geometry;
using ViewGeometry = tigre_geometry::View;

struct Params {
    std::uint32_t nx;
    std::uint32_t ny;
    std::uint32_t nz;
    std::uint32_t nu;
    std::uint32_t nv;
    std::uint32_t angleCount;
    std::uint32_t mode;
    std::uint32_t variant;
    float sx;
    float sy;
    float sz;
    float dx;
    float dy;
    float dz;
    float du;
    float dv;
    float su;
    float sv;
    float accuracy;
    float alpha;
    float ratio;
    float lambda;
    float tau1;
    float tau2;
};

static_assert(sizeof(Params) == 96, "Metal parameter layout mismatch.");

id<MTLDevice> gDevice = nil;
id<MTLCommandQueue> gQueue = nil;
id<MTLLibrary> gLibrary = nil;
NSMutableDictionary<NSString*, id<MTLComputePipelineState>>* gPipelines = nil;
bool gExitRegistered = false;

#ifdef TIGRE_METAL_TIMING
using HostClock = std::chrono::steady_clock;

double milliseconds(HostClock::time_point start,
                    HostClock::time_point finish) {
    return std::chrono::duration<double, std::milli>(
        finish - start).count();
}

struct AxTiming {
    double geometryMilliseconds = 0.0;
    double matlabOutputAllocationMilliseconds = 0.0;
    double inputSetupMilliseconds = 0.0;
    double metalOutputAllocationMilliseconds = 0.0;
    double gpuHostMilliseconds = 0.0;
    double gpuMilliseconds = 0.0;
    double outputCopyMilliseconds = 0.0;
    double totalMilliseconds = 0.0;
    NSUInteger threadExecutionWidth = 0;
    NSUInteger maximumThreadsPerThreadgroup = 0;
};
#endif

void cleanupRuntime() {
    gPipelines = nil;
    gLibrary = nil;
    gQueue = nil;
    gDevice = nil;
    gExitRegistered = false;
}

class MetalError : public std::runtime_error {
public:
    const char* id;
    MetalError(const char* errorId, const std::string& message)
        : std::runtime_error(message), id(errorId) {}
};

[[noreturn]] void fail(const char* id, const std::string& message) {
    // Unwind C++/ARC resources before reporting the error at the MEX boundary.
    throw MetalError(id, message);
}

std::string nsString(NSString* value) {
    return value == nil ? std::string() : std::string([value UTF8String]);
}

std::string getString(const mxArray* value, const char* id, const char* name) {
    if (!mxIsChar(value)) {
        fail(id, std::string(name) + " must be a character vector.");
    }
    char* raw = mxArrayToString(value);
    if (raw == nullptr) {
        fail(id, std::string("Unable to read ") + name + ".");
    }
    std::string result(raw);
    mxFree(raw);
    return result;
}

NSString* locateLibrary() {
    mxArray* input = mxCreateString("tigre_metal_mex");
    mxArray* output = nullptr;
    const int status = mexCallMATLAB(1, &output, 1, &input, "which");
    mxDestroyArray(input);
    if (status != 0 || output == nullptr) {
        if (output != nullptr) {
            mxDestroyArray(output);
        }
        fail("TIGRE:Metal:LibraryNotFound",
             "Unable to locate the loaded Metal backend binary.");
    }
    char* path = mxArrayToString(output);
    mxDestroyArray(output);
    if (path == nullptr) {
        fail("TIGRE:Metal:LibraryNotFound",
             "Unable to resolve the Metal backend folder.");
    }
    NSString* binaryPath = [NSString stringWithUTF8String:path];
    mxFree(path);
    return [[binaryPath stringByDeletingLastPathComponent]
        stringByAppendingPathComponent:@"tigre_kernels.metallib"];
}

void initializeRuntime() {
    if (gDevice != nil) {
        return;
    }
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (device == nil) {
        fail("TIGRE:Metal:DeviceUnavailable",
             "No Metal-capable Apple GPU is available.");
    }
    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (queue == nil) {
        fail("TIGRE:Metal:InitializationFailed",
             "Unable to create the Metal command queue.");
    }
    NSString* libraryPath = locateLibrary();
    NSError* error = nil;
    id<MTLLibrary> library = [device newLibraryWithURL:
        [NSURL fileURLWithPath:libraryPath] error:&error];
    if (library == nil) {
        fail("TIGRE:Metal:LibraryNotFound",
             "Unable to load " + nsString(libraryPath) + ": " +
             nsString(error.localizedDescription));
    }
    NSMutableDictionary<NSString*, id<MTLComputePipelineState>>* pipelines =
        [[NSMutableDictionary alloc] init];
    if (pipelines == nil) {
        fail("TIGRE:Metal:InitializationFailed",
             "Unable to create the Metal pipeline cache.");
    }
    if (!gExitRegistered) {
        mexAtExit(cleanupRuntime);
        gExitRegistered = true;
    }
    // A non-nil gDevice means the entire runtime is ready, never partially built.
    gQueue = queue;
    gLibrary = library;
    gPipelines = pipelines;
    gDevice = device;
}

std::size_t bufferBytes(std::size_t count, std::size_t elementSize) {
    initializeRuntime();
    return tigre_metal_limits::bytes(count, elementSize, gDevice.maxBufferLength);
}

id<MTLComputePipelineState> pipeline(NSString* name) {
    initializeRuntime();
    id<MTLComputePipelineState> result = gPipelines[name];
    if (result != nil) {
        return result;
    }
    id<MTLFunction> function = [gLibrary newFunctionWithName:name];
    if (function == nil) {
        fail("TIGRE:Metal:KernelUnavailable",
             "Metal kernel not found: " + nsString(name));
    }
    NSError* error = nil;
    result = [gDevice newComputePipelineStateWithFunction:function error:&error];
    if (result == nil) {
        fail("TIGRE:Metal:PipelineFailed",
             "Unable to create Metal pipeline " + nsString(name) + ": " +
             nsString(error.localizedDescription));
    }
    gPipelines[name] = result;
    return result;
}

id<MTLBuffer> bufferWithBytes(const void* bytes, std::size_t length) {
    bufferBytes(length, 1);
    id<MTLBuffer> result = [gDevice newBufferWithBytes:bytes
        length:length options:MTLResourceStorageModeShared];
    if (result == nil) {
        fail("TIGRE:Metal:AllocationFailed",
             "The Metal backend could not allocate its whole-array buffers.");
    }
    return result;
}

id<MTLBuffer> emptyBuffer(std::size_t length, bool clear = false) {
    bufferBytes(length, 1);
    id<MTLBuffer> result = [gDevice newBufferWithLength:length
        options:MTLResourceStorageModeShared];
    if (result == nil) {
        fail("TIGRE:Metal:AllocationFailed",
             "The Metal backend could not allocate its whole-array buffers.");
    }
    if (clear) {
        std::memset([result contents], 0, length);
    }
    return result;
}

id<MTLBuffer> noCopyBuffer(void* bytes, std::size_t length) {
    bufferBytes(length, 1);
    const std::size_t pageSize = static_cast<std::size_t>(getpagesize());
    const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(bytes);
    if (bytes == nullptr || length == 0 || address % pageSize != 0 ||
        length % pageSize != 0) {
        return nil;
    }
    return [gDevice newBufferWithBytesNoCopy:bytes
        length:length
        options:MTLResourceStorageModeShared
        deallocator:nil];
}

NSUInteger threadWidth(id<MTLComputePipelineState> state) {
    NSUInteger width = std::min<NSUInteger>(
        256, state.maxTotalThreadsPerThreadgroup);
    NSUInteger power = 1;
    while (power * 2 <= width) {
        power *= 2;
    }
    return power;
}

void dispatch(id<MTLComputeCommandEncoder> encoder, NSString* kernelName,
              std::size_t count) {
    tigre_metal_limits::checkedProduct(
        count, 1, std::numeric_limits<std::uint32_t>::max());
    id<MTLComputePipelineState> state = pipeline(kernelName);
    [encoder setComputePipelineState:state];
    const NSUInteger width = threadWidth(state);
    [encoder dispatchThreads:MTLSizeMake(count, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
}

void complete(id<MTLCommandBuffer> commandBuffer) {
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    if (commandBuffer.status == MTLCommandBufferStatusError) {
        fail("TIGRE:Metal:ExecutionFailed",
             "Metal command failed: " +
             nsString(commandBuffer.error.localizedDescription));
    }
}

const float* singleInput(const mxArray* value, const char* description) {
    if (!mxIsSingle(value) || mxIsComplex(value)) {
        fail("TIGRE:Metal:InvalidInput",
             std::string(description) + " must be a real single array.");
    }
    return static_cast<const float*>(mxGetData(value));
}

struct ImageShape {
    std::uint32_t nx;
    std::uint32_t ny;
    std::uint32_t nz;
    std::uint32_t count;
    mwSize dimensions[3];
    mwSize dimensionCount;
};

ImageShape imageShape(const mxArray* image) {
    singleInput(image, "Image");
    const mwSize dimensionCount = mxGetNumberOfDimensions(image);
    if (dimensionCount < 2 || dimensionCount > 3) {
        fail("TIGRE:Metal:TV:InvalidInput", "Image must be 2-D or 3-D.");
    }
    const mwSize* dimensions = mxGetDimensions(image);
    const std::size_t nz = dimensionCount == 3 ? dimensions[2] : 1;
    const std::size_t count =
        tigre_metal_limits::elements({dimensions[0], dimensions[1], nz});
    bufferBytes(count, sizeof(float));
    ImageShape result{};
    result.nx = static_cast<std::uint32_t>(dimensions[0]);
    result.ny = static_cast<std::uint32_t>(dimensions[1]);
    result.nz = static_cast<std::uint32_t>(nz);
    result.count = static_cast<std::uint32_t>(count);
    result.dimensionCount = dimensionCount;
    result.dimensions[0] = dimensions[0];
    result.dimensions[1] = dimensions[1];
    result.dimensions[2] = nz;
    return result;
}

Params shapeParams(const ImageShape& shape) {
    Params parameters{};
    parameters.nx = shape.nx;
    parameters.ny = shape.ny;
    parameters.nz = shape.nz;
    return parameters;
}

Params projectionParams(const Geometry& geometry) {
    const std::size_t voxels = tigre_metal_limits::elements(
        {geometry.nx, geometry.ny, geometry.nz});
    const std::size_t projections = tigre_metal_limits::elements(
        {geometry.nu, geometry.nv, geometry.views.size()});
    // Validate all whole-array allocations before creating the MATLAB output.
    bufferBytes(voxels, sizeof(float));
    bufferBytes(projections, sizeof(float));
    bufferBytes(geometry.views.size(), sizeof(ViewGeometry));
    Params parameters{};
    parameters.nx = static_cast<std::uint32_t>(geometry.nx);
    parameters.ny = static_cast<std::uint32_t>(geometry.ny);
    parameters.nz = static_cast<std::uint32_t>(geometry.nz);
    parameters.nu = static_cast<std::uint32_t>(geometry.nu);
    parameters.nv = static_cast<std::uint32_t>(geometry.nv);
    parameters.angleCount =
        static_cast<std::uint32_t>(geometry.views.size());
    parameters.mode = geometry.parallel ? 1u : 0u;
    parameters.sx = geometry.sx;
    parameters.sy = geometry.sy;
    parameters.sz = geometry.sz;
    parameters.dx = geometry.dx;
    parameters.dy = geometry.dy;
    parameters.dz = geometry.dz;
    parameters.du = geometry.du;
    parameters.dv = geometry.dv;
    parameters.su = geometry.su;
    parameters.sv = geometry.sv;
    parameters.accuracy = geometry.accuracy;
    return parameters;
}

mxArray* createLike(const ImageShape& shape) {
    return mxCreateNumericArray(
        shape.dimensionCount, shape.dimensions, mxSINGLE_CLASS, mxREAL);
}

float scalarFloat(const mxArray* value, const char* name) {
    if (!mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1) {
        fail("TIGRE:Metal:TV:InvalidInput",
             std::string(name) + " must be a real scalar.");
    }
    const double scalar = mxGetScalar(value);
    if (!std::isfinite(scalar)) {
        fail("TIGRE:Metal:TV:InvalidInput",
             std::string(name) + " must be finite.");
    }
    return static_cast<float>(scalar);
}

int scalarInteger(const mxArray* value, const char* name) {
    const float scalar = scalarFloat(value, name);
    if (scalar < 0.0f) {
        fail("TIGRE:Metal:TV:InvalidInput",
             std::string(name) + " must be nonnegative.");
    }
    return static_cast<int>(std::floor(scalar + 0.5f));
}

void copyToOutput(mxArray* output, id<MTLBuffer> buffer,
                  std::size_t bytes) {
    std::memcpy(mxGetData(output), [buffer contents], bytes);
}

void runProjectionKernel(NSString* kernelName, id<MTLBuffer> input,
                         id<MTLBuffer> views, id<MTLBuffer> output,
                         const Params& parameters, std::size_t count
#ifdef TIGRE_METAL_TIMING
                         , AxTiming* timing
#endif
                         ) {
    id<MTLCommandBuffer> command = [gQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setBuffer:input offset:0 atIndex:0];
    [encoder setBuffer:views offset:0 atIndex:1];
    [encoder setBuffer:output offset:0 atIndex:2];
    [encoder setBytes:&parameters length:sizeof(parameters) atIndex:3];
#ifdef TIGRE_METAL_TIMING
    id<MTLComputePipelineState> state = pipeline(kernelName);
    if (timing != nullptr) {
        timing->threadExecutionWidth = state.threadExecutionWidth;
        timing->maximumThreadsPerThreadgroup =
            state.maxTotalThreadsPerThreadgroup;
    }
    const HostClock::time_point gpuStart = HostClock::now();
#endif
    dispatch(encoder, kernelName, count);
    [encoder endEncoding];
    complete(command);
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point gpuFinish = HostClock::now();
    if (timing != nullptr) {
        timing->gpuHostMilliseconds = milliseconds(gpuStart, gpuFinish);
        const CFTimeInterval start = command.GPUStartTime;
        const CFTimeInterval finish = command.GPUEndTime;
        if (finish >= start && start > 0.0) {
            timing->gpuMilliseconds = (finish - start) * 1000.0;
        }
    }
#endif
}

void runAx(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5 || nlhs != 1) {
        fail("TIGRE:Metal:Ax:InvalidInput",
             "Ax requires image, geometry, angles, and projector type.");
    }
#ifdef TIGRE_METAL_TIMING
    AxTiming timing;
    const HostClock::time_point totalStart = HostClock::now();
    const HostClock::time_point geometryStart = HostClock::now();
#endif
    const float* image = singleInput(prhs[1], "Image");
    const std::string projector =
        getString(prhs[4], "TIGRE:Metal:Ax:InvalidInput", "Projector type");
    const bool siddon =
        projector == "Siddon" || projector == "ray-voxel";
    const Geometry geometry = tigre_geometry::parse(
        prhs[2], prhs[3], "TIGRE:Metal:BadGeometry", !siddon);
    const std::size_t angleCount = geometry.views.size();
    Params parameters = projectionParams(geometry);
    if (siddon) {
        parameters.variant = 1;
    } else if (projector == "interpolated") {
        parameters.variant = 0;
    } else {
        fail("TIGRE:Metal:Ax:InvalidInput", "Unknown projector type.");
    }
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point geometryFinish = HostClock::now();
    timing.geometryMilliseconds =
        milliseconds(geometryStart, geometryFinish);
#endif
    const std::size_t imageCount =
        std::size_t(parameters.nx) * parameters.ny * parameters.nz;
    if (mxGetNumberOfElements(prhs[1]) != imageCount) {
        fail("TIGRE:Metal:Ax:InvalidInput",
             "Image dimensions do not match geo.nVoxel.");
    }
    const std::size_t outputCount =
        std::size_t(parameters.nu) * parameters.nv * angleCount;
    const mwSize dimensions[3] = {
        parameters.nv, parameters.nu, static_cast<mwSize>(angleCount)};
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point matlabAllocationStart = HostClock::now();
#endif
    plhs[0] = mxCreateNumericArray(3, dimensions, mxSINGLE_CLASS, mxREAL);
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point matlabAllocationFinish = HostClock::now();
    timing.matlabOutputAllocationMilliseconds = milliseconds(
        matlabAllocationStart, matlabAllocationFinish);
    const HostClock::time_point inputSetupStart = HostClock::now();
#endif
    const std::size_t imageBytes = imageCount * sizeof(float);
    id<MTLBuffer> imageBuffer =
        noCopyBuffer(const_cast<float*>(image), imageBytes);
    const bool imageNoCopy = imageBuffer != nil;
    if (!imageNoCopy) {
        imageBuffer = bufferWithBytes(image, imageBytes);
    }
    id<MTLBuffer> viewBuffer = bufferWithBytes(
        geometry.views.data(), angleCount * sizeof(ViewGeometry));
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point inputSetupFinish = HostClock::now();
    timing.inputSetupMilliseconds =
        milliseconds(inputSetupStart, inputSetupFinish);
    const HostClock::time_point metalAllocationStart = HostClock::now();
#endif
    const std::size_t outputBytes = outputCount * sizeof(float);
    id<MTLBuffer> outputBuffer =
        noCopyBuffer(mxGetData(plhs[0]), outputBytes);
    const bool outputNoCopy = outputBuffer != nil;
    if (!outputNoCopy) {
        outputBuffer = emptyBuffer(outputBytes);
    }
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point metalAllocationFinish = HostClock::now();
    timing.metalOutputAllocationMilliseconds = milliseconds(
        metalAllocationStart, metalAllocationFinish);
#endif
    runProjectionKernel(@"tigre_ax", imageBuffer, viewBuffer, outputBuffer,
                        parameters, outputCount
#ifdef TIGRE_METAL_TIMING
                        , &timing
#endif
                        );
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point outputCopyStart = HostClock::now();
#endif
    if (!outputNoCopy) {
        copyToOutput(plhs[0], outputBuffer, outputBytes);
    }
#ifdef TIGRE_METAL_TIMING
    const HostClock::time_point outputCopyFinish = HostClock::now();
    timing.outputCopyMilliseconds =
        milliseconds(outputCopyStart, outputCopyFinish);
    timing.totalMilliseconds = milliseconds(totalStart, outputCopyFinish);
    mexPrintf(
        "[TIGRE Metal Ax timing] projector=%s rays=%zu "
        "geometry=%.3fms matlabOutputAlloc=%.3fms "
        "inputSharedSetup=%.3fms metalOutputAlloc=%.3fms "
        "gpuHost=%.3fms gpu=%.3fms outputCopy=%.3fms total=%.3fms "
        "inputNoCopy=%d outputNoCopy=%d "
        "threadExecutionWidth=%lu maxThreads=%lu\n",
        projector.c_str(), outputCount, timing.geometryMilliseconds,
        timing.matlabOutputAllocationMilliseconds,
        timing.inputSetupMilliseconds,
        timing.metalOutputAllocationMilliseconds,
        timing.gpuHostMilliseconds, timing.gpuMilliseconds,
        timing.outputCopyMilliseconds, timing.totalMilliseconds,
        imageNoCopy ? 1 : 0, outputNoCopy ? 1 : 0,
        static_cast<unsigned long>(timing.threadExecutionWidth),
        static_cast<unsigned long>(timing.maximumThreadsPerThreadgroup));
#endif
}

void runAtb(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5 || nlhs != 1) {
        fail("TIGRE:Metal:Atb:InvalidInput",
             "Atb requires projections, geometry, angles, and backprojector type.");
    }
    const float* projections = singleInput(prhs[1], "Projections");
    const Geometry geometry =
        tigre_geometry::parse(prhs[2], prhs[3], "TIGRE:Metal:BadGeometry");
    const std::size_t angleCount = geometry.views.size();
    Params parameters = projectionParams(geometry);
    const std::string projector =
        getString(prhs[4], "TIGRE:Metal:Atb:InvalidInput", "Backprojector type");
    if (projector == "matched") {
        parameters.variant = 1;
    } else if (projector == "FDK") {
        parameters.variant = 0;
    } else {
        fail("TIGRE:Metal:Atb:InvalidInput", "Unknown backprojector type.");
    }
    const std::size_t projectionCount =
        std::size_t(parameters.nu) * parameters.nv * angleCount;
    if (mxGetNumberOfElements(prhs[1]) != projectionCount) {
        fail("TIGRE:Metal:Atb:InvalidInput",
             "Projection dimensions do not match the geometry and angles.");
    }
    const std::size_t outputCount =
        std::size_t(parameters.nx) * parameters.ny * parameters.nz;
    const mwSize dimensions[3] = {
        parameters.nx, parameters.ny, parameters.nz};
    plhs[0] = mxCreateNumericArray(3, dimensions, mxSINGLE_CLASS, mxREAL);
    const std::size_t projectionBytes =
        projectionCount * sizeof(float);
    id<MTLBuffer> projectionBuffer =
        noCopyBuffer(const_cast<float*>(projections), projectionBytes);
    if (projectionBuffer == nil) {
        projectionBuffer = bufferWithBytes(projections, projectionBytes);
    }
    id<MTLBuffer> viewBuffer = bufferWithBytes(
        geometry.views.data(), angleCount * sizeof(ViewGeometry));
    const std::size_t outputBytes = outputCount * sizeof(float);
    id<MTLBuffer> outputBuffer =
        noCopyBuffer(mxGetData(plhs[0]), outputBytes);
    const bool outputNoCopy = outputBuffer != nil;
    if (!outputNoCopy) {
        outputBuffer = emptyBuffer(outputBytes);
    }
    runProjectionKernel(@"tigre_atb", projectionBuffer, viewBuffer,
                        outputBuffer, parameters, outputCount
#ifdef TIGRE_METAL_TIMING
                        , nullptr
#endif
                        );
    if (!outputNoCopy) {
        copyToOutput(plhs[0], outputBuffer, outputBytes);
    }
}

void runTvGradient(id<MTLBuffer> image, id<MTLBuffer> gradient,
                   const Params& parameters, std::size_t count) {
    id<MTLCommandBuffer> command = [gQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setBuffer:image offset:0 atIndex:0];
    [encoder setBuffer:gradient offset:0 atIndex:1];
    [encoder setBytes:&parameters length:sizeof(parameters) atIndex:3];
    dispatch(encoder, @"tigre_tv_gradient", count);
    [encoder endEncoding];
    complete(command);
}

float gpuNormSquared(id<MTLBuffer> values, std::uint32_t count) {
    id<MTLBuffer> current = values;
    std::uint32_t currentCount = count;
    bool squares = true;
    while (true) {
        id<MTLComputePipelineState> state =
            pipeline(squares ? @"tigre_reduce_squares" : @"tigre_reduce_sum");
        const NSUInteger width = threadWidth(state);
        const std::uint32_t groupCount =
            static_cast<std::uint32_t>((currentCount + width - 1) / width);
        id<MTLBuffer> output =
            emptyBuffer(std::size_t(groupCount) * sizeof(float));
        id<MTLCommandBuffer> command = [gQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder =
            [command computeCommandEncoder];
        [encoder setComputePipelineState:state];
        [encoder setBuffer:current offset:0 atIndex:0];
        [encoder setBuffer:output offset:0 atIndex:1];
        [encoder setBytes:&currentCount length:sizeof(currentCount) atIndex:2];
        [encoder dispatchThreadgroups:MTLSizeMake(groupCount, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
        [encoder endEncoding];
        complete(command);
        if (groupCount == 1) {
            return *static_cast<const float*>([output contents]);
        }
        current = output;
        currentCount = groupCount;
        squares = false;
    }
}

void scaleSubtract(id<MTLBuffer> image, id<MTLBuffer> gradient,
                   float scale, std::uint32_t count) {
    id<MTLCommandBuffer> command = [gQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setBuffer:image offset:0 atIndex:0];
    [encoder setBuffer:gradient offset:0 atIndex:1];
    [encoder setBytes:&scale length:sizeof(scale) atIndex:2];
    [encoder setBytes:&count length:sizeof(count) atIndex:3];
    dispatch(encoder, @"tigre_scale_subtract", count);
    [encoder endEncoding];
    complete(command);
}

void gradientDescentTv(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[], bool adaptive) {
    const int expectedInputs = adaptive ? 5 : 4;
    if (nrhs != expectedInputs || nlhs != 1) {
        fail("TIGRE:Metal:TV:InvalidInput", "Invalid TV operation inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const ImageShape shape = imageShape(prhs[1]);
    const float stepSize = scalarFloat(prhs[2], "Step size");
    const int iterations = scalarInteger(prhs[3], "Iteration count");
    Params parameters = shapeParams(shape);
    parameters.variant = adaptive ? 1u : 0u;
    parameters.alpha =
        adaptive ? scalarFloat(prhs[4], "Delta") : 1.0f;
    const std::size_t bytes = std::size_t(shape.count) * sizeof(float);
    id<MTLBuffer> image = bufferWithBytes(input, bytes);
    id<MTLBuffer> gradient = emptyBuffer(bytes);
    for (int iteration = 0; iteration < iterations; ++iteration) {
        runTvGradient(image, gradient, parameters, shape.count);
        const float norm =
            std::sqrt(gpuNormSquared(gradient, shape.count) + 1.0e-7f);
        scaleSubtract(image, gradient, stepSize / norm, shape.count);
    }
    plhs[0] = createLike(shape);
    copyToOutput(plhs[0], image, bytes);
}

void runPiccs(int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (nrhs != 6 || nlhs != 1) {
        fail("TIGRE:Metal:TV:InvalidInput", "Invalid minPICCS inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const float* priorInput = singleInput(prhs[2], "Prior image");
    const ImageShape shape = imageShape(prhs[1]);
    const ImageShape priorShape = imageShape(prhs[2]);
    if (shape.dimensionCount != 3 || priorShape.dimensionCount != 3 ||
        shape.nx != priorShape.nx || shape.ny != priorShape.ny ||
        shape.nz != priorShape.nz) {
        fail("TIGRE:Metal:TV:InvalidInput",
             "Image and prior must be equally sized 3-D arrays.");
    }
    const float stepSize = scalarFloat(prhs[3], "Step size");
    const int iterations = scalarInteger(prhs[4], "Iteration count");
    const float ratio = scalarFloat(prhs[5], "PICCS ratio");
    Params parameters = shapeParams(shape);
    parameters.variant = 0;
    const std::size_t bytes = std::size_t(shape.count) * sizeof(float);
    id<MTLBuffer> image = bufferWithBytes(input, bytes);
    id<MTLBuffer> prior = bufferWithBytes(priorInput, bytes);
    id<MTLBuffer> difference = emptyBuffer(bytes);
    id<MTLBuffer> gradient = emptyBuffer(bytes);
    id<MTLBuffer> priorGradient = emptyBuffer(bytes);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        runTvGradient(image, gradient, parameters, shape.count);
        id<MTLCommandBuffer> differenceCommand = [gQueue commandBuffer];
        id<MTLComputeCommandEncoder> differenceEncoder =
            [differenceCommand computeCommandEncoder];
        [differenceEncoder setBuffer:image offset:0 atIndex:0];
        [differenceEncoder setBuffer:prior offset:0 atIndex:1];
        [differenceEncoder setBuffer:difference offset:0 atIndex:2];
        [differenceEncoder setBytes:&shape.count
            length:sizeof(shape.count) atIndex:3];
        dispatch(differenceEncoder, @"tigre_difference", shape.count);
        [differenceEncoder endEncoding];
        complete(differenceCommand);

        runTvGradient(difference, priorGradient, parameters, shape.count);
        id<MTLCommandBuffer> combineCommand = [gQueue commandBuffer];
        id<MTLComputeCommandEncoder> combineEncoder =
            [combineCommand computeCommandEncoder];
        [combineEncoder setBuffer:gradient offset:0 atIndex:0];
        [combineEncoder setBuffer:priorGradient offset:0 atIndex:1];
        [combineEncoder setBytes:&ratio length:sizeof(ratio) atIndex:2];
        [combineEncoder setBytes:&shape.count
            length:sizeof(shape.count) atIndex:3];
        dispatch(combineEncoder, @"tigre_combine", shape.count);
        [combineEncoder endEncoding];
        complete(combineCommand);

        const float normSquared = gpuNormSquared(gradient, shape.count);
        if (normSquared > 0.0f) {
            scaleSubtract(image, gradient,
                          stepSize / std::sqrt(normSquared), shape.count);
        }
    }
    plhs[0] = createLike(shape);
    copyToOutput(plhs[0], image, bytes);
}

void runTvDenoise(int nlhs, mxArray* plhs[], int nrhs,
                  const mxArray* prhs[]) {
    if (nrhs != 4 || nlhs != 1) {
        fail("TIGRE:Metal:TV:InvalidInput", "Invalid tvDenoise inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const ImageShape shape = imageShape(prhs[1]);
    if (shape.dimensionCount != 3) {
        fail("TIGRE:Metal:TV:InvalidInput",
             "tvDenoise requires a 3-D image.");
    }
    const float lambda = scalarFloat(prhs[2], "Lambda");
    const int iterations = scalarInteger(prhs[3], "Iteration count");
    if (std::abs(lambda) <= 1.0e-8f) {
        fail("TIGRE:Metal:TV:InvalidInput", "Lambda must be nonzero.");
    }
    Params parameters = shapeParams(shape);
    parameters.lambda = lambda;
    const std::size_t bytes = std::size_t(shape.count) * sizeof(float);
    id<MTLBuffer> source = bufferWithBytes(input, bytes);
    id<MTLBuffer> image = bufferWithBytes(input, bytes);
    id<MTLBuffer> pz = emptyBuffer(bytes, true);
    id<MTLBuffer> py = emptyBuffer(bytes, true);
    id<MTLBuffer> px = emptyBuffer(bytes, true);
    for (int iteration = 0; iteration < iterations; ++iteration) {
        parameters.tau2 = 0.3f + 0.02f * static_cast<float>(iteration);
        parameters.tau1 =
            (1.0f / parameters.tau2) *
            (1.0f / 6.0f -
             5.0f / (15.0f + static_cast<float>(iteration)));
        id<MTLCommandBuffer> command = [gQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder =
            [command computeCommandEncoder];
        [encoder setBuffer:source offset:0 atIndex:0];
        [encoder setBuffer:pz offset:0 atIndex:1];
        [encoder setBuffer:py offset:0 atIndex:2];
        [encoder setBuffer:px offset:0 atIndex:3];
        [encoder setBuffer:image offset:0 atIndex:4];
        [encoder setBytes:&parameters length:sizeof(parameters) atIndex:5];
        dispatch(encoder, @"tigre_tv_update_u", shape.count);

        [encoder setBuffer:image offset:0 atIndex:0];
        [encoder setBuffer:pz offset:0 atIndex:1];
        [encoder setBuffer:py offset:0 atIndex:2];
        [encoder setBuffer:px offset:0 atIndex:3];
        [encoder setBytes:&parameters length:sizeof(parameters) atIndex:4];
        dispatch(encoder, @"tigre_tv_update_p", shape.count);
        [encoder endEncoding];
        complete(command);
    }
    plhs[0] = createLike(shape);
    copyToOutput(plhs[0], image, bytes);
}

}  // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]) {
    // mexErrMsgIdAndTxt need not unwind C++ locals. Keep the final diagnostic
    // in stack storage; all owning C++/ARC resources leave scope first.
    const char* errorId = nullptr;
    char errorMessage[4096] = {};
    const auto recordError = [&](const char* id, const char* message) {
        errorId = id;
        std::snprintf(errorMessage, sizeof(errorMessage), "%s", message);
    };
    @autoreleasepool {
        try {
            if (nrhs < 1) {
                fail("TIGRE:Metal:InvalidInput",
                     "An operation name is required.");
            }
            const std::string operation =
                getString(prhs[0], "TIGRE:Metal:InvalidInput", "Operation");
            if (operation == "deviceName") {
                if (nrhs != 1 || nlhs != 1) {
                    fail("TIGRE:Metal:InvalidInput",
                         "deviceName requires one output and no extra inputs.");
                }
                initializeRuntime();
                plhs[0] = mxCreateString([gDevice.name UTF8String]);
            } else if (operation == "Ax" || operation == "ax") {
                initializeRuntime();
                runAx(nlhs, plhs, nrhs, prhs);
            } else if (operation == "Atb" || operation == "atb") {
                initializeRuntime();
                runAtb(nlhs, plhs, nrhs, prhs);
            } else if (operation == "minTV" || operation == "mintv") {
                initializeRuntime();
                gradientDescentTv(nlhs, plhs, nrhs, prhs, false);
            } else if (operation == "AwminTV" || operation == "awmintv") {
                initializeRuntime();
                gradientDescentTv(nlhs, plhs, nrhs, prhs, true);
            } else if (operation == "tvDenoise" ||
                       operation == "tvdenoise") {
                initializeRuntime();
                runTvDenoise(nlhs, plhs, nrhs, prhs);
            } else if (operation == "minPICCS" ||
                       operation == "minpiccs") {
                initializeRuntime();
                runPiccs(nlhs, plhs, nrhs, prhs);
            } else {
                fail("TIGRE:Metal:UnknownOperation",
                     "Unknown Metal backend operation: " + operation);
            }
        } catch (const MetalError& exception) {
            recordError(exception.id, exception.what());
        } catch (const tigre_metal_limits::LimitExceeded& exception) {
            recordError("TIGRE:Metal:SizeLimitExceeded", exception.what());
        } catch (const std::bad_alloc&) {
            recordError("TIGRE:Metal:AllocationFailed",
                        "The Metal backend could not allocate its whole-array buffers.");
        } catch (const std::exception& exception) {
            recordError("TIGRE:Metal:RuntimeError", exception.what());
        }
    }
    if (errorId != nullptr) {
        mexErrMsgIdAndTxt(errorId, "%s", errorMessage);
    }
}
