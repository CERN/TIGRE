/*
 * Portable CPU backend for MATLAB TIGRE.
 *
 * Computation is float32 and uses only the C++ standard library for
 * threading. Per-view geometry is prepared in double precision on the host.
 */

#include "mex.h"
#include "matrix.h"
#include "../Geometry/tigre_geometry.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr float kEpsilon = 1.0e-8f;

using Geometry = tigre_geometry::Geometry;
using Vec3 = tigre_geometry::Vector;
using ViewGeometry = tigre_geometry::View;

struct Ray {
    Vec3 source;
    Vec3 detector;
};

void fail(const char* id, const std::string& message) {
    mexErrMsgIdAndTxt(id, "%s", message.c_str());
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

const float* singleInput(const mxArray* value, const char* description) {
    if (!mxIsSingle(value) || mxIsComplex(value)) {
        fail("TIGRE:CPU:InvalidInput",
             std::string(description) + " must be a real single array.");
    }
    return static_cast<const float*>(mxGetData(value));
}

template <typename Function>
void parallelFor(std::size_t count, Function function) {
    const unsigned int available = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t threadCount =
        count < 64 ? 1 : std::min<std::size_t>(available, count);
    if (threadCount == 1) {
        for (std::size_t index = 0; index < count; ++index) {
            function(index);
        }
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(threadCount);
    for (std::size_t worker = 0; worker < threadCount; ++worker) {
        const std::size_t begin = count * worker / threadCount;
        const std::size_t end = count * (worker + 1) / threadCount;
        workers.emplace_back([begin, end, &function]() {
            for (std::size_t index = begin; index < end; ++index) {
                function(index);
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
}

inline Vec3 subtract(const Vec3& left, const Vec3& right) {
    return {left.x - right.x, left.y - right.y, left.z - right.z};
}

inline Vec3 addScaled(const Vec3& point, const Vec3& direction, float scale) {
    return {point.x + direction.x * scale,
            point.y + direction.y * scale,
            point.z + direction.z * scale};
}

inline float length(const Vec3& value) {
    return std::sqrt(value.x * value.x + value.y * value.y +
                     value.z * value.z);
}

inline float dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

Ray makeRay(const Geometry& geometry, const ViewGeometry& view,
            std::size_t detectorU, std::size_t detectorV) {
    detectorV = geometry.nv - detectorV - 1;
    Ray ray{};
    ray.detector = addScaled(
        addScaled(view.detector, view.detectorU,
                  static_cast<float>(detectorU)),
        view.detectorV, static_cast<float>(detectorV));
    if (geometry.parallel) {
        ray.source = addScaled(
            addScaled(view.source, view.detectorU,
                      static_cast<float>(detectorU)),
            view.detectorV, static_cast<float>(detectorV));
    } else {
        ray.source = view.source;
    }
    return ray;
}

bool detectorCoordinates(const Geometry& geometry, const ViewGeometry& view,
                         const Vec3& voxelPosition, float& u, float& v,
                         float& sourceDistanceSquared, float& rayLength,
                         float& axialDistance) {
    const Vec3 planeDirection = subtract(view.detector, view.source);
    Vec3 detectorPoint{};
    if (geometry.parallel) {
        const float denominator = dot(planeDirection, view.detectorNormal);
        if (std::abs(denominator) <= kEpsilon) {
            return false;
        }
        const float scale =
            dot(subtract(view.detector, voxelPosition), view.detectorNormal) /
            denominator;
        detectorPoint = addScaled(voxelPosition, planeDirection, scale);
        axialDistance = 1.0f;
        sourceDistanceSquared = 1.0f;
        rayLength = length(planeDirection);
    } else {
        const Vec3 sourceToVoxel = subtract(voxelPosition, view.source);
        const float denominator = dot(sourceToVoxel, view.detectorNormal);
        const float detectorDistance =
            dot(planeDirection, view.detectorNormal);
        if (std::abs(denominator) <= kEpsilon ||
            detectorDistance / denominator <= 0.0f) {
            return false;
        }
        detectorPoint = addScaled(
            view.source, sourceToVoxel, detectorDistance / denominator);
        sourceDistanceSquared = dot(sourceToVoxel, sourceToVoxel);
        if (sourceDistanceSquared <= kEpsilon) {
            return false;
        }
        rayLength = length(subtract(detectorPoint, view.source));
        axialDistance = denominator;
        if (axialDistance <= kEpsilon) {
            return false;
        }
    }

    const Vec3 relative = subtract(detectorPoint, view.detector);
    const float uLengthSquared = dot(view.detectorU, view.detectorU);
    const float vLengthSquared = dot(view.detectorV, view.detectorV);
    if (uLengthSquared <= kEpsilon || vLengthSquared <= kEpsilon) {
        return false;
    }
    u = dot(relative, view.detectorU) / uLengthSquared;
    v = dot(relative, view.detectorV) / vLengthSquared;
    return true;
}

bool intersectAxis(float source, float direction, float lower, float upper,
                   float& entry, float& exit) {
    if (std::abs(direction) < 1.0e-20f) {
        return source >= lower && source <= upper;
    }
    float first = (lower - source) / direction;
    float second = (upper - source) / direction;
    if (first > second) {
        std::swap(first, second);
    }
    entry = std::max(entry, first);
    exit = std::min(exit, second);
    return entry <= exit;
}

bool intersectVolume(const Geometry& geometry, const Ray& ray,
                     float& entry, float& exit) {
    const Vec3 direction = subtract(ray.detector, ray.source);
    entry = 0.0f;
    exit = 1.0f;
    return intersectAxis(ray.source.x, direction.x, -geometry.sx * 0.5f,
                         geometry.sx * 0.5f, entry, exit) &&
           intersectAxis(ray.source.y, direction.y, -geometry.sy * 0.5f,
                         geometry.sy * 0.5f, entry, exit) &&
           intersectAxis(ray.source.z, direction.z, -geometry.sz * 0.5f,
                         geometry.sz * 0.5f, entry, exit) &&
           entry <= exit;
}

inline float voxel(const float* image, const Geometry& geometry,
                   int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 ||
        x >= static_cast<int>(geometry.nx) ||
        y >= static_cast<int>(geometry.ny) ||
        z >= static_cast<int>(geometry.nz)) {
        return 0.0f;
    }
    return image[static_cast<std::size_t>(x) +
                 static_cast<std::size_t>(y) * geometry.nx +
                 static_cast<std::size_t>(z) * geometry.nx * geometry.ny];
}

float trilinear(const float* image, const Geometry& geometry,
                const Vec3& point) {
    const float x = (point.x + geometry.sx * 0.5f) / geometry.dx - 0.5f;
    const float y = (point.y + geometry.sy * 0.5f) / geometry.dy - 0.5f;
    const float z = (point.z + geometry.sz * 0.5f) / geometry.dz - 0.5f;
    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));
    const int z0 = static_cast<int>(std::floor(z));
    const float fx = x - static_cast<float>(x0);
    const float fy = y - static_cast<float>(y0);
    const float fz = z - static_cast<float>(z0);

    const float c00 = voxel(image, geometry, x0, y0, z0) * (1.0f - fx) +
                      voxel(image, geometry, x0 + 1, y0, z0) * fx;
    const float c10 = voxel(image, geometry, x0, y0 + 1, z0) * (1.0f - fx) +
                      voxel(image, geometry, x0 + 1, y0 + 1, z0) * fx;
    const float c01 = voxel(image, geometry, x0, y0, z0 + 1) * (1.0f - fx) +
                      voxel(image, geometry, x0 + 1, y0, z0 + 1) * fx;
    const float c11 =
        voxel(image, geometry, x0, y0 + 1, z0 + 1) * (1.0f - fx) +
        voxel(image, geometry, x0 + 1, y0 + 1, z0 + 1) * fx;
    const float c0 = c00 * (1.0f - fy) + c10 * fy;
    const float c1 = c01 * (1.0f - fy) + c11 * fy;
    return c0 * (1.0f - fz) + c1 * fz;
}

float interpolatedRay(const float* image, const Geometry& geometry,
                      const Ray& ray) {
    float entry;
    float exit;
    if (!intersectVolume(geometry, ray, entry, exit)) {
        return 0.0f;
    }

    const Vec3 physicalDirection = subtract(ray.detector, ray.source);
    const Vec3 indexDirection{
        physicalDirection.x / geometry.dx,
        physicalDirection.y / geometry.dy,
        physicalDirection.z / geometry.dz};
    const std::size_t sampleCount = std::max<std::size_t>(
        1, static_cast<std::size_t>(
               std::ceil(length(indexDirection) / geometry.accuracy)));
    const std::size_t first = static_cast<std::size_t>(std::max(
        0.0f, std::floor(entry * static_cast<float>(sampleCount))));
    const std::size_t last = std::min<std::size_t>(
        sampleCount,
        static_cast<std::size_t>(
            std::ceil(exit * static_cast<float>(sampleCount))));

    float sum = 0.0f;
    for (std::size_t sample = first; sample <= last; ++sample) {
        const float parameter =
            static_cast<float>(sample) / static_cast<float>(sampleCount);
        sum += trilinear(image, geometry,
                         addScaled(ray.source, physicalDirection, parameter));
    }
    return sum * length(physicalDirection) /
           static_cast<float>(sampleCount);
}

int voxelIndex(float coordinate, float size, float spacing,
               std::size_t count) {
    int index =
        static_cast<int>(std::floor((coordinate + size * 0.5f) / spacing));
    return std::max(0, std::min(index, static_cast<int>(count) - 1));
}

float nextBoundaryParameter(float coordinate, float direction, float size,
                            float spacing, int index) {
    if (std::abs(direction) < 1.0e-20f) {
        return std::numeric_limits<float>::infinity();
    }
    const float boundary = direction > 0.0f
        ? -size * 0.5f + static_cast<float>(index + 1) * spacing
        : -size * 0.5f + static_cast<float>(index) * spacing;
    return (boundary - coordinate) / direction;
}

float siddonRay(const float* image, const Geometry& geometry, const Ray& ray) {
    float entry;
    float exit;
    if (!intersectVolume(geometry, ray, entry, exit) || exit <= entry) {
        return 0.0f;
    }

    const Vec3 direction = subtract(ray.detector, ray.source);
    const float rayLength = length(direction);
    const float inside = std::nextafter(entry, exit);
    const Vec3 point = addScaled(ray.source, direction, inside);
    int ix = voxelIndex(point.x, geometry.sx, geometry.dx, geometry.nx);
    int iy = voxelIndex(point.y, geometry.sy, geometry.dy, geometry.ny);
    int iz = voxelIndex(point.z, geometry.sz, geometry.dz, geometry.nz);
    const int stepX = direction.x >= 0.0f ? 1 : -1;
    const int stepY = direction.y >= 0.0f ? 1 : -1;
    const int stepZ = direction.z >= 0.0f ? 1 : -1;
    const float deltaX = std::abs(direction.x) < 1.0e-20f
        ? std::numeric_limits<float>::infinity()
        : geometry.dx / std::abs(direction.x);
    const float deltaY = std::abs(direction.y) < 1.0e-20f
        ? std::numeric_limits<float>::infinity()
        : geometry.dy / std::abs(direction.y);
    const float deltaZ = std::abs(direction.z) < 1.0e-20f
        ? std::numeric_limits<float>::infinity()
        : geometry.dz / std::abs(direction.z);
    float nextX = nextBoundaryParameter(
        ray.source.x, direction.x, geometry.sx, geometry.dx, ix);
    float nextY = nextBoundaryParameter(
        ray.source.y, direction.y, geometry.sy, geometry.dy, iy);
    float nextZ = nextBoundaryParameter(
        ray.source.z, direction.z, geometry.sz, geometry.dz, iz);

    float current = entry;
    float sum = 0.0f;
    while (current < exit &&
           ix >= 0 && ix < static_cast<int>(geometry.nx) &&
           iy >= 0 && iy < static_cast<int>(geometry.ny) &&
           iz >= 0 && iz < static_cast<int>(geometry.nz)) {
        const float next = std::min(exit, std::min(nextX, std::min(nextY, nextZ)));
        sum += voxel(image, geometry, ix, iy, iz) *
               std::max(0.0f, next - current) * rayLength;
        current = next;
        const float tolerance = 8.0f * std::numeric_limits<float>::epsilon() *
                                std::max(1.0f, std::abs(current));
        if (nextX <= current + tolerance) {
            ix += stepX;
            nextX += deltaX;
        }
        if (nextY <= current + tolerance) {
            iy += stepY;
            nextY += deltaY;
        }
        if (nextZ <= current + tolerance) {
            iz += stepZ;
            nextZ += deltaZ;
        }
    }
    return sum;
}

float projectionSample(const float* projections, const Geometry& geometry,
                       std::size_t angle, float u, float v) {
    const int u0 = static_cast<int>(std::floor(u));
    const int v0 = static_cast<int>(std::floor(v));
    const float fu = u - static_cast<float>(u0);
    const float fv = v - static_cast<float>(v0);
    const std::size_t plane = geometry.nu * geometry.nv;
    auto fetch = [&](int ui, int vi) {
        if (ui < 0 || vi < 0 || ui >= static_cast<int>(geometry.nu) ||
            vi >= static_cast<int>(geometry.nv)) {
            return 0.0f;
        }
        return projections[static_cast<std::size_t>(vi) +
                           static_cast<std::size_t>(ui) * geometry.nv +
                           angle * plane];
    };
    const float lower = fetch(u0, v0) * (1.0f - fu) +
                        fetch(u0 + 1, v0) * fu;
    const float upper = fetch(u0, v0 + 1) * (1.0f - fu) +
                        fetch(u0 + 1, v0 + 1) * fu;
    return lower * (1.0f - fv) + upper * fv;
}

void runAx(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5 || nlhs > 1) {
        fail("TIGRE:CPU:Ax:InvalidInput",
             "Ax requires image, geometry, angles, and projector type.");
    }
    const float* image = singleInput(prhs[1], "Image");
    const std::string projector =
        getString(prhs[4], "TIGRE:CPU:Ax:InvalidInput", "Projector type");
    const bool siddon = projector == "Siddon" || projector == "ray-voxel";
    if (!siddon && projector != "interpolated") {
        fail("TIGRE:CPU:Ax:InvalidInput", "Unknown projector type.");
    }
    const Geometry geometry = tigre_geometry::parse(
        prhs[2], prhs[3], "TIGRE:CPU:BadGeometry", !siddon);
    const std::size_t angleCount = geometry.views.size();
    if (mxGetNumberOfElements(prhs[1]) !=
        geometry.nx * geometry.ny * geometry.nz) {
        fail("TIGRE:CPU:Ax:InvalidInput",
             "Image dimensions do not match geo.nVoxel.");
    }

    const mwSize dimensions[3] = {
        static_cast<mwSize>(geometry.nv),
        static_cast<mwSize>(geometry.nu),
        static_cast<mwSize>(angleCount)};
    plhs[0] = mxCreateNumericArray(3, dimensions, mxSINGLE_CLASS, mxREAL);
    float* output = static_cast<float*>(mxGetData(plhs[0]));
    const std::size_t detectorPixels = geometry.nu * geometry.nv;
    parallelFor(detectorPixels * angleCount, [&](std::size_t index) {
        const std::size_t angleIndex = index / detectorPixels;
        const std::size_t within = index - angleIndex * detectorPixels;
        const std::size_t u = within / geometry.nv;
        const std::size_t v = within - u * geometry.nv;
        const Ray ray = makeRay(geometry, geometry.views[angleIndex], u, v);
        output[index] = siddon
            ? siddonRay(image, geometry, ray)
            : interpolatedRay(image, geometry, ray);
    });
}

void runAtb(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5 || nlhs > 1) {
        fail("TIGRE:CPU:Atb:InvalidInput",
             "Atb requires projections, geometry, angles, and backprojector type.");
    }
    const float* projections = singleInput(prhs[1], "Projections");
    const Geometry geometry =
        tigre_geometry::parse(prhs[2], prhs[3], "TIGRE:CPU:BadGeometry");
    const std::size_t angleCount = geometry.views.size();
    const std::string projector =
        getString(prhs[4], "TIGRE:CPU:Atb:InvalidInput", "Backprojector type");
    const bool matched = projector == "matched";
    if (!matched && projector != "FDK") {
        fail("TIGRE:CPU:Atb:InvalidInput", "Unknown backprojector type.");
    }
    if (mxGetNumberOfElements(prhs[1]) !=
        geometry.nu * geometry.nv * angleCount) {
        fail("TIGRE:CPU:Atb:InvalidInput",
             "Projection dimensions do not match the geometry and angles.");
    }

    const mwSize dimensions[3] = {
        static_cast<mwSize>(geometry.nx),
        static_cast<mwSize>(geometry.ny),
        static_cast<mwSize>(geometry.nz)};
    plhs[0] = mxCreateNumericArray(3, dimensions, mxSINGLE_CLASS, mxREAL);
    float* output = static_cast<float*>(mxGetData(plhs[0]));
    const std::size_t voxelCount = geometry.nx * geometry.ny * geometry.nz;
    const float matchedScale =
        geometry.dx * geometry.dy * geometry.dz /
        (geometry.du * geometry.dv);

    parallelFor(voxelCount, [&](std::size_t index) {
        const std::size_t plane = geometry.nx * geometry.ny;
        const std::size_t zIndex = index / plane;
        const std::size_t within = index - zIndex * plane;
        const std::size_t yIndex = within / geometry.nx;
        const std::size_t xIndex = within - yIndex * geometry.nx;
        const float x =
            -geometry.sx * 0.5f + geometry.dx * 0.5f +
            static_cast<float>(xIndex) * geometry.dx;
        const float y =
            -geometry.sy * 0.5f + geometry.dy * 0.5f +
            static_cast<float>(yIndex) * geometry.dy;
        const float z =
            -geometry.sz * 0.5f + geometry.dz * 0.5f +
            static_cast<float>(zIndex) * geometry.dz;

        float sum = 0.0f;
        for (std::size_t angleIndex = 0; angleIndex < angleCount; ++angleIndex) {
            const ViewGeometry& view = geometry.views[angleIndex];
            // CUDA Atb transforms voxels into a detector-scaled coordinate
            // system before sampling. Reproduce that path, including its
            // historical weighting formulas, rather than deriving an inverse
            // directly from the forward-projector rays.
            Vec3 point = addScaled(
                addScaled(
                    addScaled(view.backprojectionOrigin,
                              view.backprojectionDeltaX,
                              static_cast<float>(xIndex)),
                    view.backprojectionDeltaY,
                    static_cast<float>(yIndex)),
                view.backprojectionDeltaZ,
                static_cast<float>(zIndex));
            point.y -= view.originOffset.padding / geometry.du;

            float sample = 0.0f;
            if (geometry.parallel) {
                const float u =
                    point.y + static_cast<float>(geometry.nu) * 0.5f - 0.5f;
                const float v =
                    point.z + static_cast<float>(geometry.nv) * 0.5f - 0.5f;
                sample = projectionSample(
                    projections, geometry, angleIndex, u, v);
            } else {
                const Vec3 direction =
                    subtract(point, view.backprojectionSource);
                if (std::abs(direction.x) <= kEpsilon) {
                    continue;
                }
                const float t =
                    (view.dso - view.dsd - view.backprojectionSource.x) /
                    direction.x;
                const float u =
                    direction.y * t + view.backprojectionSource.y +
                    static_cast<float>(geometry.nu) * 0.5f;
                const float v =
                    direction.z * t + view.backprojectionSource.z +
                    static_cast<float>(geometry.nv) * 0.5f;
                sample = projectionSample(
                    projections, geometry, angleIndex, u - 0.5f, v - 0.5f);

                const float realX = x + view.originOffset.x;
                const float realY = y + view.originOffset.y;
                const float realZ = z + view.originOffset.z;
                const float sinAlpha = view.detectorParameters.z;
                const float cosAlpha = view.detectorParameters.padding;
                if (matched) {
                    const float realSourceX = view.dso * cosAlpha;
                    const float realSourceY = -view.dso * sinAlpha;
                    const float detectorX = -(view.dsd - view.dso);
                    const float detectorY =
                        (-geometry.su + geometry.du) * 0.5f +
                        u * geometry.du + view.detectorParameters.x;
                    const float realDetectorX =
                        detectorX * cosAlpha + detectorY * sinAlpha;
                    const float realDetectorY =
                        -detectorX * sinAlpha + detectorY * cosAlpha;
                    const float realDetectorZ =
                        (-geometry.sv + geometry.dv) * 0.5f +
                        v * geometry.dv + view.detectorParameters.y;
                    const float dx = realSourceX - realDetectorX;
                    const float dy = realSourceY - realDetectorY;
                    const float rayLength =
                        std::sqrt(dx * dx + dy * dy +
                                  realDetectorZ * realDetectorZ);
                    const float sourceX = realSourceX - realX;
                    const float sourceY = realSourceY - realY;
                    const float sourceDistanceSquared =
                        sourceX * sourceX + sourceY * sourceY +
                        realZ * realZ;
                    sample *= rayLength * rayLength * rayLength /
                              (view.dsd * sourceDistanceSquared);
                } else {
                    const float distanceRatio =
                        (view.dso +
                         (realY + view.originOffset.padding) * sinAlpha -
                         realX * cosAlpha) /
                        view.dso;
                    sample /= distanceRatio * distanceRatio;
                }
            }
            sum += sample;
        }
        output[index] =
            (!geometry.parallel && matched) ? sum * matchedScale : sum;
    });
}

struct ImageShape {
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
    std::size_t count;
    mwSize dimensions[3];
    mwSize dimensionCount;
};

ImageShape imageShape(const mxArray* image) {
    singleInput(image, "Image");
    const mwSize numberOfDimensions = mxGetNumberOfDimensions(image);
    if (numberOfDimensions < 2 || numberOfDimensions > 3) {
        fail("TIGRE:CPU:TV:InvalidInput", "Image must be 2-D or 3-D.");
    }
    const mwSize* sourceDimensions = mxGetDimensions(image);
    ImageShape shape{};
    shape.nx = sourceDimensions[0];
    shape.ny = sourceDimensions[1];
    shape.nz = numberOfDimensions == 3 ? sourceDimensions[2] : 1;
    shape.count = shape.nx * shape.ny * shape.nz;
    shape.dimensions[0] = sourceDimensions[0];
    shape.dimensions[1] = sourceDimensions[1];
    shape.dimensions[2] = shape.nz;
    shape.dimensionCount = numberOfDimensions;
    return shape;
}

mxArray* createLike(const ImageShape& shape) {
    return mxCreateNumericArray(
        shape.dimensionCount, shape.dimensions, mxSINGLE_CLASS, mxREAL);
}

inline std::size_t imageIndex(const ImageShape& shape,
                              std::size_t x, std::size_t y, std::size_t z) {
    return x + y * shape.nx + z * shape.nx * shape.ny;
}

void backwardGradient(const float* image, const ImageShape& shape,
                      int x, int y, int z, float gradient[3]) {
    gradient[0] = 0.0f;
    gradient[1] = 0.0f;
    gradient[2] = 0.0f;
    if (x < 0 || y < 0 || z < 0 ||
        x >= static_cast<int>(shape.nx) ||
        y >= static_cast<int>(shape.ny) ||
        z >= static_cast<int>(shape.nz)) {
        return;
    }
    const std::size_t index = imageIndex(
        shape, static_cast<std::size_t>(x), static_cast<std::size_t>(y),
        static_cast<std::size_t>(z));
    const float value = image[index];
    if (z > 0) {
        gradient[0] = value - image[index - shape.nx * shape.ny];
    }
    if (y > 0) {
        gradient[1] = value - image[index - shape.nx];
    }
    if (x > 0) {
        gradient[2] = value - image[index - 1];
    }
}

float tvGradientAt(const float* image, const ImageShape& shape,
                   std::size_t x, std::size_t y, std::size_t z,
                   bool adaptive, float delta, float epsilon) {
    float current[3];
    float nextX[3];
    float nextY[3];
    float nextZ[3];
    backwardGradient(image, shape, static_cast<int>(x), static_cast<int>(y),
                     static_cast<int>(z), current);
    backwardGradient(image, shape, static_cast<int>(x + 1),
                     static_cast<int>(y), static_cast<int>(z), nextX);
    backwardGradient(image, shape, static_cast<int>(x),
                     static_cast<int>(y + 1), static_cast<int>(z), nextY);
    backwardGradient(image, shape, static_cast<int>(x),
                     static_cast<int>(y), static_cast<int>(z + 1), nextZ);

    auto weights = [&](const float values[3], float result[3]) {
        if (adaptive) {
            const float denominator = std::abs(delta) > kEpsilon
                ? delta * delta
                : kEpsilon;
            for (int component = 0; component < 3; ++component) {
                result[component] =
                    std::exp(-(values[component] * values[component]) /
                             denominator);
            }
        } else {
            result[0] = 1.0f;
            result[1] = 1.0f;
            result[2] = 1.0f;
        }
    };
    auto weightedNorm = [&](const float values[3], const float weight[3]) {
        return std::sqrt(weight[0] * values[0] * values[0] +
                         weight[1] * values[1] * values[1] +
                         weight[2] * values[2] * values[2]) +
               epsilon;
    };

    float weight[3];
    float weightX[3];
    float weightY[3];
    float weightZ[3];
    weights(current, weight);
    weights(nextX, weightX);
    weights(nextY, weightY);
    weights(nextZ, weightZ);
    return (weight[0] * current[0] + weight[1] * current[1] +
            weight[2] * current[2]) /
               weightedNorm(current, weight) -
           weightX[2] * nextX[2] / weightedNorm(nextX, weightX) -
           weightY[1] * nextY[1] / weightedNorm(nextY, weightY) -
           weightZ[0] * nextZ[0] / weightedNorm(nextZ, weightZ);
}

void computeTvGradient(const float* image, float* gradient,
                       const ImageShape& shape, bool adaptive,
                       float delta, float epsilon) {
    parallelFor(shape.count, [&](std::size_t index) {
        const std::size_t plane = shape.nx * shape.ny;
        const std::size_t z = index / plane;
        const std::size_t within = index - z * plane;
        const std::size_t y = within / shape.nx;
        const std::size_t x = within - y * shape.nx;
        gradient[index] =
            tvGradientAt(image, shape, x, y, z, adaptive, delta, epsilon);
    });
}

double normSquared(const float* values, std::size_t count) {
    const unsigned int available = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t threadCount =
        count < 4096 ? 1 : std::min<std::size_t>(available, count);
    std::vector<double> partial(threadCount, 0.0);
    parallelFor(threadCount, [&](std::size_t worker) {
        const std::size_t begin = count * worker / threadCount;
        const std::size_t end = count * (worker + 1) / threadCount;
        double sum = 0.0;
        for (std::size_t index = begin; index < end; ++index) {
            const double value = values[index];
            sum += value * value;
        }
        partial[worker] = sum;
    });
    double total = 0.0;
    for (double value : partial) {
        total += value;
    }
    return total;
}

int scalarInteger(const mxArray* value, const char* name) {
    if (!mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1) {
        fail("TIGRE:CPU:TV:InvalidInput",
             std::string(name) + " must be a real scalar.");
    }
    const double scalar = mxGetScalar(value);
    if (!std::isfinite(scalar) || scalar < 0.0) {
        fail("TIGRE:CPU:TV:InvalidInput",
             std::string(name) + " must be finite and nonnegative.");
    }
    return static_cast<int>(std::floor(scalar + 0.5));
}

float scalarFloat(const mxArray* value, const char* name) {
    if (!mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1) {
        fail("TIGRE:CPU:TV:InvalidInput",
             std::string(name) + " must be a real scalar.");
    }
    const double scalar = mxGetScalar(value);
    if (!std::isfinite(scalar)) {
        fail("TIGRE:CPU:TV:InvalidInput",
             std::string(name) + " must be finite.");
    }
    return static_cast<float>(scalar);
}

void gradientDescentTv(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[], bool adaptive) {
    const int expectedInputs = adaptive ? 5 : 4;
    if (nrhs != expectedInputs || nlhs > 1) {
        fail("TIGRE:CPU:TV:InvalidInput", "Invalid TV operation inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const ImageShape shape = imageShape(prhs[1]);
    const float alpha = scalarFloat(prhs[2], "Step size");
    const int iterations = scalarInteger(prhs[3], "Iteration count");
    const float delta = adaptive ? scalarFloat(prhs[4], "Delta") : 1.0f;
    plhs[0] = createLike(shape);
    float* output = static_cast<float*>(mxGetData(plhs[0]));
    std::memcpy(output, input, shape.count * sizeof(float));
    std::vector<float> gradient(shape.count);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        computeTvGradient(
            output, gradient.data(), shape, adaptive, delta, 1.0e-8f);
        const float norm = static_cast<float>(
            std::sqrt(normSquared(gradient.data(), shape.count) + 1.0e-7));
        const float scale = alpha / norm;
        parallelFor(shape.count, [&](std::size_t index) {
            output[index] -= scale * gradient[index];
        });
    }
}

float divergenceAt(const std::vector<float>& pz,
                   const std::vector<float>& py,
                   const std::vector<float>& px,
                   const ImageShape& shape,
                   std::size_t x, std::size_t y, std::size_t z) {
    const std::size_t index = imageIndex(shape, x, y, z);
    float divergence = pz[index] + py[index] + px[index];
    if (z > 0) {
        divergence -= pz[index - shape.nx * shape.ny];
    }
    if (y > 0) {
        divergence -= py[index - shape.nx];
    }
    if (x > 0) {
        divergence -= px[index - 1];
    }
    return divergence;
}

void runTvDenoise(int nlhs, mxArray* plhs[], int nrhs,
                  const mxArray* prhs[]) {
    if (nrhs != 4 || nlhs > 1) {
        fail("TIGRE:CPU:TV:InvalidInput", "Invalid tvDenoise inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const ImageShape shape = imageShape(prhs[1]);
    if (shape.dimensionCount != 3) {
        fail("TIGRE:CPU:TV:InvalidInput", "tvDenoise requires a 3-D image.");
    }
    const float lambda = scalarFloat(prhs[2], "Lambda");
    const int iterations = scalarInteger(prhs[3], "Iteration count");
    if (std::abs(lambda) <= kEpsilon) {
        fail("TIGRE:CPU:TV:InvalidInput", "Lambda must be nonzero.");
    }

    plhs[0] = createLike(shape);
    float* output = static_cast<float*>(mxGetData(plhs[0]));
    std::memcpy(output, input, shape.count * sizeof(float));
    std::vector<float> pz(shape.count, 0.0f);
    std::vector<float> py(shape.count, 0.0f);
    std::vector<float> px(shape.count, 0.0f);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        const float tau2 = 0.3f + 0.02f * static_cast<float>(iteration);
        const float tau1 =
            (1.0f / tau2) *
            (1.0f / 6.0f - 5.0f / (15.0f + static_cast<float>(iteration)));
        parallelFor(shape.count, [&](std::size_t index) {
            const std::size_t plane = shape.nx * shape.ny;
            const std::size_t z = index / plane;
            const std::size_t within = index - z * plane;
            const std::size_t y = within / shape.nx;
            const std::size_t x = within - y * shape.nx;
            const float divergence =
                divergenceAt(pz, py, px, shape, x, y, z);
            output[index] = output[index] * (1.0f - tau1) +
                            tau1 * (input[index] + divergence / lambda);
        });
        parallelFor(shape.count, [&](std::size_t index) {
            const std::size_t plane = shape.nx * shape.ny;
            const std::size_t z = index / plane;
            const std::size_t within = index - z * plane;
            const std::size_t y = within / shape.nx;
            const std::size_t x = within - y * shape.nx;
            const float value = output[index];
            const float gz = z + 1 < shape.nz
                ? output[index + plane] - value
                : 0.0f;
            const float gy = y + 1 < shape.ny
                ? output[index + shape.nx] - value
                : 0.0f;
            const float gx = x + 1 < shape.nx
                ? output[index + 1] - value
                : 0.0f;
            const float qz = pz[index] + tau2 * gz;
            const float qy = py[index] + tau2 * gy;
            const float qx = px[index] + tau2 * gx;
            const float norm =
                std::max(1.0f, std::sqrt(qz * qz + qy * qy + qx * qx));
            pz[index] = qz / norm;
            py[index] = qy / norm;
            px[index] = qx / norm;
        });
    }
}

void runPiccs(int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (nrhs != 6 || nlhs > 1) {
        fail("TIGRE:CPU:TV:InvalidInput", "Invalid minPICCS inputs.");
    }
    const float* input = singleInput(prhs[1], "Image");
    const float* prior = singleInput(prhs[2], "Prior image");
    const ImageShape shape = imageShape(prhs[1]);
    const ImageShape priorShape = imageShape(prhs[2]);
    if (shape.dimensionCount != 3 || priorShape.dimensionCount != 3 ||
        shape.nx != priorShape.nx || shape.ny != priorShape.ny ||
        shape.nz != priorShape.nz) {
        fail("TIGRE:CPU:TV:InvalidInput",
             "Image and prior must be equally sized 3-D arrays.");
    }
    const float alpha = scalarFloat(prhs[3], "Step size");
    const int iterations = scalarInteger(prhs[4], "Iteration count");
    const float ratio = scalarFloat(prhs[5], "PICCS ratio");
    plhs[0] = createLike(shape);
    float* output = static_cast<float*>(mxGetData(plhs[0]));
    std::memcpy(output, input, shape.count * sizeof(float));
    std::vector<float> difference(shape.count);
    std::vector<float> gradient(shape.count);
    std::vector<float> priorGradient(shape.count);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        computeTvGradient(
            output, gradient.data(), shape, false, 1.0f, 1.0e-6f);
        parallelFor(shape.count, [&](std::size_t index) {
            difference[index] = output[index] - prior[index];
        });
        computeTvGradient(difference.data(), priorGradient.data(), shape,
                          false, 1.0f, 1.0e-6f);
        parallelFor(shape.count, [&](std::size_t index) {
            gradient[index] =
                (1.0f - ratio) * gradient[index] +
                ratio * priorGradient[index];
        });
        const double squaredNorm = normSquared(gradient.data(), shape.count);
        if (squaredNorm <= 0.0) {
            continue;
        }
        const float scale =
            alpha / static_cast<float>(std::sqrt(squaredNorm));
        parallelFor(shape.count, [&](std::size_t index) {
            output[index] -= scale * gradient[index];
        });
    }
}

}  // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]) {
    try {
        if (nrhs < 1) {
            fail("TIGRE:CPU:InvalidInput", "An operation name is required.");
        }
        const std::string operation =
            getString(prhs[0], "TIGRE:CPU:InvalidInput", "Operation");
        if (operation == "Ax" || operation == "ax") {
            runAx(nlhs, plhs, nrhs, prhs);
        } else if (operation == "Atb" || operation == "atb") {
            runAtb(nlhs, plhs, nrhs, prhs);
        } else if (operation == "minTV" || operation == "mintv") {
            gradientDescentTv(nlhs, plhs, nrhs, prhs, false);
        } else if (operation == "AwminTV" || operation == "awmintv") {
            gradientDescentTv(nlhs, plhs, nrhs, prhs, true);
        } else if (operation == "tvDenoise" || operation == "tvdenoise") {
            runTvDenoise(nlhs, plhs, nrhs, prhs);
        } else if (operation == "minPICCS" || operation == "minpiccs") {
            runPiccs(nlhs, plhs, nrhs, prhs);
        } else if (operation == "deviceName") {
            if (nrhs != 1 || nlhs > 1) {
                fail("TIGRE:CPU:InvalidInput",
                     "deviceName takes no additional inputs.");
            }
            plhs[0] = mxCreateString("CPU");
        } else {
            fail("TIGRE:CPU:UnknownOperation",
                 "Unknown CPU backend operation: " + operation);
        }
    } catch (const std::bad_alloc&) {
        fail("TIGRE:CPU:AllocationFailed",
             "The CPU backend could not allocate its whole-array buffers.");
    } catch (const std::exception& exception) {
        fail("TIGRE:CPU:RuntimeError", exception.what());
    }
}
