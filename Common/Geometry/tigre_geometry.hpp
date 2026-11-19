#ifndef TIGRE_BACKEND_GEOMETRY_HPP
#define TIGRE_BACKEND_GEOMETRY_HPP

#include "mex.h"
#include "matrix.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace tigre_geometry {

struct alignas(16) Vector {
    float x;
    float y;
    float z;
    float padding;
};

struct alignas(16) View {
    Vector source;
    Vector detector;
    Vector detectorU;
    Vector detectorV;
    Vector detectorNormal;
    float dso;
    float dsd;
    std::uint32_t parallel;
    float padding;
    // Precomputed in the coordinate system used by TIGRE's CUDA Atb kernels.
    // These preserve legacy sampling and weighting behavior for compatibility.
    Vector backprojectionOrigin;
    Vector backprojectionDeltaX;
    Vector backprojectionDeltaY;
    Vector backprojectionDeltaZ;
    Vector backprojectionSource;
    Vector originOffset;
    Vector detectorParameters;
};

static_assert(sizeof(Vector) == 16, "Unexpected geometry vector layout.");
static_assert(sizeof(View) == 208, "Unexpected per-view geometry layout.");

struct Geometry {
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
    std::size_t nu;
    std::size_t nv;
    float sx;
    float sy;
    float sz;
    float dx;
    float dy;
    float dz;
    float su;
    float sv;
    float du;
    float dv;
    float accuracy;
    bool parallel;
    std::vector<View> views;
};

namespace detail {

struct DoubleVector {
    double x;
    double y;
    double z;
};

inline void fail(const char* id, const std::string& message) {
    mexErrMsgIdAndTxt(id, "%s", message.c_str());
}

inline const mxArray* field(const mxArray* structure, const char* name,
                            const char* errorId) {
    const mxArray* value = mxGetField(structure, 0, name);
    if (value == nullptr) {
        fail(errorId, std::string("Missing geometry field geo.") + name + ".");
    }
    return value;
}

inline const double* doubles(const mxArray* value, const char* description,
                             const char* errorId) {
    if (!mxIsDouble(value) || mxIsComplex(value)) {
        fail(errorId, std::string(description) +
            " must be a real double array.");
    }
    return static_cast<const double*>(mxGetData(value));
}

inline std::string stringValue(const mxArray* value, const char* description,
                               const char* errorId) {
    if (!mxIsChar(value)) {
        fail(errorId, std::string(description) +
            " must be a character vector.");
    }
    char* raw = mxArrayToString(value);
    if (raw == nullptr) {
        fail(errorId, std::string("Unable to read ") + description + ".");
    }
    const std::string result(raw);
    mxFree(raw);
    return result;
}

inline std::size_t sizeValue(double value, const char* description,
                             const char* errorId) {
    if (!std::isfinite(value) || value < 1.0 ||
        value > static_cast<double>(
            std::numeric_limits<std::uint32_t>::max()) ||
        std::floor(value) != value) {
        fail(errorId, std::string(description) +
            " must contain positive integer sizes.");
    }
    return static_cast<std::size_t>(value);
}

inline double scalar(const mxArray* structure, const char* name,
                     const char* errorId) {
    const mxArray* value = field(structure, name, errorId);
    const double* data = doubles(value, name, errorId);
    if (mxGetNumberOfElements(value) != 1 || !std::isfinite(data[0])) {
        fail(errorId, std::string("geo.") + name +
            " must be a finite scalar.");
    }
    return data[0];
}

inline const double* vectorField(const mxArray* structure, const char* name,
                                 std::size_t components,
                                 std::size_t viewCount,
                                 const char* errorId) {
    const mxArray* value = field(structure, name, errorId);
    const double* data = doubles(value, name, errorId);
    if (mxGetM(value) != components ||
        (mxGetN(value) != 1 && mxGetN(value) != viewCount)) {
        fail(errorId, std::string("geo.") + name + " must be " +
            std::to_string(components) + "-by-1 or " +
            std::to_string(components) + "-by-number-of-views.");
    }
    const std::size_t count = mxGetNumberOfElements(value);
    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(data[index])) {
            fail(errorId, std::string("geo.") + name +
                " must contain finite values.");
        }
    }
    return data;
}

inline const double* scalarOrViewField(const mxArray* structure,
                                       const char* name,
                                       std::size_t viewCount,
                                       const char* errorId) {
    const mxArray* value = field(structure, name, errorId);
    const double* data = doubles(value, name, errorId);
    const std::size_t count = mxGetNumberOfElements(value);
    if (count != 1 && count != viewCount) {
        fail(errorId, std::string("geo.") + name +
            " must be scalar or contain one value per view.");
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(data[index])) {
            fail(errorId, std::string("geo.") + name +
                " must contain finite values.");
        }
    }
    return data;
}

inline double component(const mxArray* structure, const char* name,
                        const double* data, std::size_t components,
                        std::size_t componentIndex, std::size_t viewIndex,
                        const char* errorId) {
    const mxArray* value = field(structure, name, errorId);
    const std::size_t column = mxGetN(value) == 1 ? 0 : viewIndex;
    return data[componentIndex + components * column];
}

inline double viewValue(const mxArray* structure, const char* name,
                        const double* data, std::size_t viewIndex,
                        const char* errorId) {
    const mxArray* value = field(structure, name, errorId);
    return data[mxGetNumberOfElements(value) == 1 ? 0 : viewIndex];
}

inline DoubleVector add(DoubleVector first, DoubleVector second) {
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

inline DoubleVector subtract(DoubleVector first, DoubleVector second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

inline DoubleVector cross(DoubleVector first, DoubleVector second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

inline DoubleVector normalize(DoubleVector value, const char* errorId) {
    const double magnitude =
        std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
    if (!std::isfinite(magnitude) ||
        magnitude <= std::numeric_limits<double>::epsilon()) {
        fail(errorId, "Detector orientation is degenerate.");
    }
    return {value.x / magnitude, value.y / magnitude, value.z / magnitude};
}

inline DoubleVector detectorRotation(DoubleVector point, double roll,
                                     double pitch, double yaw) {
    const double sinRoll = std::sin(roll);
    const double cosRoll = std::cos(roll);
    const double sinPitch = std::sin(pitch);
    const double cosPitch = std::cos(pitch);
    const double sinYaw = std::sin(yaw);
    const double cosYaw = std::cos(yaw);
    return {
        cosRoll * cosPitch * point.x +
            (cosRoll * sinPitch * sinYaw - sinRoll * cosYaw) * point.y +
            (cosRoll * sinPitch * cosYaw + sinRoll * sinYaw) * point.z,
        sinRoll * cosPitch * point.x +
            (sinRoll * sinPitch * sinYaw + cosRoll * cosYaw) * point.y +
            (sinRoll * sinPitch * cosYaw - cosRoll * sinYaw) * point.z,
        -sinPitch * point.x + cosPitch * sinYaw * point.y +
            cosPitch * cosYaw * point.z};
}

inline DoubleVector detectorRotationTranspose(
        DoubleVector point, double roll, double pitch, double yaw) {
    const double sinRoll = std::sin(roll);
    const double cosRoll = std::cos(roll);
    const double sinPitch = std::sin(pitch);
    const double cosPitch = std::cos(pitch);
    const double sinYaw = std::sin(yaw);
    const double cosYaw = std::cos(yaw);
    return {
        cosRoll * cosPitch * point.x +
            sinRoll * cosPitch * point.y - sinPitch * point.z,
        (cosRoll * sinPitch * sinYaw - sinRoll * cosYaw) * point.x +
            (sinRoll * sinPitch * sinYaw + cosRoll * cosYaw) * point.y +
            cosPitch * sinYaw * point.z,
        (cosRoll * sinPitch * cosYaw + sinRoll * sinYaw) * point.x +
            (sinRoll * sinPitch * cosYaw - cosRoll * sinYaw) * point.y +
            cosPitch * cosYaw * point.z};
}

inline DoubleVector scanRotation(DoubleVector point, double alpha,
                                 double theta, double psi) {
    const double sinAlpha = std::sin(alpha);
    const double cosAlpha = std::cos(alpha);
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    const double sinPsi = std::sin(psi);
    const double cosPsi = std::cos(psi);
    return {
        (cosAlpha * cosTheta * cosPsi - sinAlpha * sinPsi) * point.x +
            (-cosAlpha * cosTheta * sinPsi - sinAlpha * cosPsi) * point.y +
            cosAlpha * sinTheta * point.z,
        (sinAlpha * cosTheta * cosPsi + cosAlpha * sinPsi) * point.x +
            (-sinAlpha * cosTheta * sinPsi + cosAlpha * cosPsi) * point.y +
            sinAlpha * sinTheta * point.z,
        -sinTheta * cosPsi * point.x +
            sinTheta * sinPsi * point.y + cosTheta * point.z};
}

inline DoubleVector scanRotationTranspose(DoubleVector point, double alpha,
                                          double theta, double psi) {
    const double sinAlpha = std::sin(alpha);
    const double cosAlpha = std::cos(alpha);
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    const double sinPsi = std::sin(psi);
    const double cosPsi = std::cos(psi);
    return {
        (cosPsi * cosTheta * cosAlpha - sinPsi * sinAlpha) * point.x +
            (-cosPsi * cosTheta * sinAlpha - sinPsi * cosAlpha) * point.y +
            cosPsi * sinTheta * point.z,
        (sinPsi * cosTheta * cosAlpha + cosPsi * sinAlpha) * point.x +
            (-sinPsi * cosTheta * sinAlpha + cosPsi * cosAlpha) * point.y +
            sinPsi * sinTheta * point.z,
        -sinTheta * cosAlpha * point.x +
            sinTheta * sinAlpha * point.y + cosTheta * point.z};
}

inline Vector packed(DoubleVector value, float padding = 0.0f) {
    return {static_cast<float>(value.x), static_cast<float>(value.y),
            static_cast<float>(value.z), padding};
}

}  // namespace detail

inline Geometry parse(const mxArray* value, const mxArray* angleArray,
                      const char* errorId,
                      bool cudaParallelInterpolatedOffsets = false) {
    using namespace detail;
    if (!mxIsStruct(value) || mxGetNumberOfElements(value) != 1) {
        fail(errorId, "Geometry must be a scalar structure.");
    }
    const double* angles = doubles(angleArray, "Angles", errorId);
    if (mxGetM(angleArray) != 3 || mxGetN(angleArray) == 0) {
        fail(errorId, "Angles must be a nonempty 3-by-N array.");
    }
    const std::size_t viewCount = mxGetN(angleArray);
    for (std::size_t index = 0; index < 3 * viewCount; ++index) {
        if (!std::isfinite(angles[index])) {
            fail(errorId, "Angles must contain finite values.");
        }
    }

    const mxArray* nVoxelArray = field(value, "nVoxel", errorId);
    const mxArray* sVoxelArray = field(value, "sVoxel", errorId);
    const mxArray* dVoxelArray = field(value, "dVoxel", errorId);
    const mxArray* nDetectorArray = field(value, "nDetector", errorId);
    const mxArray* sDetectorArray = field(value, "sDetector", errorId);
    const mxArray* dDetectorArray = field(value, "dDetector", errorId);
    if (mxGetNumberOfElements(nVoxelArray) != 3 ||
        mxGetNumberOfElements(sVoxelArray) != 3 ||
        mxGetNumberOfElements(dVoxelArray) != 3 ||
        mxGetNumberOfElements(nDetectorArray) != 2 ||
        mxGetNumberOfElements(sDetectorArray) != 2 ||
        mxGetNumberOfElements(dDetectorArray) != 2) {
        fail(errorId, "Geometry vector dimensions are invalid.");
    }
    const double* nVoxel = doubles(nVoxelArray, "geo.nVoxel", errorId);
    const double* sVoxel = doubles(sVoxelArray, "geo.sVoxel", errorId);
    const double* dVoxel = doubles(dVoxelArray, "geo.dVoxel", errorId);
    const double* nDetector =
        doubles(nDetectorArray, "geo.nDetector", errorId);
    const double* sDetector =
        doubles(sDetectorArray, "geo.sDetector", errorId);
    const double* dDetector =
        doubles(dDetectorArray, "geo.dDetector", errorId);

    Geometry geometry{};
    geometry.nx = sizeValue(nVoxel[0], "geo.nVoxel", errorId);
    geometry.ny = sizeValue(nVoxel[1], "geo.nVoxel", errorId);
    geometry.nz = sizeValue(nVoxel[2], "geo.nVoxel", errorId);
    geometry.nu = sizeValue(nDetector[0], "geo.nDetector", errorId);
    geometry.nv = sizeValue(nDetector[1], "geo.nDetector", errorId);
    geometry.sx = static_cast<float>(sVoxel[0]);
    geometry.sy = static_cast<float>(sVoxel[1]);
    geometry.sz = static_cast<float>(sVoxel[2]);
    geometry.dx = static_cast<float>(dVoxel[0]);
    geometry.dy = static_cast<float>(dVoxel[1]);
    geometry.dz = static_cast<float>(dVoxel[2]);
    geometry.su = static_cast<float>(sDetector[0]);
    geometry.sv = static_cast<float>(sDetector[1]);
    geometry.du = static_cast<float>(dDetector[0]);
    geometry.dv = static_cast<float>(dDetector[1]);
    geometry.accuracy = static_cast<float>(
        scalar(value, "accuracy", errorId));
    const std::string mode =
        stringValue(field(value, "mode", errorId), "geo.mode", errorId);
    if (mode != "cone" && mode != "parallel") {
        fail(errorId, "geo.mode must be 'cone' or 'parallel'.");
    }
    geometry.parallel = mode == "parallel";
    if (!std::isfinite(geometry.sx) || !std::isfinite(geometry.sy) ||
        !std::isfinite(geometry.sz) || !std::isfinite(geometry.dx) ||
        !std::isfinite(geometry.dy) || !std::isfinite(geometry.dz) ||
        !std::isfinite(geometry.su) || !std::isfinite(geometry.sv) ||
        !std::isfinite(geometry.du) || !std::isfinite(geometry.dv) ||
        !std::isfinite(geometry.accuracy) ||
        geometry.sx <= 0.0f || geometry.sy <= 0.0f ||
        geometry.sz <= 0.0f || geometry.dx <= 0.0f ||
        geometry.dy <= 0.0f || geometry.dz <= 0.0f ||
        geometry.su <= 0.0f || geometry.sv <= 0.0f ||
        geometry.du <= 0.0f || geometry.dv <= 0.0f ||
        geometry.accuracy < 0.001f) {
        fail(errorId, "Geometry dimensions and accuracy must be finite and positive.");
    }

    const double* dso =
        scalarOrViewField(value, "DSO", viewCount, errorId);
    const double* dsd =
        scalarOrViewField(value, "DSD", viewCount, errorId);
    const double* offOrigin =
        vectorField(value, "offOrigin", 3, viewCount, errorId);
    const double* offDetector =
        vectorField(value, "offDetector", 2, viewCount, errorId);
    const double* detectorRotationValues =
        vectorField(value, "rotDetector", 3, viewCount, errorId);
    const double* cor =
        scalarOrViewField(value, "COR", viewCount, errorId);

    geometry.views.reserve(viewCount);
    const double detectorU0 =
        -(static_cast<double>(geometry.nu) - 1.0) * 0.5 * geometry.du;
    const double detectorV0 =
        (static_cast<double>(geometry.nv) - 1.0) * 0.5 * geometry.dv;
    for (std::size_t index = 0; index < viewCount; ++index) {
        const double currentDso =
            viewValue(value, "DSO", dso, index, errorId);
        const double currentDsd =
            viewValue(value, "DSD", dsd, index, errorId);
        if (currentDso <= 0.0 || currentDsd < currentDso) {
            fail(errorId, "Each DSO must be positive and DSD must be at least DSO.");
        }
        const double alpha = angles[3 * index];
        const double theta = angles[3 * index + 1];
        const double psi = angles[3 * index + 2];
        // Preserve the MATLAB CUDA MEX mapping: rows 1/2/3 are stored in
        // the CUDA geometry's yaw/pitch/roll fields, respectively.
        const double yaw = component(value, "rotDetector",
            detectorRotationValues, 3, 0, index, errorId);
        const double pitch = component(value, "rotDetector",
            detectorRotationValues, 3, 1, index, errorId);
        const double roll = component(value, "rotDetector",
            detectorRotationValues, 3, 2, index, errorId);
        const DoubleVector originOffset{
            component(value, "offOrigin", offOrigin, 3, 0, index, errorId),
            component(value, "offOrigin", offOrigin, 3, 1, index, errorId),
            component(value, "offOrigin", offOrigin, 3, 2, index, errorId)};
        const double detectorOffsetU = component(value, "offDetector",
            offDetector, 2, 0, index, errorId);
        const double detectorOffsetV = component(value, "offDetector",
            offDetector, 2, 1, index, errorId);

        // The legacy CUDA parallel interpolated projector applies detector
        // offsets once before and once while forming its detector origin.
        const double detectorOffsetScale =
            geometry.parallel && cudaParallelInterpolatedOffsets ? 2.0 : 1.0;
        const DoubleVector detectorLocal = add(
            detectorRotation({0.0, detectorU0, detectorV0},
                             roll, pitch, yaw),
            {-(currentDsd - currentDso),
             detectorOffsetScale * detectorOffsetU,
             detectorOffsetScale * detectorOffsetV});
        const DoubleVector detectorU = scanRotation(
            detectorRotation({0.0, geometry.du, 0.0}, roll, pitch, yaw),
            alpha, theta, psi);
        const DoubleVector detectorV = scanRotation(
            detectorRotation({0.0, 0.0, -geometry.dv}, roll, pitch, yaw),
            alpha, theta, psi);
        DoubleVector sourceLocal{currentDso, 0.0, 0.0};
        if (geometry.parallel) {
            sourceLocal = add(
                detectorRotation({0.0, detectorU0, detectorV0},
                                 roll, pitch, yaw),
                {currentDso, 0.0, 0.0});
        }
        const double currentCor =
            viewValue(value, "COR", cor, index, errorId);
        const DoubleVector corOffset{
            -currentCor * std::sin(alpha),
            currentCor * std::cos(alpha),
            0.0};
        const DoubleVector source = add(
            subtract(scanRotation(sourceLocal, alpha, theta, psi),
                     originOffset),
            corOffset);
        const DoubleVector detector = add(
            subtract(scanRotation(detectorLocal, alpha, theta, psi),
                     originOffset),
            corOffset);
        const DoubleVector normal =
            normalize(cross(detectorU, detectorV), errorId);

        const double cudaAlpha = -alpha;
        const double cudaTheta = -theta;
        const double cudaPsi = -psi;
        const double detectorDistance = currentDsd - currentDso;
        auto backprojectionPoint = [&](DoubleVector point) {
            point = scanRotationTranspose(
                point, cudaAlpha, cudaTheta, cudaPsi);
            point.y -= detectorOffsetU;
            point.z -= detectorOffsetV;
            point.x += detectorDistance;
            point = detectorRotationTranspose(point, roll, pitch, yaw);
            point.x -= detectorDistance;
            point.y /= geometry.du;
            point.z /= geometry.dv;
            return point;
        };
        const DoubleVector backprojectionOrigin = backprojectionPoint({
            -0.5 * geometry.sx + 0.5 * geometry.dx + originOffset.x,
            -0.5 * geometry.sy + 0.5 * geometry.dy + originOffset.y,
            -0.5 * geometry.sz + 0.5 * geometry.dz + originOffset.z});
        const DoubleVector backprojectionNextX = backprojectionPoint({
            -0.5 * geometry.sx + 1.5 * geometry.dx + originOffset.x,
            -0.5 * geometry.sy + 0.5 * geometry.dy + originOffset.y,
            -0.5 * geometry.sz + 0.5 * geometry.dz + originOffset.z});
        const DoubleVector backprojectionNextY = backprojectionPoint({
            -0.5 * geometry.sx + 0.5 * geometry.dx + originOffset.x,
            -0.5 * geometry.sy + 1.5 * geometry.dy + originOffset.y,
            -0.5 * geometry.sz + 0.5 * geometry.dz + originOffset.z});
        const DoubleVector backprojectionNextZ = backprojectionPoint({
            -0.5 * geometry.sx + 0.5 * geometry.dx + originOffset.x,
            -0.5 * geometry.sy + 0.5 * geometry.dy + originOffset.y,
            -0.5 * geometry.sz + 1.5 * geometry.dz + originOffset.z});
        DoubleVector backprojectionSource{
            geometry.parallel ? 0.0 : currentDsd,
            -detectorOffsetU,
            -detectorOffsetV};
        backprojectionSource = detectorRotationTranspose(
            backprojectionSource, roll, pitch, yaw);
        backprojectionSource.x -= detectorDistance;
        backprojectionSource.y /= geometry.du;
        backprojectionSource.z /= geometry.dv;

        geometry.views.push_back({
            packed(source),
            packed(detector),
            packed(detectorU),
            packed(detectorV),
            packed(normal),
            static_cast<float>(currentDso),
            static_cast<float>(currentDsd),
            geometry.parallel ? 1u : 0u,
            0.0f,
            packed(backprojectionOrigin),
            packed(subtract(backprojectionNextX, backprojectionOrigin)),
            packed(subtract(backprojectionNextY, backprojectionOrigin)),
            packed(subtract(backprojectionNextZ, backprojectionOrigin)),
            packed(backprojectionSource),
            packed(originOffset, static_cast<float>(currentCor)),
            {static_cast<float>(detectorOffsetU),
             static_cast<float>(detectorOffsetV),
             static_cast<float>(std::sin(cudaAlpha)),
             static_cast<float>(std::cos(cudaAlpha))}});
    }
    return geometry;
}

}  // namespace tigre_geometry

#endif
