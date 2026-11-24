/*
 * Metal compute kernels for the TIGRE Metal backend.
 *
 * CODE by       Mike Croucher (MathWorks)
 */

#include <metal_stdlib>
using namespace metal;

struct Params {
    uint nx;
    uint ny;
    uint nz;
    uint nu;
    uint nv;
    uint angleCount;
    uint mode;
    uint variant;
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

struct ViewGeometry {
    float4 source;
    float4 detector;
    float4 detectorU;
    float4 detectorV;
    float4 detectorNormal;
    float dso;
    float dsd;
    uint parallel;
    float padding;
    float4 backprojectionOrigin;
    float4 backprojectionDeltaX;
    float4 backprojectionDeltaY;
    float4 backprojectionDeltaZ;
    float4 backprojectionSource;
    float4 originOffset;
    float4 detectorParameters;
};

struct Ray {
    float3 source;
    float3 detector;
};

inline ulong image_index(constant Params& p, uint x, uint y, uint z) {
    return ulong(x) + ulong(y) * p.nx + ulong(z) * p.nx * p.ny;
}

inline Ray make_ray(constant Params& p, device const ViewGeometry& view,
                    uint u, uint v) {
    v = p.nv - v - 1;
    Ray ray;
    ray.detector = view.detector.xyz + float(u) * view.detectorU.xyz +
                   float(v) * view.detectorV.xyz;
    ray.source = p.mode == 1
        ? view.source.xyz + float(u) * view.detectorU.xyz +
              float(v) * view.detectorV.xyz
        : view.source.xyz;
    return ray;
}

inline bool intersect_axis(float source, float direction,
                           float lower, float upper,
                           thread float& entry, thread float& exit) {
    if (abs(direction) < 1.0e-20f) {
        return source >= lower && source <= upper;
    }
    float first = (lower - source) / direction;
    float second = (upper - source) / direction;
    if (first > second) {
        const float temporary = first;
        first = second;
        second = temporary;
    }
    entry = max(entry, first);
    exit = min(exit, second);
    return entry <= exit;
}

inline bool intersect_volume(constant Params& p, Ray ray,
                             thread float& entry, thread float& exit) {
    const float3 direction = ray.detector - ray.source;
    entry = 0.0f;
    exit = 1.0f;
    return intersect_axis(ray.source.x, direction.x, -p.sx * 0.5f,
                          p.sx * 0.5f, entry, exit) &&
           intersect_axis(ray.source.y, direction.y, -p.sy * 0.5f,
                          p.sy * 0.5f, entry, exit) &&
           intersect_axis(ray.source.z, direction.z, -p.sz * 0.5f,
                          p.sz * 0.5f, entry, exit) &&
           entry <= exit;
}

inline float voxel(device const float* image, constant Params& p,
                   int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 ||
        x >= int(p.nx) || y >= int(p.ny) || z >= int(p.nz)) {
        return 0.0f;
    }
    return image[image_index(p, uint(x), uint(y), uint(z))];
}

inline float trilinear(device const float* image, constant Params& p,
                       float3 point) {
    const float x = (point.x + p.sx * 0.5f) / p.dx - 0.5f;
    const float y = (point.y + p.sy * 0.5f) / p.dy - 0.5f;
    const float z = (point.z + p.sz * 0.5f) / p.dz - 0.5f;
    const int x0 = int(floor(x));
    const int y0 = int(floor(y));
    const int z0 = int(floor(z));
    const float fx = x - float(x0);
    const float fy = y - float(y0);
    const float fz = z - float(z0);
    const float c00 = mix(voxel(image, p, x0, y0, z0),
                          voxel(image, p, x0 + 1, y0, z0), fx);
    const float c10 = mix(voxel(image, p, x0, y0 + 1, z0),
                          voxel(image, p, x0 + 1, y0 + 1, z0), fx);
    const float c01 = mix(voxel(image, p, x0, y0, z0 + 1),
                          voxel(image, p, x0 + 1, y0, z0 + 1), fx);
    const float c11 = mix(voxel(image, p, x0, y0 + 1, z0 + 1),
                          voxel(image, p, x0 + 1, y0 + 1, z0 + 1), fx);
    return mix(mix(c00, c10, fy), mix(c01, c11, fy), fz);
}

inline float interpolated_ray(device const float* image, constant Params& p,
                              Ray ray) {
    float entry;
    float exit;
    if (!intersect_volume(p, ray, entry, exit)) {
        return 0.0f;
    }
    const float3 direction = ray.detector - ray.source;
    const float3 indexDirection =
        float3(direction.x / p.dx, direction.y / p.dy, direction.z / p.dz);
    const uint sampleCount =
        max(1u, uint(ceil(length(indexDirection) / p.accuracy)));
    const uint first = uint(max(0.0f, floor(entry * float(sampleCount))));
    const uint last =
        min(sampleCount, uint(ceil(exit * float(sampleCount))));
    float sum = 0.0f;
    for (uint sample = first; sample <= last; ++sample) {
        const float parameter = float(sample) / float(sampleCount);
        sum += trilinear(image, p, ray.source + direction * parameter);
    }
    return sum * length(direction) / float(sampleCount);
}

inline int voxel_index(float coordinate, float size, float spacing,
                       uint count) {
    return clamp(int(floor((coordinate + size * 0.5f) / spacing)),
                 0, int(count) - 1);
}

inline float next_boundary_parameter(float source, float direction,
                                     float size, float spacing, int index) {
    if (abs(direction) < 1.0e-20f) {
        return INFINITY;
    }
    const float boundary = direction > 0.0f
        ? -size * 0.5f + float(index + 1) * spacing
        : -size * 0.5f + float(index) * spacing;
    return (boundary - source) / direction;
}

inline float siddon_ray(device const float* image, constant Params& p,
                        Ray ray) {
    float entry;
    float exit;
    if (!intersect_volume(p, ray, entry, exit) || exit <= entry) {
        return 0.0f;
    }
    const float3 direction = ray.detector - ray.source;
    const float rayLength = length(direction);
    const float3 point = ray.source + direction * (entry + 1.0e-6f);
    int ix = voxel_index(point.x, p.sx, p.dx, p.nx);
    int iy = voxel_index(point.y, p.sy, p.dy, p.ny);
    int iz = voxel_index(point.z, p.sz, p.dz, p.nz);
    const uint strideY = p.nx;
    const uint strideZ = p.nx * p.ny;
    uint linearIndex =
        uint(ix) + uint(iy) * strideY + uint(iz) * strideZ;
    const int stepX = direction.x >= 0.0f ? 1 : -1;
    const int stepY = direction.y >= 0.0f ? 1 : -1;
    const int stepZ = direction.z >= 0.0f ? 1 : -1;
    const float deltaX =
        abs(direction.x) < 1.0e-20f ? INFINITY : p.dx / abs(direction.x);
    const float deltaY =
        abs(direction.y) < 1.0e-20f ? INFINITY : p.dy / abs(direction.y);
    const float deltaZ =
        abs(direction.z) < 1.0e-20f ? INFINITY : p.dz / abs(direction.z);
    float nextX =
        next_boundary_parameter(ray.source.x, direction.x, p.sx, p.dx, ix);
    float nextY =
        next_boundary_parameter(ray.source.y, direction.y, p.sy, p.dy, iy);
    float nextZ =
        next_boundary_parameter(ray.source.z, direction.z, p.sz, p.dz, iz);
    float current = entry;
    float sum = 0.0f;
    const uint maximumSteps = p.nx + p.ny + p.nz + 3;
    for (uint step = 0; step < maximumSteps && current < exit; ++step) {
        if (ix < 0 || iy < 0 || iz < 0 ||
            ix >= int(p.nx) || iy >= int(p.ny) || iz >= int(p.nz)) {
            break;
        }
        const float next = min(exit, min(nextX, min(nextY, nextZ)));
        sum += image[linearIndex] *
               max(0.0f, next - current) * rayLength;
        current = next;
        const float tolerance =
            8.0f * FLT_EPSILON * max(1.0f, abs(current));
        if (nextX <= current + tolerance) {
            ix += stepX;
            linearIndex += stepX;
            nextX += deltaX;
        }
        if (nextY <= current + tolerance) {
            iy += stepY;
            linearIndex += stepY * strideY;
            nextY += deltaY;
        }
        if (nextZ <= current + tolerance) {
            iz += stepZ;
            linearIndex += stepZ * strideZ;
            nextZ += deltaZ;
        }
    }
    return sum;
}

kernel void tigre_ax(device const float* image [[buffer(0)]],
                     device const ViewGeometry* views [[buffer(1)]],
                     device float* output [[buffer(2)]],
                     constant Params& p [[buffer(3)]],
                     uint index [[thread_position_in_grid]]) {
    const uint detectorPixels = p.nu * p.nv;
    if (index >= detectorPixels * p.angleCount) {
        return;
    }
    const uint angleIndex = index / detectorPixels;
    const uint within = index - angleIndex * detectorPixels;
    const uint u = within / p.nv;
    const uint v = within - u * p.nv;
    const Ray ray = make_ray(p, views[angleIndex], u, v);
    output[index] =
        p.variant == 1 ? siddon_ray(image, p, ray)
                       : interpolated_ray(image, p, ray);
}

inline float projection_fetch(device const float* projections,
                              constant Params& p, uint angle,
                              int u, int v) {
    if (u < 0 || v < 0 || u >= int(p.nu) || v >= int(p.nv)) {
        return 0.0f;
    }
    return projections[ulong(v) + ulong(u) * p.nv +
                       ulong(angle) * p.nu * p.nv];
}

inline float projection_sample(device const float* projections,
                               constant Params& p, uint angle,
                               float u, float v) {
    const int u0 = int(floor(u));
    const int v0 = int(floor(v));
    const float fu = u - float(u0);
    const float fv = v - float(v0);
    const float lower =
        mix(projection_fetch(projections, p, angle, u0, v0),
            projection_fetch(projections, p, angle, u0 + 1, v0), fu);
    const float upper =
        mix(projection_fetch(projections, p, angle, u0, v0 + 1),
            projection_fetch(projections, p, angle, u0 + 1, v0 + 1), fu);
    return mix(lower, upper, fv);
}

inline bool detector_coordinates(
        constant Params& p, device const ViewGeometry& view,
        float3 voxelPosition, thread float& u, thread float& v,
        thread float& sourceDistanceSquared, thread float& rayLength,
        thread float& axialDistance) {
    const float3 source = view.source.xyz;
    const float3 detector = view.detector.xyz;
    const float3 detectorU = view.detectorU.xyz;
    const float3 detectorV = view.detectorV.xyz;
    const float3 normal = view.detectorNormal.xyz;
    const float3 planeDirection = detector - source;
    float3 detectorPoint;
    if (p.mode == 1) {
        const float denominator = dot(planeDirection, normal);
        if (abs(denominator) <= 1.0e-8f) {
            return false;
        }
        const float scale =
            dot(detector - voxelPosition, normal) / denominator;
        detectorPoint = voxelPosition + planeDirection * scale;
        axialDistance = 1.0f;
        sourceDistanceSquared = 1.0f;
        rayLength = length(planeDirection);
    } else {
        const float3 sourceToVoxel = voxelPosition - source;
        const float denominator = dot(sourceToVoxel, normal);
        const float detectorDistance = dot(planeDirection, normal);
        if (abs(denominator) <= 1.0e-8f ||
            detectorDistance / denominator <= 0.0f) {
            return false;
        }
        detectorPoint =
            source + sourceToVoxel * (detectorDistance / denominator);
        sourceDistanceSquared = dot(sourceToVoxel, sourceToVoxel);
        if (sourceDistanceSquared <= 1.0e-8f) {
            return false;
        }
        rayLength = length(detectorPoint - source);
        axialDistance = denominator;
        if (axialDistance <= 1.0e-8f) {
            return false;
        }
    }
    const float3 relative = detectorPoint - detector;
    const float uLengthSquared = dot(detectorU, detectorU);
    const float vLengthSquared = dot(detectorV, detectorV);
    if (uLengthSquared <= 1.0e-8f || vLengthSquared <= 1.0e-8f) {
        return false;
    }
    u = dot(relative, detectorU) / uLengthSquared;
    v = dot(relative, detectorV) / vLengthSquared;
    return true;
}

kernel void tigre_atb(device const float* projections [[buffer(0)]],
                      device const ViewGeometry* views [[buffer(1)]],
                      device float* output [[buffer(2)]],
                      constant Params& p [[buffer(3)]],
                      uint index [[thread_position_in_grid]]) {
    const uint voxelCount = p.nx * p.ny * p.nz;
    if (index >= voxelCount) {
        return;
    }
    const uint plane = p.nx * p.ny;
    const uint zIndex = index / plane;
    const uint within = index - zIndex * plane;
    const uint yIndex = within / p.nx;
    const uint xIndex = within - yIndex * p.nx;
    const float x = -p.sx * 0.5f + p.dx * 0.5f + float(xIndex) * p.dx;
    const float y = -p.sy * 0.5f + p.dy * 0.5f + float(yIndex) * p.dy;
    const float z = -p.sz * 0.5f + p.dz * 0.5f + float(zIndex) * p.dz;
    float sum = 0.0f;
    for (uint angleIndex = 0; angleIndex < p.angleCount; ++angleIndex) {
        device const ViewGeometry& view = views[angleIndex];
        // Match the detector-scaled coordinates and legacy weights used by
        // TIGRE's CUDA backprojectors.
        float3 point =
            view.backprojectionOrigin.xyz +
            float(xIndex) * view.backprojectionDeltaX.xyz +
            float(yIndex) * view.backprojectionDeltaY.xyz +
            float(zIndex) * view.backprojectionDeltaZ.xyz;
        point.y -= view.originOffset.w / p.du;

        float sample = 0.0f;
        if (p.mode == 1) {
            const float u = point.y + float(p.nu) * 0.5f - 0.5f;
            const float v = point.z + float(p.nv) * 0.5f - 0.5f;
            sample = projection_sample(projections, p, angleIndex, u, v);
        } else {
            const float3 direction =
                point - view.backprojectionSource.xyz;
            if (abs(direction.x) <= 1.0e-8f) {
                continue;
            }
            const float t =
                (view.dso - view.dsd - view.backprojectionSource.x) /
                direction.x;
            const float u =
                direction.y * t + view.backprojectionSource.y +
                float(p.nu) * 0.5f;
            const float v =
                direction.z * t + view.backprojectionSource.z +
                float(p.nv) * 0.5f;
            sample = projection_sample(
                projections, p, angleIndex, u - 0.5f, v - 0.5f);

            const float realX = x + view.originOffset.x;
            const float realY = y + view.originOffset.y;
            const float realZ = z + view.originOffset.z;
            const float sinAlpha = view.detectorParameters.z;
            const float cosAlpha = view.detectorParameters.w;
            if (p.variant == 1) {
                const float realSourceX = view.dso * cosAlpha;
                const float realSourceY = -view.dso * sinAlpha;
                const float detectorX = -(view.dsd - view.dso);
                const float detectorY =
                    (-p.su + p.du) * 0.5f + u * p.du +
                    view.detectorParameters.x;
                const float realDetectorX =
                    detectorX * cosAlpha + detectorY * sinAlpha;
                const float realDetectorY =
                    -detectorX * sinAlpha + detectorY * cosAlpha;
                const float realDetectorZ =
                    (-p.sv + p.dv) * 0.5f + v * p.dv +
                    view.detectorParameters.y;
                const float deltaSourceX = realSourceX - realDetectorX;
                const float deltaSourceY = realSourceY - realDetectorY;
                const float rayLength = sqrt(
                    deltaSourceX * deltaSourceX +
                    deltaSourceY * deltaSourceY +
                    realDetectorZ * realDetectorZ);
                const float voxelSourceX = realSourceX - realX;
                const float voxelSourceY = realSourceY - realY;
                const float sourceDistanceSquared =
                    voxelSourceX * voxelSourceX +
                    voxelSourceY * voxelSourceY + realZ * realZ;
                sample *= rayLength * rayLength * rayLength /
                          (view.dsd * sourceDistanceSquared);
            } else {
                const float ratio =
                    (view.dso +
                     (realY + view.originOffset.w) * sinAlpha -
                     realX * cosAlpha) /
                    view.dso;
                sample /= ratio * ratio;
            }
        }
        sum += sample;
    }
    if (p.mode == 0 && p.variant == 1) {
        sum *= p.dx * p.dy * p.dz / (p.du * p.dv);
    }
    output[index] = sum;
}

inline float3 backward_gradient(device const float* image,
                                constant Params& p,
                                int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 ||
        x >= int(p.nx) || y >= int(p.ny) || z >= int(p.nz)) {
        return float3(0.0f);
    }
    const ulong index = image_index(p, uint(x), uint(y), uint(z));
    const float value = image[index];
    float3 gradient = float3(0.0f);
    if (z > 0) {
        gradient.x = value - image[index - ulong(p.nx) * p.ny];
    }
    if (y > 0) {
        gradient.y = value - image[index - p.nx];
    }
    if (x > 0) {
        gradient.z = value - image[index - 1];
    }
    return gradient;
}

inline float3 tv_weights(float3 gradient, constant Params& p) {
    if (p.variant == 0) {
        return float3(1.0f);
    }
    const float denominator =
        max(p.alpha * p.alpha, 1.0e-8f);
    return exp(-(gradient * gradient) / denominator);
}

inline float weighted_norm(float3 gradient, float3 weight) {
    return sqrt(dot(weight, gradient * gradient)) + 1.0e-8f;
}

kernel void tigre_tv_gradient(device const float* image [[buffer(0)]],
                              device float* gradient [[buffer(1)]],
                              constant Params& p [[buffer(3)]],
                              uint index [[thread_position_in_grid]]) {
    const uint count = p.nx * p.ny * p.nz;
    if (index >= count) {
        return;
    }
    const uint plane = p.nx * p.ny;
    const uint z = index / plane;
    const uint within = index - z * plane;
    const uint y = within / p.nx;
    const uint x = within - y * p.nx;
    const float3 current =
        backward_gradient(image, p, int(x), int(y), int(z));
    const float3 nextX =
        backward_gradient(image, p, int(x + 1), int(y), int(z));
    const float3 nextY =
        backward_gradient(image, p, int(x), int(y + 1), int(z));
    const float3 nextZ =
        backward_gradient(image, p, int(x), int(y), int(z + 1));
    const float3 weight = tv_weights(current, p);
    const float3 weightX = tv_weights(nextX, p);
    const float3 weightY = tv_weights(nextY, p);
    const float3 weightZ = tv_weights(nextZ, p);
    gradient[index] =
        dot(weight, current) / weighted_norm(current, weight) -
        weightX.z * nextX.z / weighted_norm(nextX, weightX) -
        weightY.y * nextY.y / weighted_norm(nextY, weightY) -
        weightZ.x * nextZ.x / weighted_norm(nextZ, weightZ);
}

kernel void tigre_difference(device const float* first [[buffer(0)]],
                             device const float* second [[buffer(1)]],
                             device float* output [[buffer(2)]],
                             constant uint& count [[buffer(3)]],
                             uint index [[thread_position_in_grid]]) {
    if (index < count) {
        output[index] = first[index] - second[index];
    }
}

kernel void tigre_combine(device float* gradient [[buffer(0)]],
                          device const float* priorGradient [[buffer(1)]],
                          constant float& ratio [[buffer(2)]],
                          constant uint& count [[buffer(3)]],
                          uint index [[thread_position_in_grid]]) {
    if (index < count) {
        gradient[index] =
            (1.0f - ratio) * gradient[index] +
            ratio * priorGradient[index];
    }
}

kernel void tigre_scale_subtract(device float* image [[buffer(0)]],
                                 device const float* gradient [[buffer(1)]],
                                 constant float& scale [[buffer(2)]],
                                 constant uint& count [[buffer(3)]],
                                 uint index [[thread_position_in_grid]]) {
    if (index < count) {
        image[index] -= scale * gradient[index];
    }
}

kernel void tigre_reduce_squares(device const float* input [[buffer(0)]],
                                 device float* output [[buffer(1)]],
                                 constant uint& count [[buffer(2)]],
                                 uint gridIndex [[thread_position_in_grid]],
                                 uint localIndex [[thread_index_in_threadgroup]],
                                 uint groupIndex [[threadgroup_position_in_grid]],
                                 uint groupSize [[threads_per_threadgroup]]) {
    threadgroup float shared[256];
    shared[localIndex] =
        gridIndex < count ? input[gridIndex] * input[gridIndex] : 0.0f;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = groupSize / 2; stride > 0; stride /= 2) {
        if (localIndex < stride) {
            shared[localIndex] += shared[localIndex + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (localIndex == 0) {
        output[groupIndex] = shared[0];
    }
}

kernel void tigre_reduce_sum(device const float* input [[buffer(0)]],
                             device float* output [[buffer(1)]],
                             constant uint& count [[buffer(2)]],
                             uint gridIndex [[thread_position_in_grid]],
                             uint localIndex [[thread_index_in_threadgroup]],
                             uint groupIndex [[threadgroup_position_in_grid]],
                             uint groupSize [[threads_per_threadgroup]]) {
    threadgroup float shared[256];
    shared[localIndex] = gridIndex < count ? input[gridIndex] : 0.0f;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = groupSize / 2; stride > 0; stride /= 2) {
        if (localIndex < stride) {
            shared[localIndex] += shared[localIndex + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (localIndex == 0) {
        output[groupIndex] = shared[0];
    }
}

inline float divergence(device const float* pz,
                        device const float* py,
                        device const float* px,
                        constant Params& p,
                        uint x, uint y, uint z) {
    const ulong index = image_index(p, x, y, z);
    float value = pz[index] + py[index] + px[index];
    if (z > 0) {
        value -= pz[index - ulong(p.nx) * p.ny];
    }
    if (y > 0) {
        value -= py[index - p.nx];
    }
    if (x > 0) {
        value -= px[index - 1];
    }
    return value;
}

kernel void tigre_tv_update_u(device const float* source [[buffer(0)]],
                              device const float* pz [[buffer(1)]],
                              device const float* py [[buffer(2)]],
                              device const float* px [[buffer(3)]],
                              device float* image [[buffer(4)]],
                              constant Params& p [[buffer(5)]],
                              uint index [[thread_position_in_grid]]) {
    const uint count = p.nx * p.ny * p.nz;
    if (index >= count) {
        return;
    }
    const uint plane = p.nx * p.ny;
    const uint z = index / plane;
    const uint within = index - z * plane;
    const uint y = within / p.nx;
    const uint x = within - y * p.nx;
    const float div = divergence(pz, py, px, p, x, y, z);
    image[index] = image[index] * (1.0f - p.tau1) +
                   p.tau1 * (source[index] + div / p.lambda);
}

kernel void tigre_tv_update_p(device const float* image [[buffer(0)]],
                              device float* pz [[buffer(1)]],
                              device float* py [[buffer(2)]],
                              device float* px [[buffer(3)]],
                              constant Params& p [[buffer(4)]],
                              uint index [[thread_position_in_grid]]) {
    const uint count = p.nx * p.ny * p.nz;
    if (index >= count) {
        return;
    }
    const uint plane = p.nx * p.ny;
    const uint z = index / plane;
    const uint within = index - z * plane;
    const uint y = within / p.nx;
    const uint x = within - y * p.nx;
    const float value = image[index];
    const float gz = z + 1 < p.nz ? image[index + plane] - value : 0.0f;
    const float gy = y + 1 < p.ny ? image[index + p.nx] - value : 0.0f;
    const float gx = x + 1 < p.nx ? image[index + 1] - value : 0.0f;
    const float qz = pz[index] + p.tau2 * gz;
    const float qy = py[index] + p.tau2 * gy;
    const float qx = px[index] + p.tau2 * gx;
    const float norm = max(1.0f, sqrt(qz * qz + qy * qy + qx * qx));
    pz[index] = qz / norm;
    py[index] = qy / norm;
    px[index] = qx / norm;
}
