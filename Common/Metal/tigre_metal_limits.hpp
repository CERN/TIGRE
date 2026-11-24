/*
 * Size limits for Metal buffers and indexing in the TIGRE Metal backend.
 *
 * CODE by       Mike Croucher (MathWorks)
 */

#ifndef TIGRE_METAL_LIMITS_HPP
#define TIGRE_METAL_LIMITS_HPP

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace tigre_metal_limits {

class LimitExceeded : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

inline std::size_t checkedProduct(std::size_t first, std::size_t second,
                                  std::size_t limit) {
    if (first == 0 || second == 0 || first > limit / second) {
        throw LimitExceeded(
            "Array size exceeds Metal's indexing or buffer limit. "
            "Use a smaller array or the CPU/CUDA backend.");
    }
    return first * second;
}

inline std::size_t elements(std::initializer_list<std::size_t> dimensions) {
    std::size_t count = 1;
    for (std::size_t dimension : dimensions) {
        // Voxel and detector coordinates are converted to signed int in shaders.
        if (dimension > static_cast<std::size_t>(
                std::numeric_limits<std::int32_t>::max())) {
            throw LimitExceeded("Each Metal array dimension must be below 2^31.");
        }
        count = checkedProduct(
            count, dimension, std::numeric_limits<std::uint32_t>::max());
    }
    return count;
}

inline std::size_t bytes(std::size_t count, std::size_t elementSize,
                         std::size_t maximumBufferLength) {
    return checkedProduct(count, elementSize, maximumBufferLength);
}

}  // namespace tigre_metal_limits
#endif
