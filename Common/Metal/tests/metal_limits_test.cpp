/*
 * Standalone test of the Metal backend size limits.
 *
 * CODE by       Mike Croucher (MathWorks)
 */

#include "../tigre_metal_limits.hpp"
#include <cassert>
#include <iostream>

using namespace tigre_metal_limits;

template <typename Function>
void rejects(Function action) {
    bool rejected = false;
    try {
        action();
    } catch (const LimitExceeded&) {
        rejected = true;
    }
    assert(rejected);
}

int main() {
    const std::size_t uintMaximum = std::numeric_limits<std::uint32_t>::max();
    const std::size_t sizeMaximum = std::numeric_limits<std::size_t>::max();
    const std::size_t intMaximum = std::numeric_limits<std::int32_t>::max();
    assert(elements({intMaximum, 1, 1}) == intMaximum);
    assert(elements({65535, 65537}) == uintMaximum);
    rejects([] { elements({65536, 65536}); });
    rejects([] { elements({2048, 2048, 1024}); });
    rejects([] { elements({std::size_t(1) << 31, 1}); });
    rejects([] { elements({0, 16}); });
    assert(bytes(1024, 4, 4096) == 4096);
    rejects([] { bytes(1025, 4, 4096); });
    rejects([] { bytes(1024, 4, 4095); });
    rejects([] { bytes(1, 4, 0); });
    assert(checkedProduct(sizeMaximum, 1, sizeMaximum) == sizeMaximum);
    rejects([=] { bytes(sizeMaximum, 4, sizeMaximum); });
    rejects([=] { checkedProduct(uintMaximum, uintMaximum, uintMaximum); });
    std::cout << "Metal size boundary tests passed.\n";
}
