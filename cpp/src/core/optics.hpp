#ifndef CORE_OPTICS_H_
#define CORE_OPTICS_H_

#include <tuple>

#include "core/types.hpp"

namespace halo_pm {

constexpr float kMinWavelength = 350.0f;
constexpr float kMaxWavelength = 850.0f;
constexpr float kDefaultWavelength = 546.1f;  // e-line


Vec3f Refract(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1);

std::tuple<Vec3f, Mat3x3f> RefractAndDiff(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1);


// Forward declaration
struct Crystal;

Vec3f TraceDirection(const Crystal& crystal,                                           // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,                            // May be input variables
                     const std::vector<int>& raypath, float wl = kDefaultWavelength);  // Other parameter


}  // namespace halo_pm

#endif  // CORE_OPTICS_H_
