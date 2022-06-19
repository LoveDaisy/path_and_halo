#ifndef CORE_OPTICS_H_
#define CORE_OPTICS_H_

#include <tuple>
#include <vector>

#include "auto_diff/ad.hpp"
#include "core/geo.hpp"
#include "core/types.hpp"

namespace halo_pm {

constexpr float kMinWavelength = 350.0f;
constexpr float kMaxWavelength = 850.0f;
constexpr float kDefaultWavelength = 546.1f;  // e-line


Vec3f Refract(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1);

template <class VX, class VY, class VZ>
auto RefractExpr(const VX& vx, const VY& vy, const VZ& vz, const Vec3f& norm, float n0, float n1) {
  auto [vnx, vny, vnz] = NormalizeExpr(vx, vy, vz);

  auto n = n0 / n1;
  auto c = vnx * norm.x() + vny * norm.y() + vnz * norm.z();
  auto delta = n * n - (n * n - 1) / (c * c);

  auto m = Mat3x3f::Identity() - 2 * norm * norm.transpose();
  auto lx = vnx * m(0, 0) + vny * m(0, 1) + vnz * m(0, 2);
  auto ly = vnx * m(1, 0) + vny * m(1, 1) + vnz * m(1, 2);
  auto lz = vnx * m(2, 0) + vny * m(2, 1) + vnz * m(2, 2);

  auto a = (sqrt(delta) - n) * c;
  auto rx = n * vnx + a * norm.x();
  auto ry = n * vny + a * norm.y();
  auto rz = n * vnz + a * norm.z();
  return std::make_tuple(delta, rx, ry, rz, lx, ly, lz);
}


// Forward declaration
struct Crystal;

Vec3f TraceDirection(const Crystal& crystal,                                           // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,                            // May be input variables
                     const std::vector<int>& raypath, float wl = kDefaultWavelength);  // Other parameter


std::tuple<Vec3f, Mat3x4f>                               // (output vector, Jacobian wrt quaternion)
TraceDirDiffQuat(const Crystal& crystal,                 // crystal
                 const Quatf& rot, const Vec2f& ray_ll,  // input quaternion, input vector
                 const std::vector<int>& raypath, float wl = kDefaultWavelength);  // other parameter


}  // namespace halo_pm

#endif  // CORE_OPTICS_H_
