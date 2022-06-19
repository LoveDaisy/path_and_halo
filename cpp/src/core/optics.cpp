#include "core/optics.hpp"

#include <cmath>
#include <tuple>

#include "auto_diff/ad.hpp"
#include "core/crystal.hpp"
#include "core/geo.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

namespace halo_pm {

Vec3f Refract(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1) {
  auto n = n0 / n1;
  auto c = ray_in.dot(norm);
  auto delta = n * n - (n * n - 1) / (c * c);

  if (delta <= 0) {
    return (Mat3x3f::Identity() - 2 * norm * norm.transpose()) * ray_in;
  } else {
    return n * ray_in + (std::sqrt(delta) - n) * c * norm;
  }
}


std::tuple<Vec3f, Mat3x3f> RefractAndDiff(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1) {
  auto n = n0 / n1;
  ad::Var vx{ ray_in.x() };
  ad::Var vy{ ray_in.y() };
  ad::Var vz{ ray_in.z() };

  auto c = vx * norm.x() + vy * norm.y() + vz * norm.z();
  auto delta = n * n - (n * n - 1.0f) / (c * c);
  auto delta_val = ad::Eval(delta);

  // Total reflection
  if (delta_val <= 0) {
    Mat3x3f m = Mat3x3f::Identity() - 2 * norm * norm.transpose();
    return std::make_tuple(m * ray_in, m);
  }

  // Real refraction
  else {
    auto a = (sqrt(delta) - n) * c;
    auto ox = n * vx + a * norm.x();
    auto oy = n * vy + a * norm.y();
    auto oz = n * vz + a * norm.z();

    using ad::Diff;
    using ad::Eval;
    using ad::wrt;

    Vec3f out_ray{ Eval(ox), Eval(oy), Eval(oz) };
    Mat3x3f jac{ { Diff(ox, wrt(vx)), Diff(ox, wrt(vy)), Diff(ox, wrt(vz)) },
                 { Diff(oy, wrt(vx)), Diff(oy, wrt(vy)), Diff(oy, wrt(vz)) },
                 { Diff(oz, wrt(vx)), Diff(ox, wrt(vy)), Diff(oz, wrt(vz)) } };
    Mat3x3f jac_norm = Mat3x3f::Identity() - ray_in * ray_in.transpose();  // Assume |ray_in| = 1.0
    return std::make_tuple(out_ray, jac * jac_norm);
  }
}


Vec3f TraceDirection(const Crystal& crystal,                       // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,        // May be input variables
                     const std::vector<int>& raypath, float wl) {  // Other parameter
  if (raypath.empty()) {
    return Vec3f{};
  }

  auto xyz = Ll2Xyz(ray_ll);
  LOG_DEBUG("ll2xyz: %s", ObjLogFormatter<Vec3f>{ xyz }.Format());

  // Only ONE face: pure reflection
  if (raypath.size() == 1) {
    const auto* n_ptr = crystal.face_norm_.get() + (raypath[0] - 1) * 3;
    Eigen::Map<const Vec3f> n(n_ptr);
    auto m = Mat3x3f::Identity() - 2 * n * n.transpose();
    return rot * m * rot.inverse() * xyz;
  }

  // Multiple faces: reflection + refraction
  else {
    Mat3x3f m = Eigen::Matrix3f::Identity();
    for (size_t i = 1; i + 1 < raypath.size(); i++) {
      Eigen::Map<const Vec3f> n(crystal.face_norm_.get() + (raypath[i] - 1) * 3);
      m = (Mat3x3f::Identity() - 2 * n * n.transpose()) * m;
    }

    auto refractive_index = GetIceRefractiveIndex(wl);
    Eigen::Map<const Vec3f> norm1(crystal.face_norm_.get() + (raypath.front() - 1) * 3);
    Eigen::Map<const Vec3f> norm2(crystal.face_norm_.get() + (raypath.back() - 1) * 3);
    auto r = rot.inverse() * xyz;
    r = Refract(r, norm1, 1.0f, refractive_index);
    r = m * r;
    r = Refract(r, norm2, refractive_index, 1.0f);
    return rot * r;
  }
}

}  // namespace halo_pm
