#include "optics/optics.hpp"

#include <cmath>
#include <limits>

#include "core/types.hpp"
#include "optics/crystal.hpp"
#include "util/log.hpp"

namespace halo_pm {

Vec3f Refract(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1) {
  auto n = n0 / n1;
  auto c = ray_in.dot(norm);
  auto delta = n * n - (n * n - 1) / (c * c);

  if (delta <= 0) {
    Vec3f r{};
    r.fill(std::numeric_limits<float>::quiet_NaN());
    return r;
  } else {
    return n * ray_in + (std::sqrt(delta) - n) * c * norm;
  }
}


Vec3f TraceDirection(const Crystal& crystal,                       // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,        // May be input variables
                     const std::vector<int>& raypath, float wl) {  // Other parameter
  if (raypath.empty()) {
    return Vec3f{};
  }

  auto xyz = Ll2Xyz(ray_ll);

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


std::tuple<Vec3f, Mat3x4f>                                     // (output vector, Jacobian wrt quaternion)
TraceDirDiffQuat(const Crystal& crystal,                       // crystal
                 const Quatf& rot, const Vec2f& ray_ll,        // input quaternion, input vector
                 const std::vector<int>& raypath, float wl) {  // other parameter
  if (raypath.empty()) {
    return std::make_tuple(Vec3f{}, Mat3x4f{});
  }

  using ad::Diff;
  using ad::Eval;
  using ad::Var;
  using ad::wrt;

  auto xyz = Ll2Xyz(ray_ll);

  Var qw = rot.w();
  Var qx = rot.x();
  Var qy = rot.y();
  Var qz = rot.z();

  Var vx = xyz.x();
  Var vy = xyz.y();
  Var vz = xyz.z();

  // Only ONE face: pure reflection
  if (raypath.size() == 1) {
    const auto* n_ptr = crystal.face_norm_.get() + (raypath[0] - 1) * 3;
    Eigen::Map<const Vec3f> n(n_ptr);
    auto m = Mat3x3f::Identity() - 2 * n * n.transpose();

    auto [v1x, v1y, v1z] = QuatRotExpr(qw, -qx, -qy, -qz, vx, vy, vz);
    auto v2x = v1x * m(0, 0) + v1y * m(0, 1) + v1z * m(0, 2);
    auto v2y = v1x * m(1, 0) + v1y * m(1, 1) + v1z * m(1, 2);
    auto v2z = v1x * m(2, 0) + v1y * m(2, 1) + v1z * m(2, 2);
    auto [v3x, v3y, v3z] = QuatRotExpr(qw, qx, qy, qz, v2x, v2y, v2z);

    Vec3f new_xyz{ Eval(v3x), Eval(v3y), Eval(v3z) };
    Mat3x4f jac{ { Diff(v3x, wrt(qw)), Diff(v3x, wrt(qx)), Diff(v3x, wrt(qy)), Diff(v3x, wrt(qz)) },
                 { Diff(v3y, wrt(qw)), Diff(v3y, wrt(qx)), Diff(v3y, wrt(qy)), Diff(v3y, wrt(qz)) },
                 { Diff(v3z, wrt(qw)), Diff(v3z, wrt(qx)), Diff(v3z, wrt(qy)), Diff(v3z, wrt(qz)) } };
    return std::make_tuple(new_xyz, jac);
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

    Vec3f new_xyz;
    Mat3x4f jac;

    Vec3f v1;
    Mat3x4f j1;
    {
      auto [v1x, v1y, v1z] = QuatRotExpr(qw, -qx, -qy, -qz, vx, vy, vz);
      v1 = Vec3f{ Eval(v1x), Eval(v1y), Eval(v1z) };
      j1 = Mat3x4f{ { Diff(v1x, wrt(qw)), Diff(v1x, wrt(qx)), Diff(v1x, wrt(qy)), Diff(v1x, wrt(qz)) },
                    { Diff(v1y, wrt(qw)), Diff(v1y, wrt(qx)), Diff(v1y, wrt(qy)), Diff(v1y, wrt(qz)) },
                    { Diff(v1z, wrt(qw)), Diff(v1z, wrt(qx)), Diff(v1z, wrt(qy)), Diff(v1z, wrt(qz)) } };
    }

    Vec3f v2;
    Mat3x3f j2;
    {
      Var vx = v1.x();
      Var vy = v1.y();
      Var vz = v1.z();
      auto [delta, rx, ry, rz, lx, ly, lz] = RefractExpr(vx, vy, vz, norm1, 1.0f, refractive_index);
      v2 = Vec3f{ Eval(rx), Eval(ry), Eval(rz) };
      j2 = Mat3x3f{ { Diff(rx, wrt(vx)), Diff(rx, wrt(vy)), Diff(rx, wrt(vz)) },
                    { Diff(ry, wrt(vx)), Diff(ry, wrt(vy)), Diff(ry, wrt(vz)) },
                    { Diff(rz, wrt(vx)), Diff(rz, wrt(vy)), Diff(rz, wrt(vz)) } };
    }

    Vec3f v3 = m * v2;
    Mat3x3f j3 = m;

    Vec3f v4;
    Mat3x3f j4;
    {
      Var vx = v3.x();
      Var vy = v3.y();
      Var vz = v3.z();
      auto [delta, rx, ry, rz, lx, ly, lz] = RefractExpr(vx, vy, vz, norm2, refractive_index, 1.0f);
      v4 = Vec3f{ Eval(rx), Eval(ry), Eval(rz) };
      j4 = Mat3x3f{ { Diff(rx, wrt(vx)), Diff(rx, wrt(vy)), Diff(rx, wrt(vz)) },
                    { Diff(ry, wrt(vx)), Diff(ry, wrt(vy)), Diff(ry, wrt(vz)) },
                    { Diff(rz, wrt(vx)), Diff(rz, wrt(vy)), Diff(rz, wrt(vz)) } };
    }

    {
      Var vx = v4.x();
      Var vy = v4.y();
      Var vz = v4.z();
      auto [ox, oy, oz] = QuatRotExpr(qw, qx, qy, qz, vx, vy, vz);
      Mat3x4f j5_q{ { Diff(ox, wrt(qw)), Diff(ox, wrt(qx)), Diff(ox, wrt(qy)), Diff(ox, wrt(qz)) },
                    { Diff(oy, wrt(qw)), Diff(oy, wrt(qx)), Diff(oy, wrt(qy)), Diff(oy, wrt(qz)) },
                    { Diff(oz, wrt(qw)), Diff(oz, wrt(qx)), Diff(oz, wrt(qy)), Diff(oz, wrt(qz)) } };
      Mat3x3f j5_v{ { Diff(ox, wrt(vx)), Diff(ox, wrt(vy)), Diff(ox, wrt(vz)) },
                    { Diff(oy, wrt(vx)), Diff(oy, wrt(vy)), Diff(oy, wrt(vz)) },
                    { Diff(oz, wrt(vx)), Diff(oz, wrt(vy)), Diff(oz, wrt(vz)) } };

      new_xyz = Vec3f{ Eval(ox), Eval(oy), Eval(oz) };
      jac = j5_q + j5_v * j4 * j3 * j2 * j1;
    }

    return std::make_tuple(new_xyz, jac);
  }
}

}  // namespace halo_pm
