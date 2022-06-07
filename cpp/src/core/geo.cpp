#include "core/geo.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include "core/math.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

namespace halo_pm {

void Ll2Xyz(const float* ll, float* xyz,                    // input & output, ll in degree
            size_t num,                                     // data number
            size_t ll_step_bytes, size_t xyz_step_bytes) {  // step of input & output
  assert(ll_step_bytes % sizeof(float) == 0);
  assert(xyz_step_bytes % sizeof(float) == 0);

  const float* p = ll;
  float* q = xyz;
  for (size_t i = 0; i < num; i++) {
    q[0] = std::cos(p[1] * kDegree2Rad) * std::cos(p[0] * kDegree2Rad);
    q[1] = std::cos(p[1] * kDegree2Rad) * std::sin(p[0] * kDegree2Rad);
    q[2] = std::sin(p[1] * kDegree2Rad);

    p += ll_step_bytes == 0 ? 2 : ll_step_bytes / sizeof(float);
    q += xyz_step_bytes == 0 ? 3 : xyz_step_bytes / sizeof(float);
  }
}


void Xyz2Ll(const float* xyz, float* ll,                    // input & output, ll in degree, xyz normalized
            size_t num,                                     // data number
            size_t xyz_step_bytes, size_t ll_step_bytes) {  // step of input & output
  assert(xyz_step_bytes % sizeof(float) == 0);
  assert(ll_step_bytes % sizeof(float) == 0);

  const float* p = xyz;
  float* q = ll;
  for (size_t i = 0; i < num; i++) {
    q[1] = std::asin(p[2]) * kRad2Degree;
    q[0] = std::atan2(p[1], p[0]) * kRad2Degree;

    p += xyz_step_bytes == 0 ? 3 : xyz_step_bytes / sizeof(float);
    q += ll_step_bytes == 0 ? 2 : ll_step_bytes / sizeof(float);
  }
}


void Llr2Mat(const Vec3f* llr, Mat3x3f* mat,                  // input & output, llr in degree
             size_t num,                                      // data number
             size_t llr_step_bytes, size_t mat_step_bytes) {  // step of input & output
  assert(llr_step_bytes % sizeof(float) == 0);
  assert(mat_step_bytes % sizeof(float) == 0);
  if (llr_step_bytes == 0) {
    llr_step_bytes = 3 * sizeof(float);
  }
  if (mat_step_bytes == 0) {
    mat_step_bytes = 9 * sizeof(float);
  }

  for (size_t i = 0; i < num; i++) {
    const auto& p = *reinterpret_cast<const Vec3f*>(reinterpret_cast<const uint8_t*>(llr) + i * llr_step_bytes);
    auto c1 = std::cos(p.x_ * kDegree2Rad);
    auto s1 = std::sin(p.x_ * kDegree2Rad);
    auto c2 = std::cos(p.y_ * kDegree2Rad);
    auto s2 = std::sin(p.y_ * kDegree2Rad);
    auto c3 = std::cos(p.z_ * kDegree2Rad);
    auto s3 = std::sin(p.z_ * kDegree2Rad);

    auto& q = *reinterpret_cast<Mat3x3f*>(reinterpret_cast<uint8_t*>(mat) + i * mat_step_bytes);
    q[0] = -s1 * c3 - c1 * s2 * s3;
    q[1] = s1 * s3 - c1 * s2 * c3;
    q[2] = c1 * c2;
    q[3] = c1 * c3 - s1 * s2 * s3;
    q[4] = -c1 * s3 - s1 * s2 * c3;
    q[5] = s1 * c2;
    q[6] = c2 * s3;
    q[7] = c2 * c3;
    q[8] = s2;
  }
}


void RotateByQuat(const Quatf& quat, const Vec3f* xyz0, Vec3f* xyz1,  // input & output
                  size_t num,                                         // data number
                  size_t xyz0_step_bytes, size_t xyz1_step_bytes) {   // steps
  // Quaternion rotation. See wiki link for details: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
  // q =    [cos(theta/2),  [ux, uy, uz] * sin(theta/2)]
  // p =    [0,             [px, py, pz]]
  // q^-1 = [cos(theta/2), -[ux, uy, uz] * sin(theta/2)]
  // p' =   q * p * q^-1

  assert(xyz0_step_bytes % sizeof(float) == 0);
  assert(xyz1_step_bytes % sizeof(float) == 0);
  if (xyz0_step_bytes == 0) {
    xyz0_step_bytes = 3 * sizeof(float);
  }
  if (xyz1_step_bytes == 0) {
    xyz1_step_bytes = 3 * sizeof(float);
  }

  Quatf tmp{};

  for (size_t i = 0; i < num; i++) {
    const auto& p0 = *reinterpret_cast<const Vec3f*>(reinterpret_cast<const uint8_t*>(xyz0) + i * xyz0_step_bytes);
    tmp.x_ = quat.w_ * p0.x_ + quat.y_ * p0.z_ - quat.z_ * p0.y_;
    tmp.y_ = quat.w_ * p0.y_ - quat.x_ * p0.z_ + quat.z_ * p0.x_;
    tmp.z_ = quat.w_ * p0.z_ + quat.x_ * p0.y_ - quat.y_ * p0.x_;
    tmp.w_ = -quat.x_ * p0.x_ - quat.y_ * p0.y_ - quat.z_ * p0.z_;

    LOG_DEBUG("tmp=[% .4f,% .4f,% .4f,% .4f]", tmp.x_, tmp.y_, tmp.z_, tmp.w_);

    auto& p1 = *reinterpret_cast<Vec3f*>(reinterpret_cast<uint8_t*>(xyz1) + i * xyz1_step_bytes);
    p1.x_ = -tmp.w_ * quat.x_ + tmp.x_ * quat.w_ - tmp.y_ * quat.z_ + tmp.z_ * quat.y_;
    p1.y_ = -tmp.w_ * quat.y_ + tmp.x_ * quat.z_ + tmp.y_ * quat.w_ - tmp.z_ * quat.x_;
    p1.z_ = -tmp.w_ * quat.z_ - tmp.x_ * quat.y_ + tmp.y_ * quat.x_ + tmp.z_ * quat.w_;
    LOG_DEBUG("xyz1=[% .4f,% .4f,% .4f]", p1.x_, p1.y_, p1.z_);
  }
}
}  // namespace halo_pm