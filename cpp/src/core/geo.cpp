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
    auto c1 = std::cos(p(0) * kDegree2Rad);
    auto s1 = std::sin(p(0) * kDegree2Rad);
    auto c2 = std::cos(p(1) * kDegree2Rad);
    auto s2 = std::sin(p(1) * kDegree2Rad);
    auto c3 = std::cos(p(2) * kDegree2Rad);
    auto s3 = std::sin(p(2) * kDegree2Rad);

    auto& q = *reinterpret_cast<Mat3x3f*>(reinterpret_cast<uint8_t*>(mat) + i * mat_step_bytes);
    q(0, 0) = -s1 * c3 - c1 * s2 * s3;
    q(0, 1) = s1 * s3 - c1 * s2 * c3;
    q(0, 2) = c1 * c2;
    q(1, 0) = c1 * c3 - s1 * s2 * s3;
    q(1, 1) = -c1 * s3 - s1 * s2 * c3;
    q(1, 2) = s1 * c2;
    q(2, 0) = c2 * s3;
    q(2, 1) = c2 * c3;
    q(2, 2) = s2;
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
    tmp.x() = quat.w() * p0(0) + quat.y() * p0(2) - quat.z() * p0(1);
    tmp.y() = quat.w() * p0(1) - quat.x() * p0(2) + quat.z() * p0(0);
    tmp.z() = quat.w() * p0(2) + quat.x() * p0(1) - quat.y() * p0(0);
    tmp.w() = -quat.x() * p0(0) - quat.y() * p0(1) - quat.z() * p0(2);

    LOG_DEBUG("tmp=[% .4f,% .4f,% .4f,% .4f]", tmp.x(), tmp.y(), tmp.z(), tmp.w());

    auto& p1 = *reinterpret_cast<Vec3f*>(reinterpret_cast<uint8_t*>(xyz1) + i * xyz1_step_bytes);
    p1(0) = -tmp.w() * quat.x() + tmp.x() * quat.w() - tmp.y() * quat.z() + tmp.z() * quat.y();
    p1(1) = -tmp.w() * quat.y() + tmp.x() * quat.z() + tmp.y() * quat.w() - tmp.z() * quat.x();
    p1(2) = -tmp.w() * quat.z() - tmp.x() * quat.y() + tmp.y() * quat.x() + tmp.z() * quat.w();
    LOG_DEBUG("xyz1=[% .4f,% .4f,% .4f]", p1(0), p1(1), p1(2));
  }
}
}  // namespace halo_pm