#include "geo.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include "math.hpp"
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


void Llr2Mat(const float* llr, float* mat,                    // input & output, llr in degree
             size_t num,                                      // data number
             size_t llr_step_bytes, size_t mat_step_bytes) {  // step of input & output
  assert(llr_step_bytes % sizeof(float) == 0);
  assert(mat_step_bytes % sizeof(float) == 0);

  const float* p = llr;
  float* q = mat;
  for (size_t i = 0; i < num; i++) {
    auto c1 = std::cos(p[0] * kDegree2Rad);
    auto s1 = std::sin(p[0] * kDegree2Rad);
    auto c2 = std::cos(p[1] * kDegree2Rad);
    auto s2 = std::sin(p[1] * kDegree2Rad);
    auto c3 = std::cos(p[2] * kDegree2Rad);
    auto s3 = std::sin(p[2] * kDegree2Rad);

    q[0] = -s1 * c3 - c1 * s2 * s3;
    q[1] = s1 * s3 - c1 * s2 * c3;
    q[2] = c1 * c2;
    q[3] = c1 * c3 - s1 * s2 * s3;
    q[4] = -c1 * s3 - s1 * s2 * c3;
    q[5] = s1 * c2;
    q[6] = c2 * s3;
    q[7] = c2 * c3;
    q[8] = s2;

    p += llr_step_bytes == 0 ? 3 : llr_step_bytes / sizeof(float);
    q += mat_step_bytes == 0 ? 9 : mat_step_bytes / sizeof(float);
  }
}


void RotateByQuat(const float* quat, const float* xyz0, float* xyz1,  // input & output
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

  float tmp[4]{};

  for (size_t i = 0; i < num; i++) {
    const float* p0 = xyz0 + i * (xyz0_step_bytes / sizeof(float));
    tmp[0] = -quat[1] * p0[0] - quat[2] * p0[1] - quat[3] * p0[2];
    tmp[1] = quat[0] * p0[0] + quat[2] * p0[2] - quat[3] * p0[1];
    tmp[2] = quat[0] * p0[1] - quat[1] * p0[2] + quat[3] * p0[0];
    tmp[3] = quat[0] * p0[2] + quat[1] * p0[1] - quat[2] * p0[0];

    LOG_DEBUG("tmp=[% .4f,% .4f,% .4f,% .4f]", tmp[0], tmp[1], tmp[2], tmp[3]);

    float* p1 = xyz1 + i * (xyz1_step_bytes / sizeof(float));
    p1[0] = -tmp[0] * quat[1] + tmp[1] * quat[0] - tmp[2] * quat[3] + tmp[3] * quat[2];
    p1[1] = -tmp[0] * quat[2] + tmp[1] * quat[3] + tmp[2] * quat[0] - tmp[3] * quat[1];
    p1[2] = -tmp[0] * quat[3] - tmp[1] * quat[2] + tmp[2] * quat[1] + tmp[3] * quat[0];
    LOG_DEBUG("xyz1=[% .4f,% .4f,% .4f]", p1[0], p1[1], p1[2]);
  }
}
}  // namespace halo_pm