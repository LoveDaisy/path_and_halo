#include "geo.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include "math.hpp"

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
}  // namespace halo_pm