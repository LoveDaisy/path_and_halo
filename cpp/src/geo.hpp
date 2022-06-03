#ifndef GEO_HPP_
#define GEO_HPP_

#include <cstddef>

namespace halo_pm {

constexpr float kPi = 3.14159265358f;
constexpr float kDegree2Rad = kPi / 180.0f;
constexpr float kRad2Degree = 180.0f / kPi;

void Ll2Xyz(const float* ll, float* xyz,                           // input & output, ll in degree
            size_t num = 1,                                        // data number
            size_t ll_step_bytes = 0, size_t xyz_step_bytes = 0);  // step of input & output


void Llr2Mat(const float* llr, float* mat,                           // input & output, llr in degree
             size_t num = 1,                                         // data number
             size_t llr_step_bytes = 0, size_t mat_step_bytes = 0);  // step of input & output

}  // namespace halo_pm

#endif  // GEO_HPP_
