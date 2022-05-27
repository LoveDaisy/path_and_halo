#ifndef GEO_HPP_
#define GEO_HPP_

#include <cstddef>

namespace halo_pm {

void Ll2Xyz(const float* ll, float* xyz,                           // input & output
            size_t num = 1,                                        // data number
            size_t ll_step_bytes = 0, size_t xyz_step_bytes = 0);  // step of input & output


void Llr2Mat(const float* llr, float* mat,                           // input & output
             size_t num = 1,                                         // data number
             size_t llr_step_bytes = 0, size_t mat_step_bytes = 0);  // step of input & output

}  // namespace halo_pm

#endif  // GEO_HPP_
