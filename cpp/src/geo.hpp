#ifndef GEO_HPP_
#define GEO_HPP_

#include <cstddef>

namespace halo_pm {

/**
 * @brief Convert a longitude-latitude representation into a x-y-z representation.
 *
 * @param ll  [input]  (longitude, latitude) data, in degree
 * @param xyz [output] (x, y, z) data
 * @param num
 * @param ll_step_bytes  Step bytes to next data. If leave it as default 0, the function will regard input data as
 *                       contiguous and regard step = 2 * sizeof(float)
 * @param xyz_step_bytes Step bytes to next data. Similar to above.
 */
void Ll2Xyz(const float* ll, float* xyz,                           // input & output, ll in degree
            size_t num = 1,                                        // data number
            size_t ll_step_bytes = 0, size_t xyz_step_bytes = 0);  // step of input & output


void Xyz2Ll(const float* xyz, float* ll,                           // input & output, ll in degree, xyz normalized
            size_t num = 1,                                        // data number
            size_t xyz_step_bytes = 0, size_t ll_step_bytes = 0);  // step of input & output


void Llr2Mat(const float* llr, float* mat,                           // input & output, llr in degree
             size_t num = 1,                                         // data number
             size_t llr_step_bytes = 0, size_t mat_step_bytes = 0);  // step of input & output


/**
 * @brief Rotate a point by a given quaternion.
 *
 * @param quat [input]  The given quaternion. **ONLY ONE** data, 4 floats.
 * @param xyz0 [input]  Points to be rotated.
 * @param xyz1 [output] Points rotated. It can be the same as xyz0, which indicates rotate points inplace.
 * @param num Data number for xyz data.
 * @param xyz0_step_bytes
 * @param xyz1_step_bytes
 */
void RotateByQuat(const float* quat, const float* xyz0, float* xyz1,        // input & output
                  size_t num = 1,                                           // data number
                  size_t xyz0_step_bytes = 0, size_t xyz1_step_bytes = 0);  // steps

}  // namespace halo_pm

#endif  // GEO_HPP_
