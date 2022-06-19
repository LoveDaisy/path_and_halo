#ifndef CORE_GEO_H_
#define CORE_GEO_H_

#include <cstddef>
#include <tuple>

#include "auto_diff/ad.hpp"
#include "core/types.hpp"

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
void Ll2Xyz(const Vec2f* ll, Vec3f* xyz,                           // input & output, ll in degree
            size_t num = 1,                                        // data number
            size_t ll_step_bytes = 0, size_t xyz_step_bytes = 0);  // step of input & output

// A convenience overload for ONE data
Vec3f Ll2Xyz(const Vec2f& ll);


void Xyz2Ll(const Vec3f* xyz, Vec2f* ll,                           // input & output, ll in degree, xyz normalized
            size_t num = 1,                                        // data number
            size_t xyz_step_bytes = 0, size_t ll_step_bytes = 0);  // step of input & output


/**
 * @brief Convert axis-angle representation of a rotation into matrix representation.
 *
 * @param llr [input]  (longitude, latitude, roll) axis-angle representation, in degree.
 * @param mat [output] matrix representation. Row major.
 * @param num
 * @param llr_step_bytes
 * @param mat_step_bytes
 */
void Llr2Mat(const Vec3f* llr, Mat3x3f* mat,                         // input & output, llr in degree
             size_t num = 1,                                         // data number
             size_t llr_step_bytes = 0, size_t mat_step_bytes = 0);  // step of input & output


/**
 * @brief Rotate a point by a given quaternion.
 *
 * @param quat [input]  The given quaternion. **ONLY ONE** data, 4 floats [xi, yj, zk, w]
 * @param xyz0 [input]  Points to be rotated.
 * @param xyz1 [output] Points rotated. It can be the same as xyz0, which indicates rotate points inplace.
 * @param num Data number for xyz data.
 * @param xyz0_step_bytes
 * @param xyz1_step_bytes
 */
void RotateByQuat(const Quatf& quat, const Vec3f* xyz0, Vec3f* xyz1,        // input & output
                  size_t num = 1,                                           // data number
                  size_t xyz0_step_bytes = 0, size_t xyz1_step_bytes = 0);  // steps


template <class QW, class QX, class QY, class QZ,  // quaternion
          class VX, class VY, class VZ>            // input vector
auto                                               // output vector
QuatRotExpr(const ad::Var<QW>& qw, const ad::Var<QX>& qx, const ad::Var<QY>& qy, const ad::Var<QZ>& qz,  // quaternion
            const ad::Var<VX>& vx, const ad::Var<VY>& vy, const ad::Var<VZ>& vz) {                       // input vector
  auto tmp_qw = -qx * vx - qy * vy - qz * vz;
  auto tmp_qx = qw * vx + qy * vz - qz * vy;
  auto tmp_qy = qw * vy - qx * vz + qz * vx;
  auto tmp_qz = qw * vz + qx * vy - qy * vx;

  auto ox = -tmp_qw * qx + tmp_qx * qw - tmp_qy * qz + tmp_qz * qy;
  auto oy = -tmp_qw * qy + tmp_qx * qz + tmp_qy * qw - tmp_qz * qx;
  auto oz = -tmp_qw * qz - tmp_qx * qy + tmp_qy * qx + tmp_qz * qw;

  return std::make_tuple(ox, oy, oz);
}


Mat3x3f VecNormalizeDiff(const Vec3f& v);

}  // namespace halo_pm

#endif  // CORE_GEO_H_
