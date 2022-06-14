#ifndef CORE_TYPES_H_
#define CORE_TYPES_H_

#include <Eigen/Eigen>
#include <cstddef>

namespace halo_pm {

template <class T, size_t Dim>
using Vec = Eigen::Matrix<T, Dim, 1>;

template <class T, size_t R, size_t C>
using Mat = Eigen::Matrix<T, R, C>;

using Vec2f = Eigen::Vector2f;
using Vec3f = Eigen::Vector3f;
using Vec4f = Eigen::Vector4f;

using Mat3x3f = Eigen::Matrix3f;

using Quatf = Eigen::Quaternionf;

}  // namespace halo_pm

#endif  // CORE_TYPES_H_
