#ifndef AUTO_DIFF_TYPES_HPP_
#define AUTO_DIFF_TYPES_HPP_

#include <cstddef>

#include "auto_diff/common.hpp"

namespace halo_pm {
namespace ad {

// =============== Types ===============
struct NoneType {};

template <class T, size_t Len>
struct Vec {
  T data_[Len];
};

using Vec3f = Vec<float, 3>;


template <class T, size_t R, size_t C>
struct Mat {
  T data_[R * C];
};

using Mat3x3f = Mat<float, 3, 3>;


}  // namespace ad
}  // namespace halo_pm

#endif  // AUTO_DIFF_TYPES_HPP_
