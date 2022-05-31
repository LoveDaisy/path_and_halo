#ifndef AUTO_DIFF_TYPES_HPP_
#define AUTO_DIFF_TYPES_HPP_

#include <cstddef>
#include <cstring>

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

template <class T>
struct Mat<T, 1, 1> {
  T data_[1];

  operator float() const { return data_[0]; }
};

using Mat3x3f = Mat<float, 3, 3>;

template <class T, size_t R, size_t C>
Mat<T, R, C> operator+(Mat<T, R, C> a, Mat<T, R, C> b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] + b.data_[i];
  }
  return res;
}

template <class T, size_t R, size_t C>
Mat<T, R, C> operator-(Mat<T, R, C> a, Mat<T, R, C> b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] - b.data_[i];
  }
  return res;
}

}  // namespace ad
}  // namespace halo_pm

#endif  // AUTO_DIFF_TYPES_HPP_
