#ifndef CORE_MATH_H_
#define CORE_MATH_H_

#include <cmath>
#include <cstddef>

#include "core/types.hpp"

namespace halo_pm {

constexpr float kPi = 3.14159265358f;
constexpr float kDegree2Rad = kPi / 180.0f;
constexpr float kRad2Degree = 180.0f / kPi;

namespace linalgo {

template <class T, size_t N>
T Norm2(const Vec<T, N>& x) {
  double n = 0;
  for (size_t i = 0; i < N; i++) {
    n += x.data_[i] * x.data_[i];
  }
  return n;
}


template <class T, size_t N>
T Norm(const Vec<T, N>& x) {
  return std::sqrt(Norm2(x));
}

}  // namespace linalgo

}  // namespace halo_pm

#endif  // CORE_MATH_H_
