#ifndef CORE_MATH_H_
#define CORE_MATH_H_

#include <cmath>
#include <cstddef>

#include "core/types.hpp"

namespace halo_pm {

constexpr float kPi = 3.14159265358f;
constexpr float kDegree2Rad = kPi / 180.0f;
constexpr float kRad2Degree = 180.0f / kPi;

}  // namespace halo_pm

#endif  // CORE_MATH_H_
