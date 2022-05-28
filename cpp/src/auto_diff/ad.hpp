#ifndef AUTO_DIFF_AD_HPP_
#define AUTO_DIFF_AD_HPP_

#include <cstddef>
#include <tuple>

#include "auto_diff/common.hpp"
#include "auto_diff/expr.hpp"

namespace halo_pm {

namespace ad {
///////////////////////////////////////////////////////////////////////////////////////////////////
// Make auto differentiation easy and straight. Like following codes:
// ~~~c++
// VarExpr f(VarExpr x, VarExpr y, VarExpr z) { ... }
//
// VarExpr x = 0.1f;
// VarExpr y = 0.0f;
// VarExpr z{-0.3};
// VarExpr u = f(x, y, z);
//
// Vecf u_val = f(x, y, z);       // implicit conversion: VarExpr --> Vecf
// Vecf u_val = Evaluate(u);      // or explicit evaluation
// Mat3x3f u_jac = Jacobian(u);   // Jacobian() returns a MatType
// ~~~

// =============== Evaluate ===============
constexpr float Evaluate(float v) {
  return v;
}

constexpr float Evaluate(ZeroExpr /* e */) {
  return 0;
}

constexpr float Evaluate(OneExpr /* e */) {
  return 1;
}

template <class... T>
constexpr auto Evaluate(VecExpr<T...> e) -> Vec<float, sizeof...(T)> {
  Vec<float, sizeof...(T)> res;
  for (size_t i = 0; i < sizeof...(T); i++) {
    res.data_[i] = Evaluate(std::get<i>(e.val_));
  }
  return res;
}

template <class T>
constexpr auto Evaluate(VarExpr<T> e) {
  return Evaluate(e.val_);
}

template <class L, class R>
constexpr auto Evaluate(AddExpr<L, R> e) -> float {
  return Evaluate(e.l_) + Evaluate(e.r_);
}

}  // namespace ad

}  // namespace halo_pm

#endif  // AUTO_DIFF_AD_HPP_
