#ifndef AUTO_DIFF_AD_HPP_
#define AUTO_DIFF_AD_HPP_

#include <cstddef>

#include "auto_diff/expr.hpp"
#include "auto_diff/traits.hpp"
#include "auto_diff/types.hpp"

namespace halo_pm {

namespace ad {
///////////////////////////////////////////////////////////////////////////////////////////////////
// Make auto differentiation easy and straight. Like following codes:
// ~~~c++
// Var f(Var x, Var y, Var z) { ... }
//
// Var x = 0.1;     // scalar
// Var y{};         // zero-initialize
// Var z{-0.3};
// Var u = f(x, y, z);
//
// Vec u_val = f(x, y, z);    // implicit conversion: Var --> Vec, Var is an expression and Vec is real values.
// Vec u_val = u.eval();      // or explicit evaluation
// Mat u_jac = u.jacobian();  // jacobian() returns a MatType
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

template <class T>
constexpr float Evaluate(ScalarExpr<T> e) {
  return e.val_;
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
