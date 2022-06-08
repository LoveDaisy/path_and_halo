#ifndef AUTO_DIFF_AD_HPP_
#define AUTO_DIFF_AD_HPP_

#include <cstddef>

#include "auto_diff/expr.hpp"
#include "auto_diff/op.hpp"
#include "auto_diff/traits.hpp"

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
// Vecf u_val = f(x, y, z);      // implicit conversion: VarExpr --> Vecf
// Vecf u_val = Evaluate(u);     // or explicit evaluation
// auto ua = Differentiate(u, wrt(x, y));  // It is an Expr
// float ua_val = Evaluate(ua);            // explicitly convert to float
// ~~~

// =============== Evaluate ===============
constexpr float Evaluate(float v) {
  return v;
}

template <class T>
constexpr auto Evaluate(VarExpr<T> e) {
  return Evaluate(e.val_);
}

template <class L, class R>
constexpr auto Evaluate(AddExpr<L, R> e) {
  return Evaluate(e.l_) + Evaluate(e.r_);
}

template <class L, class R>
constexpr auto Evaluate(MinusExpr<L, R> e) {
  return Evaluate(e.l_) - Evaluate(e.r_);
}

template <class L, class R>
constexpr auto Evaluate(TimesExpr<L, R> e) {
  return Evaluate(e.l_) * Evaluate(e.r_);
}

template <class L, class R>
constexpr auto Evaluate(DivideExpr<L, R> e) {
  return Evaluate(e.l_) / Evaluate(e.r_);
}


// =============== Differentiate ===============
template <class V>
struct wrt {
  V val_;

  wrt(V&& v) : val_(std::move(v)) {}
  wrt(const V& v) : val_(v) {}
};


template <class Y, class X>
constexpr auto Differentiate(VarExpr<Y> v, wrt<X> x) {
  // Direct case, i.e. v.val_ is data type (non-expr type)
  if constexpr (!traits::is_expr_v<Y>) {
    return 0.0f;
  }

  // Expr
  else {
    return Differentiate(v.val_, x);
  }
};


template <class V>
constexpr auto Differentiate(VarExpr<V> v, wrt<VarExpr<V>> x) {
  // Direct case, i.e. v.val_ is data type (non-expr type)
  if constexpr (!traits::is_expr_v<V>) {
    if (v.id_ == x.val_.id_) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  }

  // Expr
  else {
    return Differentiate(v.val_, x);
  }
}


template <class L, class R, class X>
constexpr auto Differentiate(AddExpr<L, R> v, wrt<X> x) {
  return Differentiate(v.l_, x) + Differentiate(v.r_, x);
}

template <class L, class R, class X>
constexpr auto Differentiate(MinusExpr<L, R> v, wrt<X> x) {
  return Differentiate(v.l_, x) - Differentiate(v.r_, x);
}

template <class L, class R, class X>
constexpr auto Differentiate(TimesExpr<L, R> v, wrt<X> x) {
  return Differentiate(v.l_, x) * v.r_ + Differentiate(v.r_, x) * v.l_;
}

template <class L, class R, class X>
constexpr auto Differentiate(DivideExpr<L, R> v, wrt<X> x) {
  return (Differentiate(v.l_, x) * v.r_ - Differentiate(v.r_, x) * v.l_) / v.r_ / v.r_;
}

}  // namespace ad

}  // namespace halo_pm

#endif  // AUTO_DIFF_AD_HPP_
