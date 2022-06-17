#ifndef AUTO_DIFF_AD_HPP_
#define AUTO_DIFF_AD_HPP_

#include <cmath>
#include <cstddef>

#include "auto_diff/traits.hpp"
#include "expr/expr.hpp"
#include "expr/op.hpp"

namespace halo_pm {

namespace internal {
// =============== Eval ===============
constexpr float Eval(float v) {
  return v;
}

template <class T>
constexpr auto Eval(const VarExpr<T>& e) {
  return Eval(e.val_);
}

template <class T>
constexpr auto Eval(const NegExpr<T>& e) {
  return -1.0f * Eval(e.val_);
}

template <class L, class R>
constexpr auto Eval(const AddExpr<L, R>& e) {
  return Eval(e.l_) + Eval(e.r_);
}

template <class L, class R>
constexpr auto Eval(const MinusExpr<L, R>& e) {
  return Eval(e.l_) - Eval(e.r_);
}

template <class L, class R>
constexpr auto Eval(const TimesExpr<L, R>& e) {
  return Eval(e.l_) * Eval(e.r_);
}

template <class L, class R>
constexpr auto Eval(const DivideExpr<L, R>& e) {
  return Eval(e.l_) / Eval(e.r_);
}

template <class T>
constexpr auto Eval(const SqrtExpr<T>& e) {
  return std::sqrt(Eval(e.val_));
}

}  // namespace internal

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
// Vecf u_val = Eval(u);     // or explicit evaluation
// auto ua = Diff(u, wrt(x, y));  // It is an Expr
// float ua_val = Eval(ua);            // explicitly convert to float
// ~~~

using internal::Eval;
using internal::VarExpr;

}  // namespace ad


namespace internal {
// =============== Diff ===============
template <class V>
struct wrt {
  V val_;

  wrt(V&& v) : val_(std::move(v)) {}
  wrt(const V& v) : val_(v) {}
};


template <class X>
constexpr auto Diff(float /* v */, const wrt<X>& /* x */) {
  return 0.0f;
}


template <class Y, class X>
constexpr auto Diff(const VarExpr<Y>& v, const wrt<X>& x) -> float {
  // Direct case, i.e. v.val_ is data type (non-expr type)
  if constexpr (!ad::traits::is_expr_v<Y>) {
    return 0.0f;
  }

  // Expr
  else {
    return Eval(Diff(v.val_, x));
  }
};


template <class V>
constexpr auto Diff(const VarExpr<V>& v, const wrt<VarExpr<V>>& x) -> float {
  // Direct case, i.e. v.val_ is data type (non-expr type)
  if constexpr (!ad::traits::is_expr_v<V>) {
    if (v.id_ == x.val_.id_) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  }

  // Expr
  else {
    return Eval(Diff(v.val_, x));
  }
}


template <class T, class X>
constexpr auto Diff(const NegExpr<T>& v, const wrt<X>& x) {
  return Eval(-1.0f * Diff(v.val_, x));
}


template <class L, class R, class X>
constexpr auto Diff(const AddExpr<L, R>& v, const wrt<X>& x) {
  return Eval(Diff(v.l_, x) + Diff(v.r_, x));
}

template <class L, class R, class X>
constexpr auto Diff(const MinusExpr<L, R>& v, const wrt<X>& x) {
  return Eval(Diff(v.l_, x) - Diff(v.r_, x));
}

template <class L, class R, class X>
constexpr auto Diff(const TimesExpr<L, R>& v, const wrt<X>& x) {
  return Eval(Diff(v.l_, x) * v.r_ + v.l_ * Diff(v.r_, x));
}

template <class L, class R, class X>
constexpr auto Diff(const DivideExpr<L, R>& v, const wrt<X>& x) {
  return Eval(Diff(v.l_, x) / v.r_ - v.l_ * Diff(v.r_, x) / (v.r_ * v.r_));
}

template <class T, class X>
constexpr auto Diff(const SqrtExpr<T>& v, const wrt<X>& x) {
  return Eval(0.5f / sqrt(v.val_) * Diff(v.val_, x));
}

}  // namespace internal

namespace ad {

using internal::Diff;
using internal::wrt;

}  // namespace ad

}  // namespace halo_pm

#endif  // AUTO_DIFF_AD_HPP_
