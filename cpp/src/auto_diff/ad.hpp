#ifndef AUTO_DIFF_AD_HPP_
#define AUTO_DIFF_AD_HPP_

#include <cstddef>
#include <tuple>
#include <utility>

#include "auto_diff/common.hpp"
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
// Vecf u_val = f(x, y, z);    // implicit conversion: VarExpr --> Vecf
// Vecf u_val = Evaluate(u);   // or explicit evaluation
// auto jac = Differentiate(u, wrt(x, y));   // Differentiate() returns a MatExpr (Jacobian matrix)
// auto jac_val = Evaluate(jac);
// ~~~

// =============== Evaluate ===============
constexpr float Evaluate(float v) {
  return v;
}

template <class T, size_t C, size_t R>
constexpr Mat<T, C, R> Evaluate(Mat<T, R, C> v) {
  return v;
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

template <class L, class R>
constexpr auto Evaluate(MinusExpr<L, R> e) -> float {
  return Evaluate(e.l_) - Evaluate(e.r_);
}

template <class L, class R>
constexpr auto Evaluate(TimesExpr<L, R> e) -> float {
  return Evaluate(e.l_) * Evaluate(e.r_);
}

template <class L, class R>
constexpr auto Evaluate(DivideExpr<L, R> e) -> float {
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
    return Mat<float, traits::var_dim<Y>::dim, traits::var_dim<X>::dim>{};
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
      constexpr size_t n = traits::var_dim<V>::dim;
      Mat<float, n, n> res{};
      for (size_t i = 0; i < n; i++) {
        res.data_[i * n + i] = 1.0f;
      }
      return res;
    } else {
      return Mat<float, traits::var_dim<V>::dim, traits::var_dim<V>::dim>{};
    }
  }

  // Expr
  else {
    return Differentiate(v.val_, x);
  }
}

template <class T, size_t R, size_t C, class X>
constexpr auto Differentiate(const Mat<T, R, C>& m, wrt<X> /* x */) -> Mat<T, R, C> {
  return m;
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
