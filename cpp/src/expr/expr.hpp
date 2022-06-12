#ifndef AUTO_DIFF_EXPR_HPP_
#define AUTO_DIFF_EXPR_HPP_

#include <cstddef>
#include <tuple>

#include "expr/op.hpp"

namespace halo_pm {
namespace internal {

// =============== Expressions ===============
template <class V>
struct VarExpr {
  V val_;
  size_t id_;

  static size_t global_id_;

  VarExpr(V&& v) : val_(std::forward<V>(v)), id_(global_id_++) {}
  VarExpr(const V& v) : val_(v), id_(global_id_++) {}
};

template <class V>
size_t VarExpr<V>::global_id_ = 0;


template <class Op, class V>
struct UnaryExpr {
  Op op_;
  V val_;

  UnaryExpr(V&& v) : op_{}, val_{ std::forward<V>(v) } {}
};


template <class Op, class L, class R>
struct BinaryExpr {
  Op op_;
  L l_;
  R r_;

  BinaryExpr(L&& l, R&& r) : op_{}, l_{ std::forward<L>(l) }, r_{ std::forward<R>(r) } {}
};


// =============== Expression alias ===============
template <class L, class R>
using AddExpr = BinaryExpr<AddOp, L, R>;

template <class L, class R>
using MinusExpr = BinaryExpr<MinusOp, L, R>;

template <class L, class R>
using TimesExpr = BinaryExpr<TimesOp, L, R>;

template <class L, class R>
using DivideExpr = BinaryExpr<DivideOp, L, R>;


// =============== Operator overloads ===============
template <class L, class R>
AddExpr<L, R> operator+(L&& l, R&& r) {
  return AddExpr<L, R>{ std::forward<L>(l), std::forward<R>(r) };
}

template <class L, class R>
MinusExpr<L, R> operator-(L&& l, R&& r) {
  return MinusExpr<L, R>{ std::forward<L>(l), std::forward<R>(r) };
}

template <class L, class R>
TimesExpr<L, R> operator*(L&& l, R&& r) {
  return TimesExpr<L, R>{ std::forward<L>(l), std::forward<R>(r) };
}

template <class L, class R>
DivideExpr<L, R> operator/(L&& l, R&& r) {
  return DivideExpr<L, R>{ std::forward<L>(l), std::forward<R>(r) };
}

}  // namespace internal
}  // namespace halo_pm

#endif  // AUTO_DIFF_EXPR_HPP_
