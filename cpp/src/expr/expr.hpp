#ifndef AUTO_DIFF_EXPR_HPP_
#define AUTO_DIFF_EXPR_HPP_

#include <cstddef>
#include <tuple>
#include <utility>

#include "expr/op.hpp"

namespace halo_pm {
namespace internal {

// =============== Expressions ===============
template <class V>
struct Var {
  V val_;
  size_t id_;

  static size_t global_id_;

  Var(V v) : val_(v), id_(global_id_++) {}
  Var(const Var& other) : val_(other.val_), id_(other.id_) {}
  Var(Var&& other) : val_(other.val_), id_(other.id_) {}
  ~Var() = default;

  Var& operator=(const Var& other) {
    if (this != &other) {
      val_ = other.val_;
      id_ = other.id_;
    }
    return *this;
  }

  Var& operator=(Var&& other) {
    if (this != &other) {
      val_ = other.val_;
      id_ = other.id_;
    }
    return *this;
  }
};

template <class VV>
Var(VV&& v) -> Var<VV>;

template <class V>
size_t Var<V>::global_id_ = 0;


template <class Op, class V>
struct UnaryExpr {
  V val_;

  UnaryExpr(V v) : val_(v) {}
  UnaryExpr(const UnaryExpr& other) : val_(other.val_) {}
  UnaryExpr(UnaryExpr&& other) : val_(other.val_) {}
  ~UnaryExpr() = default;

  UnaryExpr& operator=(const UnaryExpr& other) {
    if (this != &other) {
      val_ = other.val_;
    }
    return *this;
  }

  UnaryExpr& operator=(UnaryExpr&& other) {
    if (this != &other) {
      val_ = std::move(other.val_);
    }
    return *this;
  }
};


template <class Op, class L, class R>
struct BinaryExpr {
  L l_;
  R r_;

  BinaryExpr(L l, R r) : l_(l), r_(r) {}
  BinaryExpr(const BinaryExpr& other) : l_(other.l_), r_(other.r_) {}
  BinaryExpr(BinaryExpr&& other) : l_(other.l_), r_(other.r_) {}
  ~BinaryExpr() = default;

  BinaryExpr& operator=(const BinaryExpr& other) {
    if (this != &other) {
      l_ = other.l_;
      r_ = other.r_;
    }
    return *this;
  }

  BinaryExpr& operator=(BinaryExpr&& other) {
    if (this != &other) {
      l_ = std::move(other.l_);
      r_ = std::move(other.r_);
    }
    return *this;
  }
};


// =============== Expression alias ===============
template <class T>
using NegExpr = UnaryExpr<NegativeOp, T>;

template <class L, class R>
using AddExpr = BinaryExpr<AddOp, L, R>;

template <class L, class R>
using MinusExpr = BinaryExpr<MinusOp, L, R>;

template <class L, class R>
using TimesExpr = BinaryExpr<TimesOp, L, R>;

template <class L, class R>
using DivideExpr = BinaryExpr<DivideOp, L, R>;

template <class T>
using SqrtExpr = UnaryExpr<SqrtOp, T>;


// =============== Operator overloads ===============
template <class T>
NegExpr<T> operator-(T v) {
  return UnaryExpr<NegativeOp, T>{ v };
}

template <class L, class R>
AddExpr<L, R> operator+(L l, R r) {
  return AddExpr<L, R>{ l, r };
}

template <class L, class R>
MinusExpr<L, R> operator-(L l, R r) {
  return MinusExpr<L, R>{ l, r };
}

template <class L, class R>
TimesExpr<L, R> operator*(L l, R r) {
  return TimesExpr<L, R>{ l, r };
}

template <class L, class R>
DivideExpr<L, R> operator/(L l, R r) {
  return DivideExpr<L, R>{ l, r };
}

// =============== Function overloads ===============
template <class T>
constexpr SqrtExpr<T> sqrt(T val) {
  return SqrtExpr<T>{ val };
}

}  // namespace internal
}  // namespace halo_pm

#endif  // AUTO_DIFF_EXPR_HPP_
