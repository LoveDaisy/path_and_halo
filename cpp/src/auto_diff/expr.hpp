#ifndef AUTO_DIFF_EXPR_HPP_
#define AUTO_DIFF_EXPR_HPP_

#include <cstddef>
#include <tuple>
#include <type_traits>

#include "auto_diff/common.hpp"
#include "auto_diff/op.hpp"
#include "auto_diff/types.hpp"

namespace halo_pm {
namespace ad {

// =============== Expressions ===============
struct ZeroExpr {};

struct OneExpr {};

template <class T>
struct ScalarExpr {
  T val_;
};

template <class... V_Expr>
struct VecExpr {
  std::tuple<V_Expr...> val_;
};

template <class Op, class V_Expr>
struct UnaryExpr {
  Op op_;
  V_Expr val_;

  UnaryExpr(V_Expr&& v) : op_{}, val_{ std::forward<V_Expr>(v) } {}
};


template <class Op, class L_Expr, class R_Expr>
struct BinaryExpr {
  Op op_;
  L_Expr l_;
  R_Expr r_;

  BinaryExpr(L_Expr&& l, R_Expr&& r) : op_{}, l_{ std::forward<L_Expr>(l) }, r_{ std::forward<R_Expr>(r) } {}

  operator Varf() { return Varf{ Evaluate(*this) }; }
};


// =============== Operator overloads ===============
template <class L_Expr, class R_Expr>
AddExpr<L_Expr, R_Expr> operator+(L_Expr&& l, R_Expr&& r) {
  return AddExpr<L_Expr, R_Expr>{ std::forward<L_Expr>(l), std::forward<R_Expr>(r) };
}

}  // namespace ad
}  // namespace halo_pm

#endif  // AUTO_DIFF_EXPR_HPP_
