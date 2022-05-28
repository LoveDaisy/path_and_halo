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

template <class V>
struct VarExpr {
  V val_;

  VarExpr(V&& v) : val_(std::forward<V>(v)) {}
};

template <class... V>
struct VecExpr {
  std::tuple<V...> val_;
};

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


}  // namespace ad
}  // namespace halo_pm

#endif  // AUTO_DIFF_EXPR_HPP_
