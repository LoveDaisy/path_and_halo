#ifndef AUTO_DIFF_TRAITS_HPP_
#define AUTO_DIFF_TRAITS_HPP_

#include <cstddef>
#include <type_traits>

#include "auto_diff/common.hpp"
#include "auto_diff/expr.hpp"
#include "auto_diff/types.hpp"

namespace halo_pm {
namespace ad {

namespace traits {

// =============== Traits ===============
// MatType
template <class T>
struct is_mat_type : public std::false_type {};

template <class T, size_t R, size_t C>
struct is_mat_type<Mat<T, R, C>> : public std::true_type {};

template <class T>
inline constexpr bool is_mat_type_v = is_mat_type<T>::value;


// VecType
template <class T>
struct is_vec_type : public std::false_type {};

template <class T, size_t Len>
struct is_vec_type<Vec<T, Len>> : public std::true_type {};

template <class... T>
struct is_vec_type<VecExpr<T...>> : public std::true_type {};

template <class T>
inline constexpr bool is_vec_type_v = is_vec_type<T>::value;


// Expr
template <class T>
struct is_expr : public std::false_type {};

template <class T>
struct is_expr<VarExpr<T>> : public std::true_type {};

template <class... T>
struct is_expr<VecExpr<T...>> : public std::true_type {};

template <class Op, class T>
struct is_expr<UnaryExpr<Op, T>> : public std::true_type {};

template <class Op, class L, class R>
struct is_expr<BinaryExpr<Op, L, R>> : public std::true_type {};

template <class T>
inline constexpr bool is_expr_v = is_expr<T>::value;


// Get underlying type of value
template <class T, class Enable = void>
struct value_type {
  using type = NoneType;
};

template <class T>
struct value_type<T, typename std::enable_if_t<std::is_arithmetic_v<T>>> {
  // Normal arithmetic types
  using type = T;
};

template <class T, size_t Len>
struct value_type<Vec<T, Len>> {
  // VecType
  using type = T;
};

template <class T, size_t R, size_t C>
struct value_type<Mat<T, R, C>> {
  // MatType
  using type = T;
};

template <class Op, class T>
struct value_type<UnaryExpr<Op, T>> {
  // UnaryExpr
  using type = typename value_type<T>::type;
};

template <class T>
struct value_type<VarExpr<T>> {
  // Termination of UnaryExpr recursion
  using type = T;
};

template <class T>
using value_type_t = typename value_type<T>::type;


// Variable dimension
template <class T, class Enable = void>
struct var_dim {
  static constexpr size_t dim = 0;
};

template <class T>
struct var_dim<VarExpr<T>, std::enable_if_t<std::is_arithmetic_v<T>>> {
  static constexpr size_t dim = 1;
};

template <class T>
struct var_dim<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
  static constexpr size_t dim = 1;
};

template <class T>
struct var_dim<wrt<T>> {
  static constexpr size_t dim = var_dim<T>::dim;
};

}  // namespace traits

}  // namespace ad
}  // namespace halo_pm

#endif
