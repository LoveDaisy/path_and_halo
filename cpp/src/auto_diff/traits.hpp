#ifndef AUTO_DIFF_TRAITS_HPP_
#define AUTO_DIFF_TRAITS_HPP_

#include <cstddef>
#include <type_traits>

#include "auto_diff/expr.hpp"
#include "auto_diff/types.hpp"

namespace halo_pm {
namespace ad {

namespace internal {

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

}  // namespace internal


template <class T>
constexpr bool IsArithmeticType() {
  // Normal arithmetic type
  if constexpr (std::is_arithmetic_v<T>) {
    return true;
  }

  // VecType
  else if constexpr (internal::is_vec_type_v<T> && IsArithmeticType<internal::value_type_t<T>>()) {
    return true;
  }

  // MatType
  else if constexpr (internal::is_mat_type_v<T> && IsArithmeticType<internal::value_type_t<T>>()) {
    return true;
  }

  // Other cases
  else {
    return false;
  }
}

template <class T>
constexpr bool IsVecType() {
  return internal::is_vec_type_v<T>;
}

}  // namespace ad
}  // namespace halo_pm

#endif
