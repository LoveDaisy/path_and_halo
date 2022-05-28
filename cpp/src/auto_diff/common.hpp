#ifndef AUTO_DIFF_COMMON_HPP_
#define AUTO_DIFF_COMMON_HPP_

#include <cstddef>

namespace halo_pm {
namespace ad {

// =============== Types ===============
struct NoneType;

template <class T, size_t Len>
struct Vec;
using Vec3f = Vec<float, 3>;

template <class T, size_t R, size_t C>
struct Mat;
using Mat3x3f = Mat<float, 3, 3>;


// =============== Special Op ===============
struct NoneOp;
struct IdenticalOp;

// =============== Unary Op ===============
struct NegativeOp;
struct PositiveOp;

// =============== Functional Op ===============
struct SinOp;
struct CosOp;
struct TanOp;
struct ExpOp;
struct LogOp;

// =============== Binary Op ===============
struct AddOp;
struct MinusOp;
struct TimesOp;
struct DivideOp;


// =============== Expressions ===============
struct ZeroExpr;
struct OneExpr;

template <class V>
struct VarExpr;

template <class... V>
struct VecExpr;


template <class Op, class V>
struct UnaryExpr;


template <class Op, class L, class R>
struct BinaryExpr;

template <class L, class R>
using AddExpr = BinaryExpr<AddOp, L, R>;

// =============== Evaluate ===============
constexpr float Evaluate(float);
constexpr float Evaluate(ZeroExpr);
constexpr float Evaluate(OneExpr);

template <class... T>
constexpr auto Evaluate(VecExpr<T...> e) -> Vec<float, sizeof...(T)>;

template <class T>
constexpr auto Evaluate(VarExpr<T> e);

template <class L, class R>
constexpr auto Evaluate(AddExpr<L, R> e) -> float;

}  // namespace ad

}  // namespace halo_pm

#endif  // AUTO_DIFF_COMMON_HPP_
