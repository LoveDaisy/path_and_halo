#ifndef AUTO_DIFF_COMMON_HPP_
#define AUTO_DIFF_COMMON_HPP_

#include <cstddef>

namespace halo_pm {
namespace ad {

// =============== Types ===============
struct NoneType;

template <class T, size_t Len>
struct Vec;

template <class T, size_t R, size_t C>
struct Mat;

using Vec3f = Vec<float, 3>;
using Mat3x3f = Mat<float, 3, 3>;


// =============== Operators ===============
// Special Op
struct NoneOp;
struct IdenticalOp;

// Unary Op
struct NegativeOp;
struct PositiveOp;

// Binary Op
struct AddOp;
struct MinusOp;
struct TimesOp;
struct DivideOp;

// Functional Op
struct SinOp;
struct CosOp;
struct TanOp;
struct ExpOp;
struct LogOp;


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


// =============== Expression alias ===============
template <class L, class R>
using AddExpr = BinaryExpr<AddOp, L, R>;

template <class L, class R>
using MinusExpr = BinaryExpr<MinusOp, L, R>;

template <class L, class R>
using TimesExpr = BinaryExpr<TimesOp, L, R>;

template <class L, class R>
using DivideExpr = BinaryExpr<DivideOp, L, R>;


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
