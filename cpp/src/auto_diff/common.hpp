#ifndef AUTO_DIFF_COMMON_HPP_
#define AUTO_DIFF_COMMON_HPP_

#include <cstddef>

namespace halo_pm {
namespace ad {

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
template <class V>
struct VarExpr;

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


// =============== Differentiate ===============
template <class V>
struct wrt;

}  // namespace ad

}  // namespace halo_pm

#endif  // AUTO_DIFF_COMMON_HPP_
