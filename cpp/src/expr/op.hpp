#ifndef AUTO_DIFF_OP_HPP_
#define AUTO_DIFF_OP_HPP_

namespace halo_pm {
namespace internal {

// =============== Special Op ===============
struct NoneOp {};
struct IdenticalOp {};

// =============== Unary Op ===============
struct NegativeOp {};
struct PositiveOp {};

// =============== Functional Op ===============
struct SinOp {};
struct CosOp {};
struct TanOp {};
struct ExpOp {};
struct LogOp {};

// =============== Binary Op ===============
struct AddOp {};
struct MinusOp {};
struct TimesOp {};
struct DivideOp {};

}  // namespace internal
}  // namespace halo_pm

#endif  // AUTO_DIFF_OP_HPP_
