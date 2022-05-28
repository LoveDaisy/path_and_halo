#include <gtest/gtest.h>

#include "auto_diff/ad.hpp"

namespace {
// NOLINTNEXTLINE
using namespace halo_pm;

class TestAD : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestAD, add) {
  ad::VarExpr a1 = 1.0f;
  ad::VarExpr b1 = 2.4f;
  ad::VarExpr c1 = -0.8f;
  ad::VarExpr x1 = a1 + b1 + c1;  // Explicitly assign result as VarExpr
  auto x2 = a1 + b1 + c1;         // or let it be auto (thus **NOT** a VarExpr type)
  float x1_val = ad::Evaluate(x1);
  float x2_val = ad::Evaluate(x2);

  ASSERT_NEAR(x1_val, 2.6f, 1e-6);
  ASSERT_NEAR(x2_val, 2.6f, 1e-6);
}

}  // namespace
