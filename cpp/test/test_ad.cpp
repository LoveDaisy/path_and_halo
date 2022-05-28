#include <gtest/gtest.h>

#include "auto_diff/ad.hpp"
#include "auto_diff/common.hpp"
#include "auto_diff/expr.hpp"

namespace {
// NOLINTNEXTLINE
using namespace halo_pm;

class TestAD : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestAD, add) {
  Varf a{ 1.0f };
  Varf b{ 2.0f };
  Varf c{ -3.5f };
  auto x = a + b + c;  // Not evaluate here
  auto x_val = ad::Evaluate(x);

  Varf y = a + b + c;  // Evaluate here

  ASSERT_NEAR(x_val, -0.5f, 1e-6);
  ASSERT_NEAR(y.val_, -0.5f, 1e-6);
}

}  // namespace
