#include <gtest/gtest.h>

#include "auto_diff/ad.hpp"
#include "auto_diff/common.hpp"

namespace {
// NOLINTNEXTLINE
using namespace halo_pm;

class TestAD : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestAD, simple_lazy) {
  ad::VarExpr a1 = 1.0f;
  ad::VarExpr b1 = 2.4f;
  ad::VarExpr c1 = -0.8f;

  {
    ad::VarExpr x1 = a1 + b1 + c1;  // Explicitly assign result as VarExpr
    auto x2 = a1 + b1 + c1;         // or let it be auto (thus **NOT** a VarExpr type)
    float x1_val = ad::Evaluate(x1);
    float x2_val = ad::Evaluate(x2);

    ASSERT_NEAR(x1_val, 2.6f, 1e-6);
    ASSERT_NEAR(x2_val, 2.6f, 1e-6);
  }

  {
    ad::VarExpr x1 = a1 - b1;
    float x1_val = ad::Evaluate(x1);

    ASSERT_NEAR(x1_val, -1.4f, 1e-6);
  }

  {
    ad::VarExpr x1 = a1 * b1 * c1;
    float x1_val = ad::Evaluate(x1);

    ASSERT_NEAR(x1_val, -1.92f, 1e-6);
  }

  {
    ad::VarExpr x1 = a1 - b1 * c1;
    float x1_val = ad::Evaluate(x1);

    ASSERT_NEAR(x1_val, 2.92f, 1e-6);
  }
}

// NOLINTNEXTLINE
TEST_F(TestAD, simple_diff) {
  ad::VarExpr a = 1.0f;
  ad::VarExpr b = 2.3f;
  ad::VarExpr c = -.3f;

  {
    ad::VarExpr u = a + b + c;
    float u_jac = ad::Differentiate(u, ad::wrt(a));
    ASSERT_NEAR(u_jac, 1.0f, 1e-5);
  }

  {
    ad::VarExpr u = a / b;
    ad::VarExpr ua = ad::Differentiate(u, ad::wrt(a));
    float ua_val = ad::Evaluate(ua);
    ASSERT_NEAR(ua_val, 1.0f / 2.3f, 1e-5);
  }

  {
    ad::VarExpr u = a / b * c;
    ad::VarExpr ua = ad::Differentiate(u, ad::wrt(a));
    float ua_val = ad::Evaluate(ua);
    ASSERT_NEAR(ua_val, -0.3f / 2.3f, 1e-5);

    ad::VarExpr ub = ad::Differentiate(u, ad::wrt(b));
    float ub_val = ad::Evaluate(ub);
    ASSERT_NEAR(ub_val, 0.3f / 2.3f / 2.3f, 1e-5);
  }
}

}  // namespace
