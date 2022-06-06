#include <gtest/gtest.h>

#include <tuple>

#include "auto_diff/ad.hpp"
#include "auto_diff/common.hpp"
#include "util/log.hpp"

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

    LOG_DEBUG("x1_val: %.4f", x1_val);

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

  LOG_DEBUG("simple_lazy done!!");
}


// NOLINTNEXTLINE
TEST_F(TestAD, simple_diff) {
  LOG_DEBUG("test!!");

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

  {
    ad::VarExpr u = a * b / (a + b);
    ad::VarExpr ua = ad::Differentiate(u, ad::wrt(a));
    float ua_val = ad::Evaluate(ua);
    ASSERT_NEAR(ua_val, b.val_ * b.val_ / (a.val_ + b.val_) / (a.val_ + b.val_), 1e-5);
  }
}


std::tuple<float, float, float> f1(float a, float b) {
  ad::VarExpr var_a = a;
  ad::VarExpr var_b = b;
  auto u = var_a * var_b / (var_a + var_b);
  auto ua = ad::Differentiate(u, ad::wrt(var_a));
  auto ub = ad::Differentiate(u, ad::wrt(var_b));
  return std::make_tuple(ad::Evaluate(u), ad::Evaluate(ua), ad::Evaluate(ub));
}


// NOLINTNEXTLINE
TEST_F(TestAD, simple_fun) {
  float a = 2.3f;
  float b = 3.4f;
  auto [u, ua, ub] = f1(a, b);
  ASSERT_NEAR(u, a * b / (a + b), 1e-5);
  ASSERT_NEAR(ua, b * b / (a + b) / (a + b), 1e-5);
  ASSERT_NEAR(ub, a * a / (a + b) / (a + b), 1e-5);
}

}  // namespace
