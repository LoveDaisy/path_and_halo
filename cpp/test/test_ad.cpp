#include <gtest/gtest.h>

#include <tuple>

#include "auto_diff/ad.hpp"
#include "util/log.hpp"

namespace {
// NOLINTNEXTLINE
using namespace halo_pm;

class TestAD : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestAD, simple_lazy) {
  ad::Var a1 = 1.0f;
  ad::Var b1 = 2.4f;
  ad::Var c1 = -0.8f;

  {
    ad::Var x1 = a1 + b1 + c1;  // Explicitly assign result as VarExpr
    auto x2 = a1 + b1 + c1;     // or let it be auto (thus **NOT** a VarExpr type)
    float x1_val = ad::Eval(x1);
    float x2_val = ad::Eval(x2);

    LOG_DEBUG("x1_val: %.4f", x1_val);

    ASSERT_NEAR(x1_val, 2.6f, 1e-6);
    ASSERT_NEAR(x2_val, 2.6f, 1e-6);
  }

  {
    ad::Var x1 = a1 - b1;
    float x1_val = ad::Eval(x1);

    ASSERT_NEAR(x1_val, -1.4f, 1e-6);
  }

  {
    ad::Var x1 = a1 * b1 * c1;
    float x1_val = ad::Eval(x1);

    ASSERT_NEAR(x1_val, -1.92f, 1e-6);
  }

  {
    ad::Var x1 = a1 - b1 * c1;
    float x1_val = ad::Eval(x1);

    ASSERT_NEAR(x1_val, 2.92f, 1e-6);
  }

  LOG_DEBUG("simple_lazy done!!");
}


// NOLINTNEXTLINE
TEST_F(TestAD, simple_diff) {
  ad::Var a = 1.0f;
  ad::Var b = 2.3f;
  ad::Var c = -.3f;

  {
    ad::Var u = a + b + c;
    float u_jac = ad::Diff(u, ad::wrt(a));
    ASSERT_NEAR(u_jac, 1.0f, 1e-5);
  }

  {
    ad::Var u = a / b;
    ad::Var ua = ad::Diff(u, ad::wrt(a));
    float ua_val = ad::Eval(ua);
    ASSERT_NEAR(ua_val, 1.0f / 2.3f, 1e-5);
  }

  {
    ad::Var u = a / b * c;
    ad::Var ua = ad::Diff(u, ad::wrt(a));
    float ua_val = ad::Eval(ua);
    ASSERT_NEAR(ua_val, -0.3f / 2.3f, 1e-5);

    ad::Var ub = ad::Diff(u, ad::wrt(b));
    float ub_val = ad::Eval(ub);
    ASSERT_NEAR(ub_val, 0.3f / 2.3f / 2.3f, 1e-5);
  }

  {
    ad::Var u = a * b / (a + b);
    ad::Var ua = ad::Diff(u, ad::wrt(a));
    float ua_val = ad::Eval(ua);
    ASSERT_NEAR(ua_val, b.val_ * b.val_ / (a.val_ + b.val_) / (a.val_ + b.val_), 1e-5);
  }

  {
    auto u = a * 3.0f;
    float u_val = ad::Eval(u);
    ASSERT_NEAR(u_val, 3.0f, 1e-5);
  }
}


// NOLINTNEXTLINE
TEST_F(TestAD, op_fun_diff) {
  ad::Var vx{ 2.0f };
  ad::Var vy{ -1.0f };
  ad::Var vz{ .3f };

  {
    auto u = sqrt(vx);
    float u_val = ad::Eval(u);
    ASSERT_NEAR(u_val, 1.4142135f, 1e-5);

    auto ua = ad::Diff(u, ad::wrt(vx));
    float ua_val = ad::Eval(ua);
    ASSERT_NEAR(ua_val, 0.353553f, 1e-5);
  }

  {
    auto u = sqrt(vx + 1.0f);
    float u_val = ad::Eval(u);
    ASSERT_NEAR(u_val, 1.7320508f, 1e-5);
  }

  {
    auto c = vx * 3.2f + vy * 0.4f - vz * 0.8f;
    auto delta = 2.0f - 0.3f / (c * c);
    auto a = sqrt(delta);
    auto o = 1.3f * vx + a * 0.2;

    auto ax_val = ad::Eval(ad::Diff(a, ad::wrt(vx)));
    auto ox_val = ad::Eval(ad::Diff(o, ad::wrt(vx)));
    auto ay_val = ad::Eval(ad::Diff(a, ad::wrt(vy)));
    auto oy_val = ad::Eval(ad::Diff(o, ad::wrt(vy)));
    LOG_DEBUG("ax_val: %.6f", ax_val);
    LOG_DEBUG("ox_val: %.6f", ox_val);
    LOG_DEBUG("ay_val: %.6f", ay_val);
    LOG_DEBUG("oy_val: %.6f", oy_val);
    EXPECT_NEAR(ax_val, 0.00356019, 1e-5);
    EXPECT_NEAR(ox_val, 1.30071, 1e-5);
    EXPECT_NEAR(ay_val, 0.000445023, 1e-7);
    EXPECT_NEAR(oy_val, 0.0000890047, 1e-7);
  }
}


std::tuple<float, float, float> f1(float a, float b) {
  ad::Var var_a = a;
  ad::Var var_b = b;
  auto u = var_a * var_b / (var_a + var_b);
  auto ua = ad::Diff(u, ad::wrt(var_a));
  auto ub = ad::Diff(u, ad::wrt(var_b));
  return std::make_tuple(ad::Eval(u), ad::Eval(ua), ad::Eval(ub));
}


template <class A, class B>
auto f2(const ad::Var<A>& a, const ad::Var<B>& b) {
  auto u = a + b;
  auto v = 2 * a - b;
  auto uv = u / v;
  return std::make_tuple(u, v, uv);
}


// NOLINTNEXTLINE
TEST_F(TestAD, simple_fun) {
  {
    float a = 2.3f;
    float b = 3.4f;
    auto [u, ua, ub] = f1(a, b);
    ASSERT_NEAR(u, a * b / (a + b), 1e-5);
    ASSERT_NEAR(ua, b * b / (a + b) / (a + b), 1e-5);
    ASSERT_NEAR(ub, a * a / (a + b) / (a + b), 1e-5);
  }

  {
    ad::Var a = 2.3f;
    ad::Var b = 3.4f;
    auto [u, v, uv] = f2(a, b);
    ASSERT_NEAR(ad::Eval(u), 5.7f, 1e-6);
    ASSERT_NEAR(ad::Eval(v), 1.2f, 1e-6);
    ASSERT_NEAR(ad::Eval(uv), 4.75f, 1e-6);

    float ua = ad::Diff(u, ad::wrt(a));
    float ub = ad::Diff(u, ad::wrt(b));
    ASSERT_NEAR(ua, 1.0f, 1e-5);
    ASSERT_NEAR(ub, 1.0f, 1e-5);

    float va = ad::Diff(v, ad::wrt(a));
    float vb = ad::Diff(v, ad::wrt(b));
    ASSERT_NEAR(va, 2.0f, 1e-5);
    ASSERT_NEAR(vb, -1.0f, 1e-5);

    float uva = ad::Diff(uv, ad::wrt(a));
    float uvb = ad::Diff(uv, ad::wrt(b));
    ASSERT_NEAR(uva, -7.08333f, 1e-5);
    ASSERT_NEAR(uvb, 4.791667f, 1e-5);
  }
}

}  // namespace
