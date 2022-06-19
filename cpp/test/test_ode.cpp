#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "core/crystal.hpp"
#include "core/geo.hpp"
#include "core/optics.hpp"
#include "core/types.hpp"
#include "ode/solver.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestOde : public ::testing::Test {};

std::tuple<Vec2f, Mat<float, 2, 2>> test_f22(const Vec2f& x) {
  // y1 = x1^2 + x2^2
  // y2 = 2 * x1 * x2
  Vec<float, 2> y{ x.x() * x.x() + x.y() * x.y(), 2 * x.y() * x.x() };
  Mat<float, 2, 2> jac{ { 2.0f * x.x(), 2.0f * x.y() }, { 2.0f * x.y(), 2.0f * x.x() } };
  return std::make_tuple(y, jac);
};


// NOLINTNEXTLINE
TEST_F(TestOde, test_f22) {
  Vec2f x0{ 1.2, 1.1 };
  Vec2f yq{ 2.96, 2.8 };
  Vec2f expected_x{ 1.4, 1.0 };

  auto [x, status] = FindSolution<float, 2, 2>(test_f22, x0, yq);
  LOG_INFO("solved: %d, x: %s", status.solved_, ObjLogFormatter<Vec2f>{ x }.Format());
  LOG_INFO("func_eval_cnt: %zu", status.func_eval_cnt_);
  ASSERT_TRUE(status.solved_);
  ASSERT_NEAR((x - expected_x).norm(), 0, 1e-5);
}


// NOLINTNEXTLINE
TEST_F(TestOde, trace_crystal) {
  auto crystal = MakePrismCrystal(1.0f);
  std::vector<int> raypath{ 1, 3, 2, 4, 5, 1 };

  Vec2f ray_in_ll{ 180, -15 };
  Vec2f ray_out_ll{ 5, 25 };
  auto out_xyz = Ll2Xyz(ray_out_ll);

  Vec4f yq;
  yq << out_xyz, 1.0f;

  Vec4f vq0{ 0.4, -0.5, 0.7, 0.4 };  // [wxyz]
  vq0.normalize();
  Vec4f expected_q{ 0.424935669519169, -0.480586073753202, 0.669476669301674, 0.374523285804727 };
  expected_q.normalize();

  auto fun = [&crystal, &ray_in_ll, &raypath](const Vec4f& q) {
    Quatf rot{ q(0), q(1), q(2), q(3) };
    auto [ray_out_xyz, j] = TraceDirDiffQuat(crystal, rot, ray_in_ll, raypath);

    return std::make_tuple(ray_out_xyz, j);
  };
  auto [x, status] = FindSolution<float, 3, 4>(fun, vq0, out_xyz);
  auto [y, j] = fun(x);

  LOG_DEBUG("solved: %d, fun_val_cnt: %zu", status.solved_, status.func_eval_cnt_);
  LOG_DEBUG("x: %s", ObjLogFormatter<Vec4f>{ x }.Format());

  ASSERT_TRUE(status.solved_);
  ASSERT_NEAR((out_xyz - y).norm(), 0, 1e-5);
}

}  // namespace
