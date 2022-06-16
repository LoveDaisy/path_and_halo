#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>

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
TEST_F(TestOde, test_find_sol_f22) {
  Vec2f x0{ 1.2, 1.1 };
  Vec2f yq{ 2.96, 2.8 };
  Vec2f expected_x{ 1.4, 1.0 };

  auto [x, status] = FindSolution<float, 2, 2>(test_f22, x0, yq);
  LOG_INFO("solved: %d, x: %s", status.solved_, ObjLogFormatter<Vec2f>{ x }.Format());
  LOG_INFO("func_eval_cnt: %zu", status.func_eval_cnt_);
  ASSERT_TRUE(status.solved_);
  ASSERT_NEAR((x - expected_x).norm(), 0, 1e-5);
}

}  // namespace
