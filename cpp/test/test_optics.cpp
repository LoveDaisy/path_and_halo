#include <gtest/gtest.h>

#include <vector>

#include "core/crystal.hpp"
#include "core/geo.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestOptics : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestOptics, trace_direction_0) {
  auto crystal = MakePrismCrystal(1.0f);
  Quatf q0{ 1.0f, 0.0f, 0.0f, 0.0f };
  Vec2f ray_in_ll{ -40.0f + 180.0f, 0 };
  Vec2f ray_out_ll{ 161.934798401879f, 0 };
  auto ray_out_xyz0 = Ll2Xyz(ray_out_ll);
  std::vector<int> raypath = { 3, 5 };

  auto ray_out_xyz = TraceDirection(crystal, q0, ray_in_ll, raypath);
  LOG_DEBUG("ray_out_xyz: %s", ObjLogFormatter<Vec3f>{ ray_out_xyz }.Format());
  ASSERT_NEAR((ray_out_xyz - ray_out_xyz0).norm(), 0, 1e-5);
}


// NOLINTNEXTLINE
TEST_F(TestOptics, trace_direction_1) {
  auto crystal = MakePrismCrystal(1.0f);
  Quatf q{ 0.424935669519169, -0.480586073753202, 0.669476669301674, 0.374523285804727 };
  Vec2f ray_in_ll{ 180, -15 };
  std::vector<int> raypath = { 1, 3, 2, 4, 5, 1 };

  Vec2f ray_out_ll{ 5, 25 };
  auto ray_out_xyz0 = Ll2Xyz(ray_out_ll);

  auto ray_out_xyz = TraceDirection(crystal, q, ray_in_ll, raypath);
  LOG_DEBUG("ray_out_xyz: %s", ObjLogFormatter<Vec3f>{ ray_out_xyz }.Format());
  ASSERT_NEAR((ray_out_xyz - ray_out_xyz0).norm(), 0, 1e-5);
}

}  // namespace
