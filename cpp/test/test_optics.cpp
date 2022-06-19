#include <gtest/gtest.h>

#include <vector>

#include "auto_diff/ad.hpp"
#include "core/crystal.hpp"
#include "core/geo.hpp"
#include "core/optics.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestOptics : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestOptics, trace_direction_1) {
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
TEST_F(TestOptics, trace_direction_2) {
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


// NOLINTNEXTLINE
TEST_F(TestOptics, trace_direction_diff) {
  auto crystal = MakePrismCrystal(1.0f);
  Quatf q{ 0.424935669519169, -0.480586073753202, 0.669476669301674, 0.374523285804727 };
  Vec2f ray_in_ll{ 180, -15 };
  std::vector<int> raypath = { 1, 3, 2, 4, 5, 1 };

  Vec2f ray_out_ll{ 5, 25 };
  auto ray_out_xyz0 = Ll2Xyz(ray_out_ll);
  LOG_DEBUG("expected: %s", ObjLogFormatter<Vec3f>{ ray_out_xyz0 }.Format());

  auto ray_out_xyz1 = TraceDirection(crystal, q, ray_in_ll, raypath);
  LOG_DEBUG("normal trace: %s", ObjLogFormatter<Vec3f>{ ray_out_xyz1 }.Format());

  auto [ray_out_xyz2, jac2] = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);
  LOG_DEBUG("ad trace: %s", ObjLogFormatter<Vec3f>{ ray_out_xyz2 }.Format());

  EXPECT_NEAR((ray_out_xyz0 - ray_out_xyz1).norm(), 0, 1e-5);
  EXPECT_NEAR((ray_out_xyz2 - ray_out_xyz1).norm(), 0, 1e-5);

  constexpr float kD = 3e-4;
  Quatf q_d[4]{ q, q, q, q };
  q_d[0].w() += kD;
  q_d[1].x() += kD;
  q_d[2].y() += kD;
  q_d[3].z() += kD;
  for (auto& qq : q_d) {
    qq.normalize();
  }

  for (int i = 0; i < 4; i++) {
    auto tmp_xyz = TraceDirection(crystal, q_d[i], ray_in_ll, raypath);
    auto j = (tmp_xyz - ray_out_xyz0) / kD;
    LOG_DEBUG("j: %s", ObjLogFormatter<Vec3f>{ j }.Format());
    LOG_DEBUG("jac.col(%d): %s", i, ObjLogFormatter<Vec3f>{ jac2.col(i) }.Format());
    EXPECT_NEAR((j - jac2.col(i)).norm(), 0, 1e-2);
  }
}


// NOLINTNEXTLINE
TEST_F(TestOptics, refract_diff) {
  float n = 1.31f;
  Vec3f ray_in{ -0.4423, 0.1474, -0.8847 };
  ray_in.normalize();
  Vec3f norm{ 0, 0, 1 };

  auto ray_out_1 = Refract(ray_in, norm, 1.0f, n);

  ad::Var vx = ray_in.x();
  ad::Var vy = ray_in.y();
  ad::Var vz = ray_in.z();

  auto [delta, rx, ry, rz, lx, ly, lz] = RefractExpr(vx, vy, vz, norm, 1.0f, n);
  ASSERT_GT(ad::Eval(delta), 0.0f);

  Vec3f ray_out_2{ ad::Eval(rx), ad::Eval(ry), ad::Eval(rz) };
  Mat3x3f jac_2{ { ad::Diff(rx, ad::wrt(vx)), ad::Diff(rx, ad::wrt(vy)), ad::Diff(rx, ad::wrt(vz)) },
                 { ad::Diff(ry, ad::wrt(vx)), ad::Diff(ry, ad::wrt(vy)), ad::Diff(ry, ad::wrt(vz)) },
                 { ad::Diff(rz, ad::wrt(vx)), ad::Diff(rz, ad::wrt(vy)), ad::Diff(rz, ad::wrt(vz)) } };

  LOG_DEBUG("r1: %s", ObjLogFormatter<Vec3f>{ ray_out_1 }.Format());
  LOG_DEBUG("r2: %s", ObjLogFormatter<Vec3f>{ ray_out_2 }.Format());
  ASSERT_NEAR((ray_out_1 - ray_out_2).norm(), 0, 1e-5);

  constexpr float kD = 1e-4;
  Vec3f ray_in_diff[3]{ ray_in + Vec3f{ kD, 0, 0 }, ray_in + Vec3f{ 0, kD, 0 }, ray_in + Vec3f{ 0, 0, kD } };
  for (auto& v : ray_in_diff) {
    v.normalize();
  }

  Vec3f ray_out_d[3];
  for (int i = 0; i < 3; i++) {
    ray_out_d[i] = Refract(ray_in_diff[i], norm, 1.0f, n);
    Vec3f j = (ray_out_d[i] - ray_out_1) / kD;
    LOG_DEBUG("j: %s", ObjLogFormatter<Vec3f>{ j }.Format());
    LOG_DEBUG("jac.col(%d): %s", i, ObjLogFormatter<Vec3f>{ jac_2.col(i) }.Format());
    EXPECT_NEAR((j - jac_2.col(i)).norm(), 0.0f, 1e-3);
  }
}

}  // namespace
