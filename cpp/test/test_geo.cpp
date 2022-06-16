#include <gtest/gtest.h>

#include <cstddef>

#include "core/geo.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestGeo : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestGeo, llr2mat) {
  constexpr size_t kNum = 1;
  Vec3f llr[kNum]{ { 20, 35, -12 } };
  Mat3x3f expect_mat[kNum]{ Mat3x3f{ { -0.222484786672610f, -0.598317403666930f, 0.769751131320057f },  //
                                     { 0.959945094975663f, 0.00348527442794544f, 0.280166499593236f },  //
                                     { -0.170311286564948f, 0.801251606757469f, 0.573576436351046f } } };

  Mat3x3f mat[kNum];
  Llr2Mat(llr, mat);
  for (size_t i = 0; i < kNum; i++) {
    ASSERT_NEAR((mat[i] - expect_mat[i]).norm(), 0.0, 1e-6);  // Eigen matrix is col-major
  }
}


// NOLINTNEXTLINE
TEST_F(TestGeo, quat_rotate) {
  constexpr size_t kNum = 3;
  Quatf quat{
    0.581931465918965,   // w
    -0.223860169831750,  // xi
    -0.403854436879648,  // yj
    -0.669435573560987,  // zk
  };
  LOG_DEBUG("quaternion: [%.6f,%.6f,%.6f,%.6f]", quat.w(), quat.x(), quat.y(), quat.z());
  Vec3f basis0[kNum] = { { 1.0f, 0.0f, 0.0f },  //
                         { 0.0f, 1.0f, 0.0f },  //
                         { 0.0f, 0.0f, 1.0f } };
  Vec3f expect_basis[kNum]{ { -0.222484786672610, -0.598317403666930, 0.769751131320057 },  //
                            { 0.959945094975663, 0.00348527442794544, 0.280166499593236 },  //
                            { -0.170311286564948, 0.801251606757469, 0.573576436351046 } };
  Vec3f basis1[kNum];
  Vec3f basis2[kNum];

  RotateByQuat(quat, basis0, basis1, kNum);
  for (size_t i = 0; i < kNum; i++) {
    LOG_DEBUG("original (%zu): %s", i, ObjLogFormatter<Vec3f>{ basis0[i] }.Format());
    basis2[i] = quat * basis0[i];
    LOG_DEBUG("rotated (%zu): %s", i, ObjLogFormatter<Vec3f>{ basis2[i] }.Format());
  }

  for (size_t i = 0; i < kNum; i++) {
    EXPECT_NEAR((expect_basis[i] - basis1[i]).norm(), 0.0f, 1e-5);
    EXPECT_NEAR((expect_basis[i] - basis2[i]).norm(), 0.0f, 1e-5);
  }
}


// NOLINTNEXTLINE
TEST_F(TestGeo, normalize) {
  Vec3f v{ -0.4423, 0.1474, -0.8847 };
  v.normalize();
  LOG_DEBUG("v_norm: %.6f", v.norm());

  constexpr float kD = 1e-4;
  Vec3f vd[3]{ v + Vec3f{ kD, 0, 0 }, v + Vec3f{ 0, kD, 0 }, v + Vec3f{ 0, 0, kD } };
  for (auto& x : vd) {
    x.normalize();
  }

  Mat3x3f m = VecNormalizeDiff(v);
  for (int i = 0; i < 3; i++) {
    auto j = (vd[i] - v) / kD;
    LOG_DEBUG("j: %s", ObjLogFormatter<Vec3f>{ j }.Format());
    LOG_DEBUG("m.col(%d): %s", i, ObjLogFormatter<Vec3f>{ m.col(i) }.Format());
    EXPECT_NEAR((j - m.col(i)).norm(), 0, 5e-4);
  }
}


}  // namespace
