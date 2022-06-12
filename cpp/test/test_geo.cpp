#include <gtest/gtest.h>

#include <cstddef>

#include "core/geo.hpp"
#include "core/types.hpp"

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
    basis2[i] = quat * basis0[i];
  }

  for (size_t i = 0; i < kNum; i++) {
    EXPECT_NEAR((expect_basis[i] - basis1[i]).norm(), 0.0f, 1e-5);
    EXPECT_NEAR((expect_basis[i] - basis2[i]).norm(), 0.0f, 1e-5);
  }
}

}  // namespace
