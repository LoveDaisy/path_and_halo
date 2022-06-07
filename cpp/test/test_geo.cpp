#include <gtest/gtest.h>

#include <cstddef>

#include "core/geo.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestGeo : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestGeo, llr2mat) {
  float llr[3]{ 20, 35, -12 };
  float expect_mat[9]{ -0.222484786672610, -0.598317403666930,  0.769751131320057,  //
                       0.959945094975663,  0.00348527442794544, 0.280166499593236,  //
                       -0.170311286564948, 0.801251606757469,   0.573576436351046 };

  float mat[9];
  Llr2Mat(llr, mat);
  for (size_t i = 0; i < 9; i++) {
    ASSERT_NEAR(mat[i], expect_mat[i], 1e-6);
  }
}


// NOLINTNEXTLINE
TEST_F(TestGeo, quat_rotate) {
  float quat[4]{ 0.581931465918965, -0.223860169831750, -0.403854436879648, -0.669435573560987 };
  float basis0[9] = { 1.0f, 0.0f, 0.0f,  //
                      0.0f, 1.0f, 0.0f,  //
                      0.0f, 0.0f, 1.0f };
  float expect_mat[9]{ -0.222484786672610, -0.598317403666930,  0.769751131320057,  //
                       0.959945094975663,  0.00348527442794544, 0.280166499593236,  //
                       -0.170311286564948, 0.801251606757469,   0.573576436351046 };

  float basis1[9];
  RotateByQuat(quat, basis0, basis1, 3);
  for (size_t i = 0; i < 9; i++) {
    EXPECT_NEAR(basis1[i], expect_mat[i], 1e-6);
  }
}

}  // namespace
