#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "util/log.hpp"

namespace {

class TestEigen : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestEigen, test_size) {
  constexpr int kNum = 4;
  Eigen::Vector3f vector[kNum];
  LOG_INFO("sizeof Vector3f[%d]: %zu", kNum, sizeof(vector));
  ASSERT_EQ(sizeof(vector), kNum * 3 * sizeof(float));
}


// NOLINTNEXTLINE
TEST_F(TestEigen, test_svd) {
  Eigen::Matrix3f m;
  m << 0.0333, 0.1983, 0.1908,  //
      0.5724, 0.0207, 0.9521,   //
      0.0093, 0.7923, 0.7702;

  Eigen::Vector3f expected_s;
  expected_s << 1.43129, 0.695169, 0.0238823;
  Eigen::Matrix3f expected_u;
  expected_u << 0.179768, 0.144515, 0.973036,  //
      0.69289, -0.720738, -0.0209677,          //
      0.698274, 0.677976, -0.229698;
  Eigen::Matrix3f expected_v;
  expected_v << 0.285819, -0.577461, 0.764753,  //
      0.42146, 0.792467, 0.440872,              //
      0.860627, -0.196303, -0.469878;

  Eigen::JacobiSVD svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

  const auto& s = svd.singularValues();
  for (int i = 0; i < s.rows(); i++) {
    ASSERT_NEAR(s(i), expected_s(i), 1e-5);
  }

  const auto& u = svd.matrixU();
  const auto& v = svd.matrixV();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_NEAR(u(i, j), expected_u(i, j), 1e-5);
      ASSERT_NEAR(v(i, j), expected_v(i, j), 1e-5);
    }
  }
}

}  // namespace
