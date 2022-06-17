#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "core/types.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

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

  Eigen::JacobiSVD svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);

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


// NOLINTNEXTLINE
TEST_F(TestEigen, test_quaternion) {
  // NOTE: quaternion in Eigen (and wiki) seems different to that in MATLAB. The vector parts are opposite.
  Quatf q{ 0.424935669519169, -0.480586073753202, 0.669476669301674, 0.374523285804727 };
  Vec3f rotated_v[3]{
    { -0.176933404857992, -0.325185721574301, -0.928950384428339 },  // [1,0,0]
    { -0.961778934476391, 0.257538668108388, 0.0930328739017340 },   // [0,1,0]
    { 0.208987682514584, 0.909905534060080, -0.358323970233677 },    // [0,0,1]
  };

  Vec3f v[3]{ q * Vec3f{ 1, 0, 0 }, q * Vec3f{ 0, 1, 0 }, q * Vec3f{ 0, 0, 1 } };
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR((rotated_v[i] - v[i]).norm(), 0, 1e-5);
  }
}

}  // namespace
