#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "auto_diff/ad.hpp"
#include "core/geo.hpp"
#include "core/math.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestGeo : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestGeo, llr_rot) {
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
TEST_F(TestGeo, llr_quat) {
  Vec3f llr{ 20, 35, -12 };
  Quatf q = Llr2Quat(llr);
  LOG_DEBUG("q: %s", ObjLogFormatter<Quatf>{ q }.Format());

  Vec3f basis[3]{ { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f } };
  Vec3f expect_basis[3]{ { 0.598317, -0.00348527, -0.801252 },  //
                         { -0.222485, 0.959945, -0.170311 },    //
                         { 0.769751, 0.280166, 0.573576 } };
  for (int i = 0; i < 3; i++) {
    auto b = q * basis[i];
    LOG_DEBUG("expect_basis(%d): %s", i, ObjLogFormatter<Vec3f>{ expect_basis[i] }.Format());
    LOG_DEBUG("b: %s", ObjLogFormatter<Vec3f>{ b }.Format());
    EXPECT_NEAR((b - expect_basis[i]).norm(), 0, 1e-5);
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
TEST_F(TestGeo, normalize_Vec3f) {
  Vec3f v{ -0.4423, 0.1474, -0.8847 };
  v.normalize();
  LOG_DEBUG("v_norm: %.6f", v.norm());

  constexpr float kD = 1e-4;
  Vec3f vd[3]{ v + Vec3f{ kD, 0, 0 }, v + Vec3f{ 0, kD, 0 }, v + Vec3f{ 0, 0, kD } };
  for (auto& x : vd) {
    x.normalize();
  }

  ad::Var vx = v.x();
  ad::Var vy = v.y();
  ad::Var vz = v.z();

  auto [ox, oy, oz] = NormalizeExpr(vx, vy, vz);
  Mat3x3f m{ { ad::Diff(ox, ad::wrt(vx)), ad::Diff(ox, ad::wrt(vy)), ad::Diff(ox, ad::wrt(vz)) },
             { ad::Diff(oy, ad::wrt(vx)), ad::Diff(oy, ad::wrt(vy)), ad::Diff(oy, ad::wrt(vz)) },
             { ad::Diff(oz, ad::wrt(vx)), ad::Diff(oz, ad::wrt(vy)), ad::Diff(oz, ad::wrt(vz)) } };

  for (int i = 0; i < 3; i++) {
    auto j = (vd[i] - v) / kD;
    LOG_DEBUG("j: %s", ObjLogFormatter<Vec3f>{ j }.Format());
    LOG_DEBUG("m.col(%d): %s", i, ObjLogFormatter<Vec3f>{ m.col(i) }.Format());
    EXPECT_NEAR((j - m.col(i)).norm(), 0, 6e-4);
  }
}


// NOLINTNEXTLINE
TEST_F(TestGeo, normalize_Quatf) {
  Quatf q{ -0.4423, 0.1474, -0.8847, 0.2345 };
  q.normalize();

  using ad::Diff;
  using ad::Var;
  using ad::wrt;

  Var qw = q.w();
  Var qx = q.x();
  Var qy = q.y();
  Var qz = q.z();

  auto [ow, ox, oy, oz] = NormalizeExpr(qw, qx, qy, qz);
  LOG_DEBUG("d(qx)/d_qx: %.6f", Diff(ox, wrt(qx)));

  Mat4x4f m{ { Diff(ow, wrt(qw)), Diff(ow, wrt(qx)), Diff(ow, wrt(qy)), Diff(ow, wrt(qz)) },
             { Diff(ox, wrt(qw)), Diff(ox, wrt(qx)), Diff(ox, wrt(qy)), Diff(ox, wrt(qz)) },
             { Diff(oy, wrt(qw)), Diff(oy, wrt(qx)), Diff(oy, wrt(qy)), Diff(oy, wrt(qz)) },
             { Diff(oz, wrt(qw)), Diff(oz, wrt(qx)), Diff(oz, wrt(qy)), Diff(oz, wrt(qz)) } };

  LOG_DEBUG("jac:\n%s", ObjLogFormatter<Mat4x4f>{ m }.Format());
}


// NOLINTNEXTLINE
TEST_F(TestGeo, quat_rot_expr) {
  Vec3f v{ 1.2f, 3.2f, -0.3f };
  Quatf q{ 1.0f, 2.0f, -2.0f, 0.2f };
  q.normalize();

  LOG_DEBUG("v: %s", ObjLogFormatter<Vec3f>{ v }.Format());
  LOG_DEBUG("q: %s", ObjLogFormatter<Quatf>{ q }.Format());

  using ad::Diff;
  using ad::Eval;
  using ad::Var;
  using ad::wrt;

  Var vx = v.x();
  Var vy = v.y();
  Var vz = v.z();
  Var qw = q.w();
  Var qx = q.x();
  Var qy = q.y();
  Var qz = q.z();

  auto [ox, oy, oz] = QuatRotExpr(qw, qx, qy, qz, vx, vy, vz);
  float ox_val = Eval(ox);
  float oy_val = Eval(oy);
  float oz_val = Eval(oz);
  EXPECT_NEAR(ox_val, -2.739823, 1e-6);
  EXPECT_NEAR(oy_val, -0.509735, 1e-6);
  EXPECT_NEAR(oz_val, 2.000885, 1e-6);
  LOG_DEBUG("v_out: [%.6f,%.6f,%.6f]", Eval(ox), Eval(oy), Eval(oz));
  LOG_DEBUG("d/d_vx: [%.6f,%.6f,%.6f]", Diff(ox, wrt(vx)), Diff(oy, wrt(vx)), Diff(oz, wrt(vx)));
  LOG_DEBUG("d/d_qw: [%.6f,%.6f,%.6f]", Diff(ox, wrt(qw)), Diff(oy, wrt(qw)), Diff(oz, wrt(qw)));
}


// NOLINTNEXTLINE
TEST_F(TestGeo, interp_spline) {
  constexpr size_t kXNum = 8;
  constexpr size_t kQNum = 64;

  std::vector<float> x(kXNum);
  std::vector<Vec2f> y(kXNum);
  for (size_t i = 0; i < kXNum; i++) {
    x[i] = 2 * 3.1415926 / (kXNum - 1) * i;
    y[i] = Vec2f{ std::sin(x[i]), std::cos(x[i]) };
  }

  std::vector<float> xq(kQNum);
  for (size_t i = 0; i < kQNum; i++) {
    xq[i] = 2 * 3.1415926 / (kQNum - 1) * i;
  }

  auto yq = InterpSpline(x, y, xq);

  ASSERT_EQ(yq.size(), xq.size());
  for (size_t i = 0; i < xq.size(); i++) {
    LOG_DEBUG("xq: %.6f, yq: %s", xq[i], ObjLogFormatter<Vec2f>{ yq[i] }.Format());
    EXPECT_NEAR(yq[i].x(), std::sin(xq[i]), 0.05);
    EXPECT_NEAR(yq[i].y(), std::cos(xq[i]), 0.05);
  }


  auto y2q = InterpCurve(y, 0.05);
  for (const auto& pt : y2q) {
    LOG_DEBUG("pt: %s", ObjLogFormatter<Vec2f>{ pt }.Format());
    EXPECT_NEAR(pt.norm(), 1.0, 0.06);
  }
}

}  // namespace
