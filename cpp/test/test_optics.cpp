#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <tuple>
#include <vector>

#include "auto_diff/ad.hpp"
#include "core/types.hpp"
#include "geo/geo.hpp"
#include "optics/crystal.hpp"
#include "optics/system.hpp"
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
  q.normalize();
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
    EXPECT_NEAR((j - jac2.col(i)).norm(), 0, 3e-3);
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


// NOLINTNEXTLINE
TEST_F(TestOptics, all_contours) {
  auto crystal = MakePrismCrystal(1.0f);
  std::vector<int> raypath = { 3, 5 };
  Vec2f ray_in_ll{ 0, -15 };  // sun at (180, 15)
  auto config = MakeConfigData(crystal, ray_in_ll, raypath, 3);

  auto optics_system = [&crystal, &ray_in_ll, &raypath](const Vec4f& rot) -> std::tuple<Vec4f, Mat4x4f> {
    Quatf q{ rot(0), rot(1), rot(2), rot(3) };  // w, x, y, z
    auto [xyz, j] = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);
    Vec4f out;
    Mat4x4f jac;
    out << xyz, rot.squaredNorm();
    jac << j, 2 * rot.transpose();
    return std::make_tuple(out, jac);
  };

  Vec2f target_ll{ 24.5, -15 };
  auto t0 = std::chrono::high_resolution_clock::now();
  auto [contours, status] = FindAllCrystalPoses(optics_system, target_ll, config);
  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
  LOG_INFO("time elapsed: %dus = %.3fms", dt.count(), dt.count() / 1000.0);

  LOG_INFO("contours: %zu, func_eval_cnt: %zu", contours.size(), status.func_eval_cnt_);
  for (size_t i = 0; i < contours.size(); i++) {
    const auto& c = contours[i];
    LOG_INFO("contour[%zu].size: %zu", i, c.size());
    for (const auto& q : c) {
      LOG_DEBUG("q: %s", ObjLogFormatter<Vec4f>{ q }.Format());
    }
  }
}


// NOLINTNEXTLINE
TEST_F(TestOptics, weight_component) {
  Vec4f qv{ -0.273748983855590, -0.0131274424612287, 0.000834851418692518, 0.961711219661580 };
  qv.normalize();
  auto crystal = MakePrismCrystal(1.0f);
  std::vector<int> raypath = { 3, 5 };
  Vec2f ray_in_ll{ 0, -15 };  // sun at (180, 15)
  AxisPdf a(0, 0.5);
  Func<float, 1, 3> axis_pdf = [&a](const Vec3f& llr) { return Vec<float, 1>{ a(llr) }; };

  auto data = ComputeWeightComponents(qv, crystal, axis_pdf, ray_in_ll, raypath);
  LOG_DEBUG("s: %.6f, axis_p: %.4e, jac: %.4f, geo: %.4f, transit: %.4f",  //
            data.s_, data.axis_prob_, data.jac_factor_, data.geo_factor_, data.transit_factor_);

  EXPECT_NEAR(data.axis_prob_, 3.5348, 5e-4);
  EXPECT_NEAR(data.jac_factor_, 2.3640, 0.05);
  EXPECT_NEAR(data.geo_factor_, 0.2460, 5e-4);
  EXPECT_NEAR(data.transit_factor_, 1.1558, 2e-3);
}


// NOLINTNEXTLINE
TEST_F(TestOptics, contour_weight) {
  // Data from above find_contour case
  auto crystal = MakePrismCrystal(1.0f);
  std::vector<int> raypath = { 3, 5 };
  Vec2f ray_in_ll{ 0, -15 };  // sun at (180, 15)
  AxisPdf a(0, 0.5);
  Func<float, 1, 3> axis_pdf = [&a](const Vec3f& llr) { return Vec<float, 1>{ a(llr) }; };
  Curve4f rot_contour{
    { -0.252592, 0.056621, -0.002378, 0.965912 },  //
    { -0.244650, 0.105874, -0.003684, 0.963807 },  //
    { -0.242699, 0.155373, -0.005351, 0.957563 },  //
    { -0.246603, 0.204025, -0.008339, 0.947360 },  //
    { -0.257002, 0.250470, -0.013767, 0.933287 },  //
    { -0.275563, 0.292183, -0.023002, 0.915513 },  //
    { -0.288781, 0.309657, -0.029484, 0.905454 },  //
    { -0.305153, 0.323225, -0.037266, 0.894996 },  //
    { -0.324682, 0.330994, -0.045930, 0.884824 },  //
    { -0.346360, 0.331203, -0.054373, 0.876003 },  //
    { -0.368140, 0.323677, -0.061148, 0.869463 },  //
    { -0.388190, 0.310088, -0.065401, 0.865376 },  //
    { -0.405717, 0.292523, -0.067109, 0.863318 },  //
    { -0.420651, 0.272519, -0.066645, 0.862753 },  //
    { -0.433174, 0.251019, -0.064437, 0.863248 },  //
    { -0.443501, 0.228583, -0.060856, 0.864495 },  //
    { -0.451816, 0.205547, -0.056209, 0.866286 },  //
    { -0.462896, 0.158466, -0.044657, 0.870988 },  //
    { -0.463884, 0.063405, -0.017528, 0.883450 },  //
    { -0.452911, 0.017145, -0.004493, 0.891379 },  //
    { -0.432388, -0.025908, 0.006162, 0.901294 },  //
    { -0.417762, -0.044833, 0.009916, 0.907395 },  //
    { -0.399830, -0.060488, 0.012158, 0.914510 },  //
    { -0.378753, -0.070942, 0.012595, 0.922688 },  //
    { -0.355758, -0.074059, 0.011269, 0.931570 },  //
    { -0.333149, -0.068977, 0.008786, 0.940306 },  //
    { -0.312938, -0.056934, 0.006002, 0.948046 },  //
    { -0.295880, -0.040044, 0.003486, 0.954379 },  //
    { -0.281867, -0.020058, 0.001447, 0.959242 },  //
    { -0.270526, 0.001875, -0.000112, 0.962710 },  //
    { -0.261484, 0.025039, -0.001274, 0.964882 },  //
    { -0.254436, 0.048981, -0.002145, 0.965846 },  //
    { -0.252592, 0.056621, -0.002378, 0.965912 },  //
  };

  auto [weight, data] = ComputPoseWeight(rot_contour, crystal, axis_pdf, ray_in_ll, raypath);
  LOG_INFO("contour weight: %.6f", weight);
}

}  // namespace
