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
TEST_F(TestOptics, llr_output) {
  auto crystal = MakePrismCrystal(1.0f);
  Vec2f ray_in_ll{ 0, -15 };
  std::vector<int> raypath = { 3, 5 };

  {
    LOG_INFO("case 1...");
    Vec3f llr{ 12.2499996679285, 51.9668075825233, 125.527501860455 };
    Quatf q = Llr2Quat(llr);

    auto [ray_out, jac] = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);
    Vec2f out_ll = Xyz2Ll(ray_out);
    LOG_DEBUG("llr: %s, quat: %s", ObjLogFormatter<Vec3f>{ llr }.Format(), ObjLogFormatter<Quatf>{ q }.Format());
    LOG_DEBUG("out_xyz: %s, out_ll: %s", ObjLogFormatter<Vec3f>{ ray_out }.Format(),
              ObjLogFormatter<Vec2f>{ out_ll }.Format());

    Vec2f expected_out_ll{ 24.5, -15 };

    EXPECT_NEAR((out_ll - expected_out_ll).norm(), 0, 0.1);
  }

  {
    LOG_INFO("case 2...");
    Vec3f llr{ -63.7804477035974, -2.17764029100887, 209.573936971985 };
    Quatf q = Llr2Quat(llr);

    auto [ray_out, jac] = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);
    Vec2f out_ll = Xyz2Ll(ray_out);
    LOG_DEBUG("llr: %s, quat: %s", ObjLogFormatter<Vec3f>{ llr }.Format(), ObjLogFormatter<Quatf>{ q }.Format());
    LOG_DEBUG("out_xyz: %s, out_ll: %s", ObjLogFormatter<Vec3f>{ ray_out }.Format(),
              ObjLogFormatter<Vec2f>{ out_ll }.Format());

    Vec2f expected_out_ll{ -0.4, 9 };

    EXPECT_NEAR((out_ll - expected_out_ll).norm(), 0, 0.1);
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

  {
    LOG_INFO("test case 1...");
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
      Vec4f res;
      for (const auto& q : c) {
        std::tie(res, std::ignore) = optics_system(q);
        auto ll = Xyz2Ll(res.head(3));
        EXPECT_NEAR((ll - target_ll).norm(), 0, 5e-4);
        LOG_DEBUG("q: %s, ll: %s", ObjLogFormatter<Vec4f>{ q }.Format(), ObjLogFormatter<Vec<float, 2>>{ ll }.Format());
      }
    }
  }

  {
    LOG_INFO("test case 2...");
    Vec2f target_ll{ 23.2, -14.4 };
    auto t0 = std::chrono::high_resolution_clock::now();
    auto [contours, status] = FindAllCrystalPoses(optics_system, target_ll, config);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    LOG_INFO("time elapsed: %dus = %.3fms", dt.count(), dt.count() / 1000.0);

    LOG_INFO("contours: %zu, func_eval_cnt: %zu", contours.size(), status.func_eval_cnt_);
    for (size_t i = 0; i < contours.size(); i++) {
      const auto& c = contours[i];
      LOG_INFO("contour[%zu].size: %zu", i, c.size());
      Vec4f res;
      for (const auto& q : c) {
        std::tie(res, std::ignore) = optics_system(q);
        auto ll = Xyz2Ll(res.head(3));
        EXPECT_NEAR((ll - target_ll).norm(), 0, 5e-4);
        LOG_DEBUG("q: %s, ll: %s", ObjLogFormatter<Vec4f>{ q }.Format(), ObjLogFormatter<Vec<float, 2>>{ ll }.Format());
      }
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
    { -0.179994, 0.045840, -0.015003, 0.982486 },   //
    { -0.175987, 0.095466, -0.011578, 0.979686 },   //
    { -0.175069, 0.144978, -0.008201, 0.973791 },   //
    { -0.177136, 0.194034, -0.005419, 0.964857 },   //
    { -0.190652, 0.289243, -0.004053, 0.938071 },   //
    { -0.202794, 0.334275, -0.006934, 0.920375 },   //
    { -0.219377, 0.376286, -0.013427, 0.900059 },   //
    { -0.241351, 0.413436, -0.024604, 0.877621 },   //
    { -0.269678, 0.442690, -0.041247, 0.854166 },   //
    { -0.304575, 0.459832, -0.062801, 0.831774 },   //
    { -0.324060, 0.462570, -0.074546, 0.821864 },   //
    { -0.344379, 0.461051, -0.086137, 0.813277 },   //
    { -0.365010, 0.455344, -0.096896, 0.806254 },   //
    { -0.385417, 0.445800, -0.106235, 0.800896 },   //
    { -0.405131, 0.432968, -0.113748, 0.797165 },   //
    { -0.423804, 0.417478, -0.119246, 0.794913 },   //
    { -0.441217, 0.399933, -0.122725, 0.793929 },   //
    { -0.471926, 0.360628, -0.124176, 0.794869 },   //
    { -0.517329, 0.273395, -0.110554, 0.803372 },   //
    { -0.543785, 0.181687, -0.083730, 0.815035 },   //
    { -0.552989, 0.088579, -0.050868, 0.826906 },   //
    { -0.544657, -0.004361, -0.017248, 0.838473 },  //
    { -0.516278, -0.094680, 0.011522, 0.851095 },   //
    { -0.493040, -0.136992, 0.021671, 0.858882 },   //
    { -0.462752, -0.175033, 0.027213, 0.868612 },   //
    { -0.425198, -0.205309, 0.026779, 0.881102 },   //
    { -0.404065, -0.216023, 0.024134, 0.888530 },   //
    { -0.381883, -0.222935, 0.020032, 0.896698 },   //
    { -0.359247, -0.225564, 0.014807, 0.905454 },   //
    { -0.336837, -0.223690, 0.008930, 0.914564 },   //
    { -0.315307, -0.217414, 0.002917, 0.923746 },   //
    { -0.295185, -0.207125, -0.002768, 0.932717 },  //
    { -0.276810, -0.193397, -0.007789, 0.941232 },  //
    { -0.245752, -0.158148, -0.015191, 0.956226 },  //
    { -0.221869, -0.116084, -0.019039, 0.967957 },  //
    { -0.204042, -0.070171, -0.019963, 0.976243 },  //
    { -0.182390, 0.026951, -0.016208, 0.982725 },   //
    { -0.179994, 0.045840, -0.015003, 0.982486 },   //
  };

  auto [weight, data] = ComputPoseWeight(rot_contour, crystal, axis_pdf, ray_in_ll, raypath);
  LOG_INFO("contour weight: %.6f", weight);
}


// NOLINTNEXTLINE
TEST_F(TestOptics, contour_and_weight) {
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
  AxisPdf a(0, 0.5);
  Func<float, 1, 3> axis_pdf = [&a](const Vec3f& llr) { return Vec<float, 1>{ a(llr) }; };

  Vec2f target_ll{ 26, -15.1 };
  auto [contours, status] = FindAllCrystalPoses(optics_system, target_ll, config);

  LOG_INFO("contours: %zu, func_eval_cnt: %zu", contours.size(), status.func_eval_cnt_);
  for (size_t i = 0; i < contours.size(); i++) {
    const auto& c = contours[i];
    LOG_INFO("contour[%zu].size: %zu", i, c.size());
    Vec4f res;
    for (const auto& q : c) {
      std::tie(res, std::ignore) = optics_system(q);
      auto ll = Xyz2Ll(res.head(3));
      EXPECT_NEAR((ll - target_ll).norm(), 0, 5e-4);
      LOG_DEBUG("q: %s, ll: %s", ObjLogFormatter<Vec4f>{ q }.Format(), ObjLogFormatter<Vec<float, 2>>{ ll }.Format());
    }

    auto [weight, data] = ComputPoseWeight(c, crystal, axis_pdf, ray_in_ll, raypath);
    LOG_DEBUG("--- s, axis_pdf, jac_factor, geo_factor, transit_factor: w ---");
    for (const auto& d : data) {
      LOG_DEBUG("%.6f,%.4e,%.4e,%.6f,%.4e:%.6f",  //
                d.s_, d.axis_prob_, d.jac_factor_, d.geo_factor_, d.transit_factor_, d.w_);
    }
    LOG_INFO("weight: %.6e", weight);
  }
}

}  // namespace
