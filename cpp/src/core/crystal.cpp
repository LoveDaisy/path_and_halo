#include "core/crystal.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "core/geo.hpp"
#include "core/math.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

namespace halo_pm {

Crystal MakePrismCrystal(float h) {
  constexpr size_t kVtxCnt = 12;

  std::unique_ptr<float[]> vtx{ new float[kVtxCnt * 3]{} };

  std::vector<std::vector<int>> face_id{
    { 0, 1, 2, 3, 4, 5 },    // fn1
    { 11, 10, 9, 8, 7, 6 },  // fn2
    { 0, 6, 7, 1 },          // fn3
    { 1, 7, 8, 2 },          // fn4
    { 2, 8, 9, 3 },          // fn5
    { 3, 9, 10, 4 },         // fn6
    { 4, 10, 11, 5 },        // fn7
    { 5, 11, 12, 6 },        // fn8
  };
  size_t face_cnt = face_id.size();
  std::unique_ptr<float[]> face_norm{ new float[face_cnt * 3]{} };
  std::unique_ptr<float[]> face_area{ new float[face_cnt * 1]{} };

  for (int i = 0; i < 6; i++) {
    auto q = i * kPi / 3 - kPi / 6;
    vtx[i * 3 + 0] = 0.5 * std::cos(q);
    vtx[i * 3 + 1] = 0.5 * std::sin(q);
    vtx[i * 3 + 2] = 0.5 * h;
    vtx[i * 3 + 18] = 0.5 * std::cos(q);
    vtx[i * 3 + 19] = 0.5 * std::sin(q);
    vtx[i * 3 + 20] = -0.5 * h;
  }

  for (size_t i = 0; i < face_cnt; i++) {
    Eigen::Map<Vec3f> vtx0(vtx.get() + face_id[i][0] * 3);
    Eigen::Map<Vec3f> vtx1(vtx.get() + face_id[i][1] * 3);
    Eigen::Map<Vec3f> vtx2(vtx.get() + face_id[i][2] * 3);

    Eigen::Map<Vec3f> n(face_norm.get() + i * 3);
    n = (vtx1 - vtx0).cross(vtx2 - vtx0);
    n.normalize();
  }

  for (size_t i = 0; i < face_cnt; i++) {
    const auto& id = face_id[i];
    Eigen::Map<Vec3f> vtx0(vtx.get() + id[0] * 3);
    for (size_t j = 1; j + 1 < id.size(); j++) {
      Eigen::Map<Vec3f> vtx1(vtx.get() + id[j] * 3);
      Eigen::Map<Vec3f> vtx2(vtx.get() + id[j + 1] * 3);

      auto n = (vtx1 - vtx0).cross(vtx2 - vtx0);
      face_area[i] += n.norm() / 2.0f;
    }
  }

  return { kVtxCnt, face_cnt, std::move(vtx), std::move(face_id), std::move(face_norm), std::move(face_area) };
}


constexpr float GetIceRefractiveIndex(float lambda) {
  if (lambda < kMinWavelength || lambda > kMaxWavelength) {
    return -1.0f;
  }

  lambda /= 1000.0f;  // Convert to micrometer (um)
  float lambda2 = lambda * lambda;

  float n2 = 1.0f;
  n2 += 0.701777f / (1.0f - 0.884400e-2 / lambda2);
  n2 += 1.091144 / (1.0f - 0.796950e2 / lambda2);
  return std::sqrt(n2);
}


Vec3f Refract(const Vec3f& ray_in, const Vec3f& norm, float n0, float n1) {
  auto n = n0 / n1;
  auto c = ray_in.dot(norm);
  auto delta = n * n - (n * n - 1) / (c * c);

  if (delta <= 0) {
    return (Mat3x3f::Identity() - 2 * norm * norm.transpose()) * ray_in;
  } else {
    return n * ray_in + (std::sqrt(delta) - n) * c * norm;
  }
}


Vec3f TraceDirection(const Crystal& crystal,                       // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,        // May be input variables
                     const std::vector<int>& raypath, float wl) {  // Other parameter
  if (raypath.empty()) {
    return Vec3f{};
  }

  auto xyz = Ll2Xyz(ray_ll);
  LOG_DEBUG("ll2xyz: [%.6f,%.6f,%.6f]", xyz.x(), xyz.y(), xyz.z());

  // Only ONE face: pure reflection
  if (raypath.size() == 1) {
    const auto* n_ptr = crystal.face_norm_.get() + (raypath[0] - 1) * 3;
    Eigen::Map<const Vec3f> n(n_ptr);
    auto m = Mat3x3f::Identity() - 2 * n * n.transpose();
    return rot * m * rot.inverse() * xyz;
  }

  // Multiple faces: reflection + refraction
  else {
    Mat3x3f m = Eigen::Matrix3f::Identity();
    for (size_t i = 1; i + 1 < raypath.size(); i++) {
      Eigen::Map<const Vec3f> n(crystal.face_norm_.get() + (raypath[i] - 1) * 3);
      m = (Mat3x3f::Identity() - 2 * n * n.transpose()) * m;
    }

    auto refractive_index = GetIceRefractiveIndex(wl);
    Eigen::Map<const Vec3f> norm1(crystal.face_norm_.get() + (raypath.front() - 1) * 3);
    Eigen::Map<const Vec3f> norm2(crystal.face_norm_.get() + (raypath.back() - 1) * 3);
    auto r = rot.inverse() * xyz;
    LOG_DEBUG("r: [%.6f,%.6f,%.6f]", r.x(), r.y(), r.z());
    r = Refract(r, norm1, 1.0f, refractive_index);
    LOG_DEBUG("r: [%.6f,%.6f,%.6f]", r.x(), r.y(), r.z());
    r = m * r;
    LOG_DEBUG("r: [%.6f,%.6f,%.6f]", r.x(), r.y(), r.z());
    r = Refract(r, norm2, refractive_index, 1.0f);
    LOG_DEBUG("r: [%.6f,%.6f,%.6f]", r.x(), r.y(), r.z());
    return rot * r;
  }
}

}  // namespace halo_pm
