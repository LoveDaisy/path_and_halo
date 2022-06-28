#include "optics/system.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

#include "core/math.hpp"
#include "core/solver.hpp"
#include "core/types.hpp"
#include "geo/geo.hpp"
#include "geo/grid.hpp"
#include "geo/spline.hpp"
#include "optics/crystal.hpp"
#include "optics/optics.hpp"
#include "util/log.hpp"

namespace halo_pm {

ConfigData MakeConfigData(const Crystal& crystal, const Vec2f& ray_in_ll, const std::vector<int>& raypath, int level) {
  auto grid = GenerateGrid(level);
  if (grid.xyz_.empty()) {
    return {};
  }

  const auto& xyz = grid.xyz_;
  auto dr = grid.dr_;
  auto n_pix = xyz.size();

  auto roll_cnt = static_cast<size_t>(360 / dr / 2);
  std::unique_ptr<float[]> roll_list{ new float[roll_cnt]{} };
  for (size_t i = 0; i < roll_cnt; i++) {
    roll_list[i] = dr * i * 2;
  }

  std::vector<std::tuple<Quatf, Vec3f>> res;
  res.reserve(n_pix * roll_cnt);
  for (size_t i = 0; i < n_pix; i++) {
    Vec2f ll = Xyz2Ll(xyz[i]);
    for (size_t j = 0; j < roll_cnt; j++) {
      auto curr_rot = Llr2Quat(Vec3f{ ll.x(), ll.y(), roll_list[j] });
      auto curr_xyz = TraceDirection(crystal, curr_rot, ray_in_ll, raypath);
      res.emplace_back(std::make_tuple(curr_rot, curr_xyz));
    }
  }

  return ConfigData{ dr, Ll2Xyz(ray_in_ll), std::move(res) };
}


// =============== Find pose contour ===============
std::tuple<std::vector<Curve4f>, FindPoseStatus>  // (all contours, status)
FindAllCrystalPoses(const FuncAndDiff<float, 4, 4>& optics_system, const Vec2f& target_ll, const ConfigData& config) {
  // Find candidate rotations (rotation seeds) from config data
  Vec3f target_xyz = Ll2Xyz(target_ll);
  std::vector<std::tuple<Quatf, float>> sorted_rot_da;
  for (const auto& [q, out_xyz] : config.rot_xyz_) {
    auto da = std::acos(out_xyz.dot(target_xyz)) * kRad2Degree;
    sorted_rot_da.emplace_back(std::make_tuple(q, da));
  }
  std::sort(sorted_rot_da.begin(), sorted_rot_da.end(), [](const auto& x1, const auto& x2) {
    auto a1 = std::get<1>(x1);
    auto a2 = std::get<1>(x2);
    return a1 < a2;
  });

  FindPoseStatus status;
  std::vector<Curve4f> res_contours;

  auto seed_threshold = config.dr * 1.5f;
  for (const auto& [q, a] : sorted_rot_da) {
    if (a < seed_threshold) {
      // NaNs will not included, because NaN < limited_value always be false.
      status.rot_seeds_.emplace_back(q);
    }
  }
  LOG_DEBUG("find %zu seeds", status.rot_seeds_.size());

  if (status.rot_seeds_.empty()) {
    return { res_contours, status };
  }

  Curve4f rot_cand;
  for (const auto& q : status.rot_seeds_) {
    Vec4f qv{ q.w(), q.x(), q.y(), q.z() };
    rot_cand.emplace_back(qv);
  }
  std::reverse(rot_cand.begin(), rot_cand.end());  // from largest distance to smallest distance

  SolverOption solver_option;
  solver_option.h_ = 0.05;
  constexpr float kReduceTh = 0.05;
  Vec4f target_out;
  target_out << target_xyz, 1.0f;
  while (!rot_cand.empty()) {
    // Find initial solution
    auto [rot0, rot0_status] = FindSolution(optics_system, rot_cand.back(), target_out);
    rot_cand.pop_back();
    status.func_eval_cnt_ += rot0_status.func_eval_cnt_;
    LOG_DEBUG("find solution: %zu func_eval_cnt", rot0_status.func_eval_cnt_);
    if (!rot0_status.solved_) {
      continue;
    }

    // Check if it locates on/close to previous contour
    bool reduced = false;
    for (const auto& c : res_contours) {
      float d = 0;
      std::tie(d, std::ignore) = DistanceToPolyLine(rot0, c);
      if (d < kReduceTh) {
        reduced = true;
        break;
      }
      std::tie(d, std::ignore) = DistanceToPolyLine(Vec4f(-rot0), c);
      if (d < kReduceTh) {
        reduced = true;
        break;
      }
    }
    if (reduced) {
      continue;
    }

    // Find a contour
    auto [contour, contour_status] = FindContour(optics_system, rot0, solver_option);
    status.func_eval_cnt_ += contour_status.func_eval_cnt_;
    LOG_DEBUG("find contour: %zu func_eval_cnt", contour_status.func_eval_cnt_);

    if (contour.empty()) {
      continue;
    } else {
      res_contours.emplace_back(contour);
    }

    // Reduce seeds
    // NOTE: is it neccessary to use Jacobian? Maybe we can directly calculate difference in output space.
    std::vector<Mat4x4f> contour_jac;
    for (const auto& q : contour) {
      Mat4x4f jac;
      std::tie(std::ignore, jac) = optics_system(q);
      status.func_eval_cnt_++;
      contour_jac.emplace_back(jac);
    }
    Curve4f tmp_rot_cand;
    for (const auto& q : rot_cand) {
      PointLineDistanceInfo<float, 4> info{};
      size_t idx = 0;
      std::tie(info, idx) = DistanceToPolyLine(q, contour);
      Mat4x4f jac = contour_jac[idx] * (1 - info.t_) + contour_jac[idx + 1] * info.t_;
      Vec4f dx = q - info.nearest_point_;
      auto dy = (jac * dx).norm();
      if (dy < config.dr * 2.5) {
        continue;
      }

      std::tie(info, idx) = DistanceToPolyLine(Vec4f(-q), contour);
      jac = contour_jac[idx] * (1 - info.t_) + contour_jac[idx + 1] * info.t_;
      dx = -q - info.nearest_point_;
      dy = (jac * dx).norm();
      if (dy < config.dr * 2.5) {
        continue;
      }

      tmp_rot_cand.emplace_back(q);
    }
    rot_cand.swap(tmp_rot_cand);
  }

  return { res_contours, status };
}


// =============== Compute weight ===============
float ComputeJacobianFactor(const Vec4f& rot, const Crystal& crystal, const Vec2f& ray_in_ll,
                            const std::vector<int>& raypath) {
  Quatf q{ rot(0), rot(1), rot(2), rot(3) };
  q.normalize();
  Mat3x4f jac;
  std::tie(std::ignore, jac) = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);

  constexpr float kMinVal = 1e-8;
  Eigen::JacobiSVD svd(jac);
  svd.setThreshold(kDefaultFloatEps);

  float res = 0.0f;
  auto rank = svd.rank();
  if (rank > 0) {
    res = 1.0f;
    const auto& s = svd.singularValues();
    for (int i = 0; i < rank; i++) {
      res *= s(i);
    }
  }
  return 1.0f / std::max(res, kMinVal);
}

float ComputeEntryFactor(const Vec4f& rot, const Crystal& crystal, const Vec2f& ray_in_ll,
                         const std::vector<int>& raypath) {
  Quatf q{ rot(0), rot(1), rot(2), rot(3) };
  Vec3f ray_in_xyz = q.inverse() * Ll2Xyz(ray_in_ll);

  auto entry_fid = raypath.front() - 1;
  const auto* norm_ptr = crystal.face_norm_.get() + entry_fid * 3;
  Eigen::Map<const Vec3f> entry_norm(norm_ptr);

  float entry_area = -1 * crystal.face_area_[entry_fid] * (entry_norm.dot(ray_in_xyz));
  if (entry_area < 0) {
    return 0;
  } else {
    float total_area = 0.0f;
    for (size_t i = 0; i < crystal.face_cnt_; i++) {
      norm_ptr = crystal.face_norm_.get() + i * 3;
      Eigen::Map<const Vec3f> curr_norm(norm_ptr);
      auto curr_area = std::max(-1 * crystal.face_area_[i] * curr_norm.dot(ray_in_xyz), 0.0f);
      total_area += curr_area;
    }
    return entry_area / total_area;
  }
}

std::tuple<Vec3f, float, float>  // ray_out, transmissivity, cos_factor
TransitFace(const Vec3f& ray, const Vec3f& norm, float n0, float n1) {
  Vec3f ray_out;
  if (n0 * n1 < 0) {
    auto m = Mat3x3f::Identity() - 2 * norm * norm.transpose();
    ray_out = m * ray;
  } else {
    ray_out = Refract(ray, norm, n0, n1);
  }

  n0 = std::abs(n0);
  n1 = std::abs(n1);
  auto cos_qi = std::abs(ray.dot(norm));
  auto cos_qt = std::abs(ray_out.dot(norm));
  auto Rs1 = std::abs((n0 * cos_qi - n1 * cos_qt) / (n0 * cos_qi + n1 * cos_qt));
  Rs1 *= Rs1;
  auto Rp1 = std::abs((n0 * cos_qt - n1 * cos_qi) / (n0 * cos_qt + n1 * cos_qi));
  Rp1 *= Rp1;
  auto t = (1 - (Rs1 + Rp1) / 2);
  auto q_factor = cos_qi / cos_qt;

  return { ray_out, t, q_factor };
}

std::tuple<float, float>  // transit_factor, geo_factor
ComputeTransitGeoFactor(const Vec4f& rot, const Crystal& crystal, const Vec2f& ray_in_ll,
                        const std::vector<int>& raypath) {
  Quatf rot_q{ rot(0), rot(1), rot(2), rot(3) };
  Vec3f r = Ll2Xyz(ray_in_ll);
  r = rot_q.inverse() * r;  // convert to crystal-local frame

  // generate refractive index list
  std::vector<float> rp_n(raypath.size() + 1, 1.0f);
  auto crystal_n = GetIceRefractiveIndex(kDefaultWavelength);
  for (size_t i = 1; i + 1 < rp_n.size(); i++) {
    rp_n[i] = crystal_n;  // TODO:
  }
  for (size_t i = 2; i + 1 < rp_n.size(); i++) {
    for (size_t j = i; j < rp_n.size(); j++) {
      rp_n[j] *= -1;
    }
  }

  float transit_factor = 1.0f;
  float geo_factor = 1.0f;
  Curve3f poly0;
  for (auto idx : crystal.face_id_[raypath.front() - 1]) {
    poly0.emplace_back(Eigen::Map<const Vec3f>(crystal.vtx_.get() + idx * 3));
  }

  for (size_t i = 0; i < raypath.size(); i++) {
    auto curr_fid = raypath[i] - 1;
    auto n0 = rp_n[i];
    auto n1 = rp_n[i + 1];
    const auto* norm_ptr = crystal.face_norm_.get() + curr_fid * 3;
    Eigen::Map<const Vec3f> curr_norm(norm_ptr);

    float t = 0;
    float q = 0;

    // Refraction
    if (n0 * n1 > 0) {
      std::tie(r, t, q) = TransitFace(r, curr_norm, n0, n1);
    }

    // Reflection in crystal
    else if (std::abs(n0) > 1 + kDefaultFloatEps) {
      std::tie(std::ignore, t, std::ignore) = TransitFace(r, curr_norm, std::abs(n0), 1.0f);
      if (std::isnan(t)) {
        t = 1;
      } else {
        t = 1 - t;
      }
      r = (Mat3x3f::Identity() - 2 * curr_norm * curr_norm.transpose()) * r;
      q = 1;
    }

    // Reflection out crystal
    else {
      std::tie(std::ignore, t, std::ignore) = TransitFace(r, curr_norm, 1.0f, crystal_n);
      t = 1 - t;
      q = 1;
      r = (Mat3x3f::Identity() - 2 * curr_norm * curr_norm.transpose()) * r;
    }

    transit_factor *= (t * q);

    if (i + 1 < raypath.size() && geo_factor > 1e-8) {
      Curve3f poly1;
      auto next_fid = raypath[i + 1] - 1;
      for (auto idx : crystal.face_id_[next_fid]) {
        poly1.emplace_back(Eigen::Map<const Vec3f>(crystal.vtx_.get() + idx * 3));
      }
      float a = 0;
      std::tie(poly0, a) = ProjectPolygonIntersection(poly0, poly1, r);
      geo_factor *= a;
    }
  }

  if (std::isnan(transit_factor)) {
    transit_factor = 0;
  }
  if (std::isnan(geo_factor)) {
    geo_factor = 0;
  }
  return { transit_factor, geo_factor };
}

PoseWeightData  // one data sample
ComputeWeightComponents(const Vec4f& rot, const Crystal& crystal, const Func<float, 1, 3>& axis_pdf,
                        const Vec2f& ray_in_ll, const std::vector<int>& raypath) {
  PoseWeightData data{};
  data.jac_factor_ = ComputeJacobianFactor(rot, crystal, ray_in_ll, raypath);
  auto entry_factor = ComputeEntryFactor(rot, crystal, ray_in_ll, raypath);
  std::tie(data.transit_factor_, data.geo_factor_) = ComputeTransitGeoFactor(rot, crystal, ray_in_ll, raypath);
  data.geo_factor_ *= entry_factor;
  data.axis_llr_ = Quat2Llr(rot);
  data.axis_prob_ = axis_pdf(data.axis_llr_).x();
  data.w_ = data.axis_prob_ * data.jac_factor_ * data.geo_factor_ * data.transit_factor_;
  data.s_ = 0;
  return data;
}

std::vector<PoseWeightData>  // only the components, wrt to input rotations
ComputeWeightComponents(const Curve4f& rots, const Crystal& crystal, const Func<float, 1, 3>& axis_pdf,
                        const Vec2f& ray_in_ll, const std::vector<int>& raypath) {
  std::vector<PoseWeightData> res;
  res.reserve(rots.size());

  for (const auto& r : rots) {
    res.emplace_back(ComputeWeightComponents(r, crystal, axis_pdf, ray_in_ll, raypath));
  }

  return res;
}

std::tuple<float, std::vector<PoseWeightData>>  // (total weight, interpolated weight components)
ComputPoseWeight(const Curve4f& rots, const Crystal& crystal, const Func<float, 1, 3>& axis_pdf, const Vec2f& ray_in_ll,
                 const std::vector<int>& raypath) {
  auto origin_data_num = rots.size();
  auto arc_len = CumulatedArcLength(rots);
  auto interp_step = arc_len.back() * 0.01f;
  constexpr float kPdfTh = 1e-10;

  // Original components
  auto cmp0 = ComputeWeightComponents(rots, crystal, axis_pdf, ray_in_ll, raypath);

  // Check and make sure there are enough points whthin high PDF area
  Curve4f interp_rot;
  std::vector<float> interp_pdf;
  std::vector<float> interp_s;
  std::vector<float> s0;
  while (true) {
    std::tie(interp_rot, interp_s, s0) = InterpCurve(rots, interp_step);

    float max_pdf = std::numeric_limits<float>::lowest();
    float high_pdf_cnt = 0;
    interp_pdf.clear();
    for (const auto& r : interp_rot) {
      auto p = axis_pdf(Quat2Llr(r));
      interp_pdf.emplace_back(p(0));
      if (p(0) > max_pdf) {
        max_pdf = p(0);
      }
      if (p(0) >= kPdfTh) {
        high_pdf_cnt++;
      }
    }

    if (max_pdf > kPdfTh && high_pdf_cnt < 50 && interp_step > arc_len.back() * 0.001) {
      interp_step *= 0.5;
    } else {
      break;
    }
  }

  auto interp_data_num = interp_rot.size();
  std::vector<PoseWeightData> w_cmp(interp_data_num);
  bool refine = false;
  for (size_t i = 0, j = 0; i < interp_data_num; i++) {
    auto& curr_cmp = w_cmp[i];
    while (j + 1 < origin_data_num && s0[j + 1] < interp_s[i]) {
      j++;
    }

    if (!refine && (((i == 0 || i + 1 == interp_data_num) && interp_pdf[i] > kPdfTh) ||
                    (i + 1 < interp_data_num && interp_pdf[i] < kPdfTh && interp_pdf[i + 1] > kPdfTh))) {
      refine = true;
    }

    // If no refine needed
    if (!refine) {
      // interp jac_factor
      float y = 0.0f;
      float y0 = std::log(cmp0[j].jac_factor_);
      float y1 = std::log(cmp0[j + 1].jac_factor_);
      y = InterpLinear(s0[j], s0[j + 1], y0, y1, interp_s[i]);
      curr_cmp.jac_factor_ = std::exp(y);

      // interp geo_factor
      y0 = cmp0[j].geo_factor_;
      y1 = cmp0[j + 1].geo_factor_;
      y = InterpLinear(s0[j], s0[j + 1], y0, y1, interp_s[i]);
      curr_cmp.geo_factor_ = y;

      // interp transit_factor
      y0 = std::log(cmp0[j].transit_factor_);
      y1 = std::log(cmp0[j + 1].transit_factor_);
      y = InterpLinear(s0[j], s0[j + 1], y0, y1, interp_s[i]);
      curr_cmp.transit_factor_ = std::exp(y);
    }
    // Need refine
    else {
      curr_cmp = ComputeWeightComponents(interp_rot[i], crystal, axis_pdf, ray_in_ll, raypath);
    }

    if (refine && i > 0 && interp_pdf[i - 1] > kPdfTh && interp_pdf[i] < kPdfTh) {
      refine = false;
    }
    curr_cmp.axis_prob_ = interp_pdf[i];
    curr_cmp.s_ = interp_s[i];
    if (std::isnan(curr_cmp.jac_factor_)) {
      curr_cmp.jac_factor_ = 0.0f;
    }
    if (std::isnan(curr_cmp.geo_factor_)) {
      curr_cmp.geo_factor_ = 0.0f;
    }
    if (std::isnan(curr_cmp.transit_factor_)) {
      curr_cmp.transit_factor_ = 0.0f;
    }
    curr_cmp.w_ = curr_cmp.axis_prob_ * curr_cmp.jac_factor_ * curr_cmp.geo_factor_ * curr_cmp.transit_factor_;
    curr_cmp.w_ = std::max(curr_cmp.w_, 0.0f);
  }

  float total_weight = 0;
  for (size_t i = 0; i + 1 < interp_data_num; i++) {
    total_weight += (w_cmp[i].w_ + w_cmp[i + 1].w_) * (w_cmp[i + 1].s_ - w_cmp[i].s_) * 0.5;
  }
  return { total_weight, w_cmp };
}

}  // namespace halo_pm
