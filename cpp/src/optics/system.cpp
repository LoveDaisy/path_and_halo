#include "optics/system.hpp"

#include <algorithm>
#include <tuple>
#include <vector>

#include "core/math.hpp"
#include "core/solver.hpp"
#include "core/types.hpp"
#include "geo/geo.hpp"
#include "geo/grid.hpp"
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


std::tuple<std::vector<Curve4f>, PoseContourStatus>  // (all contours, status)
FindAllPoseContour(const FuncAndDiff<float, 4, 4>& optics_system, const Vec2f& target_ll, const ConfigData& config) {
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

  PoseContourStatus status;
  std::vector<Curve4f> res_contours;

  auto seed_threshold = config.dr * 1.5f;
  for (const auto& [q, a] : sorted_rot_da) {
    if (a < seed_threshold) {
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
      float d = DistanceToPolyLine(rot0, c);
      if (d < kReduceTh) {
        reduced = true;
        break;
      }
      d = DistanceToPolyLine(Vec4f(-rot0), c);
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
    // NOTE: is it neccessary to use Jacobian? Maybe directly calculate difference in output space.
    Curve4f tmp_rot_cand;
    for (const auto& q : rot_cand) {
      Mat4x4f jac;
      auto info = DistanceToPolyLine(q, contour);
      std::tie(std::ignore, jac) = optics_system(info.nearest_point_);
      status.func_eval_cnt_++;
      Vec4f dx = q - info.nearest_point_;
      auto dy = (jac * dx).norm();
      if (dy < config.dr * 2.5) {
        continue;
      }

      info = DistanceToPolyLine(Vec4f(-q), contour);
      std::tie(std::ignore, jac) = optics_system(info.nearest_point_);
      status.func_eval_cnt_++;
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

}  // namespace halo_pm
