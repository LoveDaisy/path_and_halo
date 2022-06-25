#ifndef OPTICS_SYSTEM_H_
#define OPTICS_SYSTEM_H_

#include <cstddef>
#include <tuple>
#include <vector>

#include "core/types.hpp"
#include "optics/optics.hpp"

namespace halo_pm {

// Forward declaration
struct Crystal;

// =============== Pre-computed cache ===============
struct ConfigData {
  float dr;
  Vec3f in_xyz_;
  std::vector<std::tuple<Quatf, Vec3f>> rot_xyz_;
};

ConfigData MakeConfigData(const Crystal& crystal, const Vec2f& ray_in_ll, const std::vector<int>& raypath, int level);


// =============== Trace throgh the crystal optics system ===============
Vec3f TraceDirection(const Crystal& crystal,                                           // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,                            // May be input variables
                     const std::vector<int>& raypath, float wl = kDefaultWavelength);  // Other parameter


std::tuple<Vec3f, Mat3x4f>                               // (output vector, Jacobian wrt quaternion)
TraceDirDiffQuat(const Crystal& crystal,                 // crystal
                 const Quatf& rot, const Vec2f& ray_ll,  // input quaternion, input vector
                 const std::vector<int>& raypath, float wl = kDefaultWavelength);  // other parameter


// =============== Find all pose contours ===============
struct PoseContourStatus {
  size_t func_eval_cnt_;
  std::vector<Quatf> rot_seeds_;  // quaternion [wxyz]

  PoseContourStatus() : func_eval_cnt_(0) {}
};

std::tuple<std::vector<Curve4f>, PoseContourStatus>  // (all contours of rotation, status)
FindAllPoseContour(const FuncAndDiff<float, 4, 4>& optics_system, const Vec2f& target_ll, const ConfigData& config);

}  // namespace halo_pm

#endif  // OPTICS_SYSTEM_H_
