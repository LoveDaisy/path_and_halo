#ifndef OPTICS_SYSTEM_H_
#define OPTICS_SYSTEM_H_

#include <cstddef>
#include <tuple>
#include <vector>

#include "core/types.hpp"
#include "optics/crystal.hpp"
#include "optics/optics.hpp"

namespace halo_pm {

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
struct FindPoseStatus {
  size_t func_eval_cnt_;
  std::vector<Quatf> rot_seeds_;  // quaternion [wxyz]

  FindPoseStatus() : func_eval_cnt_(0) {}
};

std::tuple<std::vector<Curve4f>, FindPoseStatus>                    // (all contours of rotation, status)
FindAllCrystalPoses(const FuncAndDiff<float, 4, 4>& optics_system,  // quaternion -> (out_xyz, q_norm2)
                    const Vec2f& target_ll, const ConfigData& config);


struct PoseWeightData {
  float s_;               // curve arc length parameter
  float w_;               // total weight
  float axis_prob_;       // probability of axis orientation
  float jac_factor_;      // determinator of Jacobian
  float geo_factor_;      // geometry factor
  float transit_factor_;  // transit factor
  Vec3f axis_llr_;        // axis pose
};


PoseWeightData  // one data sample
ComputeWeightComponents(const Vec4f& rot, const Crystal& crystal, const Func<float, 1, 3>& axis_pdf,
                        const Vec2f& ray_in_ll, const std::vector<int>& raypath);

std::tuple<float, std::vector<PoseWeightData>>  // (total weight, interpolated weight components)
ComputPoseWeight(const Curve4f& rots, const Crystal& crystal, const Func<float, 1, 3>& axis_pdf, const Vec2f& ray_in_ll,
                 const std::vector<int>& raypath);

}  // namespace halo_pm

#endif  // OPTICS_SYSTEM_H_
