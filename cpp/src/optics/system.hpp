#ifndef OPTICS_SYSTEM_H_
#define OPTICS_SYSTEM_H_

#include "core/types.hpp"
#include "optics/optics.hpp"

namespace halo_pm {

// Forward declaration
struct Crystal;

Vec3f TraceDirection(const Crystal& crystal,                                           // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,                            // May be input variables
                     const std::vector<int>& raypath, float wl = kDefaultWavelength);  // Other parameter


std::tuple<Vec3f, Mat3x4f>                               // (output vector, Jacobian wrt quaternion)
TraceDirDiffQuat(const Crystal& crystal,                 // crystal
                 const Quatf& rot, const Vec2f& ray_ll,  // input quaternion, input vector
                 const std::vector<int>& raypath, float wl = kDefaultWavelength);  // other parameter


struct MappingGridData {
  size_t size_;
  Vec3f in_xyz_;
  std::unique_ptr<Vec3f[]> out_xyz_;
  std::unique_ptr<Quatf[]> rot_;
};

MappingGridData GenerateMappingGridData(const Crystal& crystal, const Vec2f& ray_in_ll, const std::vector<int>& raypath,
                                        int level);

}  // namespace halo_pm

#endif  // OPTICS_SYSTEM_H_
