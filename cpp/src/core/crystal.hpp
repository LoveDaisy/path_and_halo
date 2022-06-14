#ifndef CORE_CRYSTAL_H_
#define CORE_CRYSTAL_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "core/types.hpp"

namespace halo_pm {

struct Crystal {
  size_t vtx_cnt_;
  size_t face_cnt_;
  std::unique_ptr<float[]> vtx_;
  std::vector<std::vector<int>> face_id_;
  std::unique_ptr<float[]> face_norm_;
  std::unique_ptr<float[]> face_area_;
};


Crystal MakePrismCrystal(float h);

constexpr float kMinWavelength = 350.0f;
constexpr float kMaxWavelength = 850.0f;
constexpr float kDefaultWavelength = 546.1f;  // e-line

constexpr inline float GetIceRefractiveIndex(float lambda);


Vec3f TraceDirection(const Crystal& crystal,                                           // Crystal
                     const Quatf& rot, const Vec2f& ray_ll,                            // May be input variables
                     const std::vector<int>& raypath, float wl = kDefaultWavelength);  // Other parameter

}  // namespace halo_pm

#endif  // CORE_CRYSTAL_H_
