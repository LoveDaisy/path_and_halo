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

float GetIceRefractiveIndex(float lambda);

}  // namespace halo_pm

#endif  // CORE_CRYSTAL_H_
