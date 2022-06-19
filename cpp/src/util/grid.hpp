#ifndef SPH_GRID_H_
#define SPH_GRID_H_

#include <cstddef>
#include <memory>
#include <tuple>

#include "core/types.hpp"

namespace halo_pm {

struct SphGrid {
  size_t n_pix_;  // pix number
  float dr_;      // in degree
  std::unique_ptr<Vec3f[]> v_;
};


SphGrid GenerateGrid(int level);

}  // namespace halo_pm

#endif  // SPH_GRID_H_
