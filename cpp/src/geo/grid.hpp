#ifndef GEO_GRID_H_
#define GEO_GRID_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "core/types.hpp"

namespace halo_pm {

struct SphGrid {
  float dr_;  // in degree
  std::vector<Vec3f> xyz_;
};


SphGrid GenerateGrid(int level);

}  // namespace halo_pm

#endif  // GEO_GRID_H_
