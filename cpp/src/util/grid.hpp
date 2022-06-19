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


struct MappingGridData {
  size_t size_;
  Vec3f in_xyz_;
  std::unique_ptr<Vec3f[]> out_xyz_;
  std::unique_ptr<Quatf[]> rot_;
};

// Forward declaration
struct Crystal;

MappingGridData GenerateMappingGridData(const Crystal& crystal, const Vec2f& ray_in_ll, const std::vector<int>& raypath,
                                        int level);

}  // namespace halo_pm

#endif  // SPH_GRID_H_
