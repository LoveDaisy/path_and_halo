#ifndef SPH_GRID_H_
#define SPH_GRID_H_

#include <cstddef>
#include <memory>
#include <tuple>

namespace halo_pm {

std::tuple<std::unique_ptr<float[]>, size_t> GenerateGrid(int level);

}  // namespace halo_pm

#endif  // SPH_GRID_H_
