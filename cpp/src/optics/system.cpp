#include "optics/system.hpp"

#include "geo/grid.hpp"

namespace halo_pm {

MappingGridData GenerateMappingGridData(const Crystal& crystal, const Vec2f& ray_in_ll, const std::vector<int>& raypath,
                                        int level) {
  auto [n_pix, dr, ll] = GenerateGrid(level);
  if (n_pix == 0) {
    return {};
  }

  auto roll_cnt = static_cast<size_t>(360 / dr / 2);
  std::unique_ptr<float[]> roll_list{ new float[roll_cnt]{} };
  for (size_t i = 0; i < roll_cnt; i++) {
    roll_list[i] = dr * i * 2;
  }

  std::unique_ptr<Quatf[]> rot{ new Quatf[n_pix * roll_cnt]{} };
  for (size_t i = 0; i < n_pix; i++) {
    for (size_t j = 0; j < roll_cnt; j++) {
      rot[i * roll_cnt + j] = Llr2Quat(Vec3f{ ll[i].x(), ll[i].y(), roll_list[j] });
    }
  }

  std::unique_ptr<Vec3f[]> ray_out_xyz{ new Vec3f[n_pix * roll_cnt]{} };
  for (size_t i = 0; i < n_pix * roll_cnt; i++) {
    ray_out_xyz[i] = TraceDirection(crystal, rot[i], ray_in_ll, raypath);
  }

  return MappingGridData{ n_pix * roll_cnt, Ll2Xyz(ray_in_ll), std::move(ray_out_xyz), std::move(rot) };
}

}  // namespace halo_pm
