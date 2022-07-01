#include <chrono>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <vector>

#include "core/types.hpp"
#include "optics/crystal.hpp"
#include "optics/system.hpp"
#include "util/log.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

int main(int argc, char** argv) {
  if (argc != 2) {
    return -1;
  }

  auto close_file = [](std::FILE* file) { std::fclose(file); };
  std::unique_ptr<std::FILE, decltype(close_file)> file_ptr(std::fopen(argv[1], "wb"), close_file);

  auto crystal = MakePrismCrystal(1.0f);
  std::vector<int> raypath{ 3, 5 };
  Vec2f sun_ll{ 180, 15 };
  Vec2f ray_in_ll{ sun_ll.x() + 180, -sun_ll.y() };
  auto config = MakeConfigData(crystal, ray_in_ll, raypath, 3);

  AxisPdf a(90, 0.5);
  Func<float, 1, 3> axis_pdf = [&a](const Vec3f& llr) { return Vec<float, 1>{ a(llr) }; };

  auto optics_system = [&crystal, &ray_in_ll, &raypath](const Vec4f& rot) -> std::tuple<Vec4f, Mat4x4f> {
    Quatf q{ rot(0), rot(1), rot(2), rot(3) };  // w, x, y, z
    auto [xyz, j] = TraceDirDiffQuat(crystal, q, ray_in_ll, raypath);
    Vec4f out;
    Mat4x4f jac;
    out << xyz, rot.squaredNorm();
    jac << j, 2 * rot.transpose();
    return std::make_tuple(out, jac);
  };

  std::vector<float> lon_list;
  std::vector<float> lat_list;
  float dx = 0.02f;
  for (float lon = -2.5f; lon < 2.5f + dx / 2; lon += dx) {
    lon_list.emplace_back(lon);
  }
  for (float lat = 6.0f; lat < 22.0f + dx / 2; lat += dx) {
    lat_list.emplace_back(lat);
  }
  size_t img_wid = lon_list.size();
  size_t img_hei = lat_list.size();
  std::unique_ptr<float[]> img_data_buf{ new float[img_wid * img_hei]{} };
  size_t total_pix = img_wid * img_hei;

  auto t0 = std::chrono::high_resolution_clock::now();
  for (size_t x = 0; x < img_wid; x++) {
    auto lon = lon_list[x];
    LOG_INFO("progress %.2f%%", x * img_hei * 100.0 / total_pix);
    for (size_t y = 0; y < img_hei; y++) {
      auto lat = lat_list[y];
      Vec2f target_ll{ lon, lat };
      auto [contours, status] = FindAllCrystalPoses(optics_system, target_ll, config);
      float weight = 0;
      for (const auto& c : contours) {
        auto [w, data] = ComputPoseWeight(c, crystal, axis_pdf, ray_in_ll, raypath);
        weight += w;
      }
      img_data_buf[y * img_wid + x] = weight;
      LOG_DEBUG("[%.3f,%.3f]:%.6e", lon, lat, weight);
    }
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
  std::fwrite(img_data_buf.get(), sizeof(float), total_pix, file_ptr.get());
  LOG_INFO("time elapsed: %dms = %.3fs", dt.count(), dt.count() / 1000.0);

  return 0;
}