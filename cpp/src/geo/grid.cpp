#include "geo/grid.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <tuple>

#include "core/math.hpp"
#include "core/types.hpp"
#include "geo/geo.hpp"
#include "optics/optics.hpp"
#include "util/log.hpp"

namespace halo_pm {

void FillVec(double z, double phi, Vec3f* xyz) {
  auto sth = std::sqrt((1 - z * z));
  auto x = sth * std::cos(phi);
  auto y = sth * std::sin(phi);

  xyz->x() = x;
  xyz->y() = y;
  xyz->z() = z;
}

SphGrid GenerateGrid(int level) {
  if (level <= 0) {
    return SphGrid{};
  }

  size_t n_side = 1 << level;
  size_t n_pix = 12 * n_side * n_side;
  std::vector<Vec3f> data(n_pix);


  size_t n_cap = 2 * n_side * (n_side - 1);
  size_t n_eq = 2 * n_side * (5 * n_side + 1);
  auto fact1 = 1.5 * n_side;
  auto fact2 = 3.0 * n_side * n_side;
  size_t nl2 = 2 * n_side;
  size_t nl4 = 4 * n_side;

  // North cap
  for (size_t i = 0; i < n_cap; i++) {
    auto hip = (i + 1) / 2.0;
    auto fihip = std::trunc(hip);
    auto i_ring = std::trunc(std::sqrt(hip - std::sqrt(fihip))) + 1;
    auto i_phi = (i + 1) - 2 * i_ring * (i_ring - 1);

    auto z = 1 - i_ring * i_ring / fact2;
    auto phi = (i_phi - 0.5) * kPi / (2 * i_ring);

    FillVec(z, phi, data.data() + i);
  }

  // Equatorial region
  for (size_t i = n_cap; i < n_eq; i++) {
    size_t ip = i - n_cap;
    auto i_ring = ip / nl4 + n_side;
    auto i_phi = ip % nl4;
    auto f_odd = 0.5 * ((i_ring + n_side) % 2 - 1.0);

    auto z = (nl2 * 1.0 - i_ring) / fact1;
    auto phi = (i_phi - f_odd) * kPi / nl2;

    FillVec(z, phi, data.data() + i);
  }

  // South cap
  for (size_t i = n_eq; i < n_pix; i++) {
    auto ip = n_pix - i;
    auto hip = ip / 2.0;
    auto fihip = std::trunc(hip);
    auto i_ring = std::trunc(std::sqrt(hip - std::sqrt(fihip))) + 1;
    auto i_phi = 4 * i_ring + 1 - (ip - 2 * i_ring * (i_ring - 1));

    auto z = -1 + i_ring * i_ring / fact2;
    auto phi = (i_phi - 0.5) * kPi / (2 * i_ring);

    FillVec(z, phi, data.data() + i);
  }

  float dr = std::sqrt(4 * kPi / n_pix) * 90 / kPi;  // degree

  return SphGrid{ dr, data };
}

}  // namespace halo_pm
