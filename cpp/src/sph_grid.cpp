#include "sph_grid.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <tuple>

#include "math.hpp"
#include "util/log.hpp"

namespace halo_pm {

void FillVec(double z, double phi, float* xyz) {
  auto sth = std::sqrt(1 - z) * std::sqrt(1 + z);
  auto x = sth * std::cos(phi);
  auto y = sth * std::sin(phi);

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
}

std::tuple<std::unique_ptr<float[]>, size_t> GenerateGrid(int level) {
  if (level <= 0) {
    return std::make_tuple(nullptr, 0);
  }

  size_t n_side = 1 << level;
  size_t n_pix = 12 * n_side * n_side;
  std::unique_ptr<float[]> data{ new float[n_pix * 3]{} };


  size_t n_cap = 2 * n_side * (n_side - 1);
  size_t n_eq = 2 * n_side * (5 * n_side + 1);

  // North cap
  for (size_t i = 0; i < n_cap; i++) {
    auto fact2 = 3.0 * n_side * n_side;

    auto hip = (i + 1) / 2.0;
    auto fihip = std::trunc(hip);
    auto i_ring = std::trunc(std::sqrt(hip - std::sqrt(fihip))) + 1;
    auto i_phi = (i + 1) - 2 * i_ring * (i_ring - 1);

    auto z = 1 - i_ring * i_ring / fact2;
    auto phi = (i_phi - 0.5) * kPi / (2 * i_ring);

    FillVec(z, phi, data.get() + i * 3);
  }

  // Equatorial region
  for (size_t i = n_cap; i < n_eq; i++) {
    size_t nl2 = 2 * n_side;
    size_t nl4 = 4 * n_side;
    auto fact1 = 1.5 * n_side;

    size_t ip = i - n_cap;
    auto i_ring = static_cast<size_t>(std::trunc(ip * 1.0 / nl4)) + n_side;
    auto i_phi = ip % nl4 + 1;
    auto f_odd = 0.5 * (1 + (i_ring + n_side) % 2);

    auto z = (nl2 * 1.0 - i_ring) / fact1;
    auto phi = (i_phi - f_odd) * kPi / (2.0 * n_side);

    FillVec(z, phi, data.get() + i * 3);
  }

  // South cap
  for (size_t i = n_eq; i < n_pix; i++) {
    auto fact2 = 3.0 * n_side * n_side;

    auto ip = n_pix - i;
    auto hip = ip / 2.0;
    auto fihip = std::trunc(hip);
    auto i_ring = std::trunc(std::sqrt(hip - std::sqrt(fihip))) + 1;
    auto i_phi = 4 * i_ring + 1 - (ip - 2 * i_ring * (i_ring - 1));

    auto z = -1 + i_ring * i_ring / fact2;
    auto phi = (i_phi - 0.5) * kPi / (2 * i_ring);

    FillVec(z, phi, data.get() + i * 3);
  }

  return std::make_tuple(std::move(data), n_pix);
}

}  // namespace halo_pm
