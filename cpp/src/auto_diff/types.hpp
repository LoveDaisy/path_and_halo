#ifndef AUTO_DIFF_TYPES_HPP_
#define AUTO_DIFF_TYPES_HPP_

#include <cstddef>
#include <cstring>

#include "auto_diff/common.hpp"

namespace halo_pm {
namespace ad {

struct NoneType {};

// =============== Vec (Data type, not Expr type) ===============
template <class T, size_t Len>
struct Vec {
  T data_[Len];
};

template <class T>
struct Vec<T, 1> {
  T data_[1];

  operator float() const { return data_[0]; }
};

using Vec3f = Vec<float, 3>;

// ---------- Vec + Vec ----------
template <class T, size_t Len>
Vec<T, Len> operator+(Vec<T, Len> a, Vec<T, Len> b) {
  Vec<T, Len> res{};
  for (size_t i = 0; i < Len; i++) {
    res.data_[i] = a.data_[i] + b.data_[i];
  }
  return res;
}

// ---------- Vec + Scalar ----------
template <class T, size_t Len>
Vec<T, Len> operator+(Vec<T, Len> a, T b) {
  Vec<T, Len> res{};
  for (size_t i = 0; i < Len; i++) {
    res.data_[i] = a.data_[i] + b;
  }
  return res;
}

// ---------- Scalar + Vec ----------
template <class T, size_t Len>
Vec<T, Len> operator+(T b, Vec<T, Len> a) {
  Vec<T, Len> res{};
  for (size_t i = 0; i < Len; i++) {
    res.data_[i] = b + a.data_[i];
  }
  return res;
}


// =============== Mat (Data type, not Expr type) ===============
template <class T, size_t R, size_t C>
struct Mat {
  T data_[R * C];
};

template <class T>
struct Mat<T, 1, 1> {
  T data_[1];

  operator float() const { return data_[0]; }
};

using Mat3x3f = Mat<float, 3, 3>;

// ---------- Mat + Mat ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator+(Mat<T, R, C> a, Mat<T, R, C> b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] + b.data_[i];
  }
  return res;
}

// ---------- Mat + Scalar ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator+(Mat<T, R, C> a, T b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] + b;
  }
  return res;
}

// ---------- Scalar + Mat ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator+(T b, Mat<T, R, C> a) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = b + a.data_[i];
  }
  return res;
}

// ---------- Mat - Mat ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator-(Mat<T, R, C> a, Mat<T, R, C> b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] - b.data_[i];
  }
  return res;
}

// ---------- Mat - Scalar ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator-(Mat<T, R, C> a, T b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] - b;
  }
  return res;
}

// ---------- Scalar - Mat ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator-(T b, Mat<T, R, C> a) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = b - a.data_[i];
  }
  return res;
}

// ---------- Mat * Scalar ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator*(Mat<T, R, C> a, T b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] * b;
  }
  return res;
}

// ---------- Scalar * Mat ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator*(T b, Mat<T, R, C> a) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = b * a.data_[i];
  }
  return res;
}

// ---------- Mat / Scalar ----------
template <class T, size_t R, size_t C>
Mat<T, R, C> operator/(Mat<T, R, C> a, T b) {
  Mat<T, R, C> res{};
  for (size_t i = 0; i < R * C; i++) {
    res.data_[i] = a.data_[i] / b;
  }
  return res;
}

}  // namespace ad
}  // namespace halo_pm

#endif  // AUTO_DIFF_TYPES_HPP_
