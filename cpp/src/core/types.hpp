#ifndef CORE_TYPES_H_
#define CORE_TYPES_H_

#include <cstddef>
#include <initializer_list>

namespace halo_pm {

template <class T, size_t N>
struct Vec {
  T data_[N];

  Vec() : data_{} {};
};


template <class T>
struct Vec<T, 1> {
  union {
    T x_;
    T data_[1];
  };

  Vec() : x_(0){};
  Vec(T x) : x_(x){};
};


template <class T>
struct Vec<T, 2> {
  union {
    struct {
      T x_;
      T y_;
    };
    T data_[2];
  };

  Vec() : x_(0), y_(0){};
  Vec(T x, T y) : x_(x), y_(y){};
};


template <class T>
struct Vec<T, 3> {
  union {
    struct {
      T x_;
      T y_;
      T z_;
    };
    T data_[3];
  };

  Vec() : x_(0), y_(0), z_(0){};
  Vec(T x, T y, T z) : x_(x), y_(y), z_(z){};
};


template <class T>
struct Vec<T, 4> {
  union {
    struct {
      T x_;
      T y_;
      T z_;
      T w_;
    };
    T data_[4];
  };

  Vec() : x_(0), y_(0), z_(0), w_(0){};
  Vec(T x, T y, T z, T w) : x_(x), y_(y), z_(z), w_(w){};
};


using Vec2f = Vec<float, 2>;
using Vec3f = Vec<float, 3>;
using Quatf = Vec<float, 4>;


template <class T, size_t R, size_t C>
struct Mat {
  T data_[R * C];

  Mat() : data_{} {}

  template <class... V>
  Mat(V... v) : data_{ static_cast<T>(v)... } {}

  T& operator[](size_t i) { return data_[i]; }
};


using Mat3x3f = Mat<float, 3, 3>;


// =============== operators ===============
// Negative
template <class T, size_t N>
Vec<T, N> operator-(const Vec<T, N>& x) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = -x.data_[i];
  }
  return res;
}

// Add
template <class T, size_t N>
Vec<T, N> operator+(const Vec<T, N>& l, const Vec<T, N>& r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l.data_[i] + r.data_[i];
  }
  return res;
}

template <class T, size_t N>
Vec<T, N> operator+(const Vec<T, N>& l, T r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l.data_[i] + r;
  }
  return res;
}

template <class T, size_t N>
Vec<T, N> operator+(T l, const Vec<T, N>& r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l + r.data_[i];
  }
  return res;
}

// Minus
template <class T, size_t N>
Vec<T, N> operator-(const Vec<T, N>& l, const Vec<T, N>& r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l.data_[i] - r.data_[i];
  }
  return res;
}

template <class T, size_t N>
Vec<T, N> operator-(const Vec<T, N>& l, T r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l.data_[i] - r;
  }
  return res;
}

template <class T, size_t N>
Vec<T, N> operator-(T l, const Vec<T, N>& r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l - r.data_[i];
  }
  return res;
}

// Scalar multiplication
template <class T, size_t N>
Vec<T, N> operator*(const Vec<T, N>& l, T r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l.data_[i] * r;
  }
  return res;
}

template <class T, size_t N>
Vec<T, N> operator*(T l, const Vec<T, N>& r) {
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res.data_[i] = l * r.data_[i];
  }
  return res;
}


}  // namespace halo_pm

#endif  // CORE_TYPES_H_
