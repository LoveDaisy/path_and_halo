#ifndef CORE_TYPES_H_
#define CORE_TYPES_H_

#include <Eigen/Eigen>
#include <cstddef>
#include <cstdio>

namespace halo_pm {

template <class T, int Dim>
using Vec = Eigen::Matrix<T, Dim, 1>;

template <class T, int R, int C>
using Mat = Eigen::Matrix<T, R, C>;

using Vec2f = Eigen::Vector2f;
using Vec3f = Eigen::Vector3f;
using Vec4f = Eigen::Vector4f;
using Vec2d = Eigen::Vector2d;
using Vec3d = Eigen::Vector3d;
using Vec4d = Eigen::Vector4d;

using Mat3x3f = Eigen::Matrix3f;
using Mat3x4f = Eigen::Matrix<float, 3, 4>;
using Mat4x4f = Eigen::Matrix4f;
using Mat3x3d = Eigen::Matrix3d;
using Mat3x4d = Eigen::Matrix<double, 3, 4>;
using Mat4x4d = Eigen::Matrix4d;

using Quatf = Eigen::Quaternionf;
using Quatd = Eigen::Quaterniond;


template <class T>
struct ObjLogFormatter {};


template <class T, int Dim>
struct ObjLogFormatter<Vec<T, Dim>> {
  static constexpr size_t kBufLen = 1024;
  char obj_buf_[kBufLen];
  const Vec<T, Dim>& vec_;

  ObjLogFormatter(const Vec<T, Dim>& vec) : vec_(vec){};

  operator const char*() {
    int offset = std::snprintf(obj_buf_, kBufLen, "[");
    for (auto x : vec_) {
      offset += std::snprintf(obj_buf_ + offset, kBufLen, "%.6f,", x);
    }
    std::snprintf(obj_buf_ + offset, kBufLen, "]");
    return obj_buf_;
  }

  const char* Format() { return static_cast<const char*>(*this); }
};


template <class T, int R, int C>
struct ObjLogFormatter<Mat<T, R, C>> {
  static constexpr size_t kBufLen = 1024;
  char obj_buf_[kBufLen];
  const Mat<T, R, C>& mat_;

  ObjLogFormatter(const Mat<T, R, C>& m) : mat_(m){};

  operator const char*() {
    int offset = 0;
    for (int r = 0; r < R - 1; r++) {
      for (auto x : mat_.row(r)) {
        offset += std::snprintf(obj_buf_ + offset, kBufLen, "%.6f,", x);
      }
      offset += std::snprintf(obj_buf_ + offset, kBufLen, "\n");
    }
    for (auto x : mat_.row(R - 1)) {
      offset += std::snprintf(obj_buf_ + offset, kBufLen, "%.6f,", x);
    }
    return obj_buf_;
  }

  const char* Format() { return static_cast<const char*>(*this); }
};


template <class T>
struct ObjLogFormatter<Eigen::Quaternion<T>> {
  static constexpr size_t kBufLen = 1024;
  char obj_buf_[kBufLen];
  const Eigen::Quaternion<T>& q_;

  ObjLogFormatter(const Eigen::Quaternion<T>& q) : q_(q){};

  operator const char*() {
    std::snprintf(obj_buf_, kBufLen, "[%.6f,%.6f,%.6f,%.6f]", q_.w(), q_.x(), q_.y(), q_.z());
    return obj_buf_;
  }

  const char* Format() { return static_cast<const char*>(*this); }
};

}  // namespace halo_pm

#endif  // CORE_TYPES_H_
