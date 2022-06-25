#ifndef GEO_GEO_H_
#define GEO_GEO_H_

#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

#include "core/math.hpp"
#include "core/types.hpp"
#include "geo/spline.hpp"

namespace halo_pm {

/**
 * @brief Convert a longitude-latitude representation into a x-y-z representation.
 *
 * @param ll  [input]  (longitude, latitude) data, in degree
 * @param xyz [output] (x, y, z) data
 * @param num
 * @param ll_step_bytes  Step bytes to next data. If leave it as default 0, the function will regard input data as
 *                       contiguous and regard step = 2 * sizeof(float)
 * @param xyz_step_bytes Step bytes to next data. Similar to above.
 */
void Ll2Xyz(const Vec2f* ll, Vec3f* xyz,                           // input & output, ll in degree
            size_t num = 1,                                        // data number
            size_t ll_step_bytes = 0, size_t xyz_step_bytes = 0);  // step of input & output

// A convenience overload for ONE data
Vec3f Ll2Xyz(const Vec2f& ll);


void Xyz2Ll(const Vec3f* xyz, Vec2f* ll,                           // input & output, ll in degree, xyz normalized
            size_t num = 1,                                        // data number
            size_t xyz_step_bytes = 0, size_t ll_step_bytes = 0);  // step of input & output

Vec2f Xyz2Ll(const Vec3f& xyz);


/**
 * @brief Convert axis-angle representation of a rotation into matrix representation.
 *
 * @param llr [input]  (longitude, latitude, roll) axis-angle representation, in degree.
 * @param mat [output] matrix representation. Row major.
 * @param num
 * @param llr_step_bytes
 * @param mat_step_bytes
 */
void Llr2Mat(const Vec3f* llr, Mat3x3f* mat,                         // input & output, llr in degree
             size_t num = 1,                                         // data number
             size_t llr_step_bytes = 0, size_t mat_step_bytes = 0);  // step of input & output


Quatf Llr2Quat(const Vec3f& llr);


Vec3f Quat2Llr(const Quatf& q);

Vec3f Quat2Llr(const Vec4f& qv);


class AxisPdf {
 public:
  AxisPdf();
  AxisPdf(float zenith_mean, float zenith_std);
  AxisPdf(float zenith_mean, float zenith_std, float roll_mean, float roll_std);

  float operator()(const Vec3f& llr) const;

 private:
  float zenith_mean_;  // all in degree
  float zenith_std_;   // std < 0 means it is uniform distributed
  float roll_mean_;
  float roll_std_;

  float zenith_int_c_;
  float roll_int_c_;
};


template <class... VX>
auto NormalizeExpr(VX... vx) {
  auto f = sqrt((... + (vx * vx)));
  return std::make_tuple(vx / f...);
}


/**
 * @brief Rotate a point by a given quaternion.
 *
 * @param quat [input]  The given quaternion. **ONLY ONE** data, 4 floats [xi, yj, zk, w]
 * @param xyz0 [input]  Points to be rotated.
 * @param xyz1 [output] Points rotated. It can be the same as xyz0, which indicates rotate points inplace.
 * @param num Data number for xyz data.
 * @param xyz0_step_bytes
 * @param xyz1_step_bytes
 */
void RotateByQuat(const Quatf& quat, const Vec3f* xyz0, Vec3f* xyz1,        // input & output
                  size_t num = 1,                                           // data number
                  size_t xyz0_step_bytes = 0, size_t xyz1_step_bytes = 0);  // steps


template <class QW, class QX, class QY, class QZ,  // quaternion
          class VX, class VY, class VZ>            // input vector
auto                                               // output vector
QuatRotExpr(QW qw0, QX qx0, QY qy0, QZ qz0,        // quaternion
            VX vx, VY vy, VZ vz) {                 // input vector
  auto [qw, qx, qy, qz] = NormalizeExpr(qw0, qx0, qy0, qz0);

  auto tmp_qw = -qx * vx - qy * vy - qz * vz;
  auto tmp_qx = qw * vx + qy * vz - qz * vy;
  auto tmp_qy = qw * vy - qx * vz + qz * vx;
  auto tmp_qz = qw * vz + qx * vy - qy * vx;

  auto ox = -tmp_qw * qx + tmp_qx * qw - tmp_qy * qz + tmp_qz * qy;
  auto oy = -tmp_qw * qy + tmp_qx * qz + tmp_qy * qw - tmp_qz * qx;
  auto oz = -tmp_qw * qz - tmp_qx * qy + tmp_qy * qx + tmp_qz * qw;

  return std::make_tuple(ox, oy, oz);
}


template <class T, int Dim>
struct PointLineDistanceInfo {
  T distance_;
  Vec<T, Dim> nearest_point_;

  operator T() const { return distance_; }
};


template <class T, int Dim>
PointLineDistanceInfo<T, Dim>  // info: (distance, nearest_point)
Point2LineDistance(const Vec<T, Dim>& p, const Vec<T, Dim>& p1, const Vec<T, Dim>& p2) {
  auto d = p2 - p1;
  auto d2 = d.squaredNorm();

  PointLineDistanceInfo<T, Dim> info{};
  if (d2 < 1e-8) {
    info.distance_ = (p - p1).norm();
    info.nearest_point_ = p1;
  } else {
    T t = ((p - p1).dot(d)) / d2;
    t = std::clamp(t, static_cast<T>(0), static_cast<T>(1));
    info.nearest_point_ = p1 + d * t;
    info.distance_ = (p - info.nearest_point_).norm();
  }
  return info;
}


template <class T, int Dim>
PointLineDistanceInfo<T, Dim>  // info
DistanceToPolyLine(const Vec<T, Dim>& p, const Curve<T, Dim>& poly_line) {
  PointLineDistanceInfo<T, Dim> res{ std::numeric_limits<T>::max(), Vec<T, Dim>{} };
  for (size_t i = 0; i + 1 < poly_line.size(); i++) {
    auto info = Point2LineDistance(p, poly_line[i], poly_line[i + 1]);
    if (info.distance_ < res.distance_) {
      res = info;
    }
  }
  return res;
}


template <class T, int Dim>
std::tuple<bool, Curve<T, Dim>>  // Check result & reduced curve
CheckLoopAndReduce(const Curve<T, Dim>& pts, double eps, double ds) {
  assert(pts.size() > 1);

  Curve<T, Dim> res = pts;
  bool closed = false;
  for (double d = 0; d < eps && !res.empty(); /* nothing */) {
    Curve<T, Dim> interp_pts;
    for (auto it = res.rbegin() + 1; it != res.rend(); it++) {
      if (interp_pts.empty() || (interp_pts.back() - *it).norm() > kDefaultFloatEps) {
        interp_pts.emplace_back(*it);
      }
    }
    if (ds > 0) {
      std::tie(interp_pts, std::ignore, std::ignore) = InterpCurve(interp_pts, ds);
    }

    const auto& p = res.back();
    d = DistanceToPolyLine(p, interp_pts);
    if (d < eps) {
      closed = true;
      res.pop_back();
    }
  }

  return { closed, res };
}

template <class T, int Dim>
T BendingAngle(const Vec<T, Dim>& x0, const Vec<T, Dim>& x1, const Vec<T, Dim>& x2) {
  auto d1 = x1 - x0;
  auto d2 = x2 - x1;

  auto d1_norm = d1.norm();
  auto d2_norm = d2.norm();
  if (d1_norm < kDefaultFloatEps || d2_norm < kDefaultFloatEps) {
    return static_cast<T>(0);
  } else {
    return std::acos(std::clamp(d1.dot(d2) / d1_norm / d2_norm, static_cast<T>(0), static_cast<T>(1)));
  }
}


template <class T, int Dim>
std::vector<T> CumulatedArcLength(const Curve<T, Dim>& pts) {
  std::vector<T> s;
  s.emplace_back(0);

  for (size_t i = 1; i < pts.size(); i++) {
    T len = (pts[i] - pts[i - 1]).norm();
    s.emplace_back(len + s.back());
  }

  return s;
}


std::tuple<Curve3f, float>  // (projected intersect polygon, area factor)
ProjectPolygonIntersection(const Curve3f& p1, const Curve3f& p2, const Vec3f& dir);

}  // namespace halo_pm

#endif  // GEO_GEO_H_
