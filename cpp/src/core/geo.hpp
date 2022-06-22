#ifndef CORE_GEO_H_
#define CORE_GEO_H_

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

#include "core/math.hpp"
#include "core/types.hpp"

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


// =============== Spline ===============
template <class T>
struct Spline {
  std::vector<T> x0_;
  Eigen::MatrixX<T> coef_;
};


template <class T, int Dim>
Curve<T, Dim> SampleSplinePoints(const Spline<T>& spline, const std::vector<T>& xq) {
  int output_num = static_cast<int>(xq.size());
  int spline_num = spline.coef_.rows() / 4;
  Curve<T, Dim> res(output_num);

  int idx = 0;
  for (int i = 0; i < output_num; i++) {
    while (xq[i] > spline.x0_[idx] && idx < spline_num) {
      idx++;
    }
    idx = std::max(idx - 1, 0);
    auto dx = xq[i] - spline.x0_[idx];
    Vec<T, 4> xx{ 1.0, dx, dx * dx, dx * dx * dx };
    const Eigen::MatrixX<T>& cc = spline.coef_.middleRows(idx * 4, 4);
    res[i] = cc.transpose() * xx;
  }
  return res;
}


template <class T, int YDim>
Spline<T> InterpSpline(const std::vector<T>& x, const Curve<T, YDim>& y) {
  assert(x.size() == y.size());
  assert(x.size() > 1);

  int input_num = x.size();
  int spline_num = input_num - 1;

  bool periodic = (y.front() - y.back()).norm() < kDefaultFloatEps;

  Eigen::MatrixX<T> mat = Eigen::MatrixX<T>::Zero(4 * spline_num, 4 * spline_num);
  Eigen::MatrixX<T> b = Eigen::MatrixX<T>::Zero(4 * spline_num, YDim);
  size_t offset = 0;
  for (int i = 0; i < spline_num; i++) {
    mat(i + offset, i * 4 + 0) = 1;
    b.row(i + offset) = y[i].transpose();
  }
  offset += spline_num;

  for (int i = 0; i < spline_num; i++) {
    auto dx = x[i + 1] - x[i];
    mat(i + offset, i * 4 + 0) = 1;
    mat(i + offset, i * 4 + 1) = dx;
    mat(i + offset, i * 4 + 2) = dx * dx;
    mat(i + offset, i * 4 + 3) = dx * dx * dx;
    b.row(i + offset) = y[i + 1].transpose();
  }
  offset += spline_num;

  for (int i = 0; i + 1 < spline_num; i++) {
    auto dx = x[i + 1] - x[i];
    mat(i + offset, i * 4 + 0) = 0;
    mat(i + offset, i * 4 + 1) = 1;
    mat(i + offset, i * 4 + 2) = 2 * dx;
    mat(i + offset, i * 4 + 3) = 3 * dx * dx;
    mat(i + offset, i * 4 + 4) = 0;
    mat(i + offset, i * 4 + 5) = -1;
  }
  offset += (spline_num - 1);

  for (int i = 0; i + 1 < spline_num; i++) {
    auto dx = x[i + 1] - x[i];
    mat(i + offset, i * 4 + 0) = 0;
    mat(i + offset, i * 4 + 1) = 0;
    mat(i + offset, i * 4 + 2) = 2;
    mat(i + offset, i * 4 + 3) = 6 * dx;
    mat(i + offset, i * 4 + 4) = 0;
    mat(i + offset, i * 4 + 5) = 0;
    mat(i + offset, i * 4 + 6) = -2;
  }
  offset += (spline_num - 1);

  if (!periodic) {
    mat(offset + 0, 0) = 0;
    mat(offset + 0, 1) = 0;
    mat(offset + 0, 2) = 2;

    mat(offset + 1, (spline_num - 1) * 4 + 0) = 0;
    mat(offset + 1, (spline_num - 1) * 4 + 1) = 0;
    mat(offset + 1, (spline_num - 1) * 4 + 2) = 2;
    mat(offset + 1, (spline_num - 1) * 4 + 3) = 6 * (x[spline_num] - x[spline_num - 1]);
  } else {
    auto dx = x[spline_num] - x[spline_num - 1];
    mat(offset + 0, 0) = 0;
    mat(offset + 0, 1) = 1;
    mat(offset + 0, (spline_num - 1) * 4 + 0) = 0;
    mat(offset + 0, (spline_num - 1) * 4 + 1) = -1;
    mat(offset + 0, (spline_num - 1) * 4 + 2) = -2 * dx;
    mat(offset + 0, (spline_num - 1) * 4 + 3) = -3 * dx * dx;

    mat(offset + 1, 0) = 0;
    mat(offset + 1, 1) = 0;
    mat(offset + 1, 2) = 2;
    mat(offset + 1, (spline_num - 1) * 4 + 0) = 0;
    mat(offset + 1, (spline_num - 1) * 4 + 1) = 0;
    mat(offset + 1, (spline_num - 1) * 4 + 2) = -2;
    mat(offset + 1, (spline_num - 1) * 4 + 3) = -6 * dx;
  }

  Eigen::MatrixX<T> coef = mat.partialPivLu().solve(b);
  return Spline<T>{ x, coef };
}


template <class T, int Dim>
Curve<T, Dim> InterpSpline(const std::vector<T>& x, const Curve<T, Dim>& y, const std::vector<T>& xq) {
  assert(x.size() == y.size());
  assert(x.size() > 1);

  auto spline = InterpSpline(x, y);
  return SampleSplinePoints<T, Dim>(spline, xq);
}


template <class T, int Dim>
Curve<T, Dim> InterpCurve(const Curve<T, Dim>& pts, double ds) {
  assert(pts.size() > 1);

  int pt_num = static_cast<int>(pts.size());
  std::vector<T> s(pt_num, 0);

  for (int i = 1; i < pt_num; i++) {
    auto len = (pts[i] - pts[i - 1]).norm();
    s[i] = s[i - 1] + len;
  }

  Curve<T, Dim> res;
  double diff = ds;
  while (diff > ds * 0.1) {
    auto last_len = s.back();
    auto spline = InterpSpline(s, pts);

    std::vector<T> ss;
    std::vector<int> s_idx;
    {
      double x = ds;
      for (auto si : s) {
        while (x < si && x < last_len) {
          if (x < si - ds * 0.5 && x > ss.back() + ds * 0.5) {
            ss.emplace_back(x);
          }
          x += ds;
        }
        s_idx.emplace_back(ss.size());
        ss.emplace_back(si);
      }
    }
    assert(s_idx.size() == pts.size());

    auto tmp_res = SampleSplinePoints<T, Dim>(spline, ss);
    res.clear();
    for (const auto& r : tmp_res) {
      if (res.empty() || (r - res.back()).norm() > kDefaultFloatEps) {
        res.emplace_back(r);
      }
    }

    s.clear();
    s.resize(ss.size());
    s[0] = 0;
    for (int i = 1; i < static_cast<int>(res.size()); i++) {
      auto len = (res[i] - res[i - 1]).norm();
      s[i] = s[i - 1] + len;
    }
    diff = std::abs(last_len - s.back());

    std::vector<T> tmp_s;
    for (auto i : s_idx) {
      tmp_s.emplace_back(s[i]);
    }
    s.swap(tmp_s);
  }
  return res;
}


template <class T, int Dim>
T Point2LineDistance(const Vec<T, Dim>& p, const Vec<T, Dim>& p1, const Vec<T, Dim>& p2) {
  auto d = p2 - p1;
  auto d2 = d.squaredNorm();

  if (d2 < 1e-8) {
    return (p - p1).norm();
  } else {
    T t = ((p - p1).dot(d)) / d2;
    t = std::clamp(t, static_cast<T>(0), static_cast<T>(1));
    return ((p1 + d * t) - p).norm();
  }
}


template <class T, int Dim>
T DistanceToPolyLine(const Vec<T, Dim>& p, const Curve<T, Dim>& poly_line) {
  T d = std::numeric_limits<T>::max();
  for (size_t i = 0; i + 1 < poly_line.size(); i++) {
    T curr_d = Point2LineDistance(p, poly_line[i], poly_line[i + 1]);
    d = std::min(curr_d, d);
  }
  return d;
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
      interp_pts = InterpCurve(interp_pts, ds);
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

}  // namespace halo_pm

#endif  // CORE_GEO_H_
