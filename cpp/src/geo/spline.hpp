#ifndef GEO_SPLINE_H_
#define GEO_SPLINE_H_

#include <cassert>
#include <cstddef>
#include <tuple>
#include <vector>

#include "core/math.hpp"
#include "core/types.hpp"

namespace halo_pm {
template <class X, class Y>
Y InterpLinear(const X& x0, const X& x1, const Y& y0, const Y& y1, const X& xq) {
  return (xq * (y0 - y1) + x0 * y1 - x1 * y0) / (x0 - x1);
}

template <class X, class Y>
Y InterpLinear(const std::vector<X>& x, const std::vector<Y>& y, X xq) {
  assert(x.size() == y.size());
  assert(x.size() > 1);

  // NOTE: Assume x is sorted in increasing order
  size_t i = 0;
  for (; i + 1 < x.size(); i++) {
    if (x[i] <= xq && xq < x[i + 1]) {
      break;
    }
  }

  X x0 = x[i];
  X x1 = x[i + 1];
  const auto& y0 = y[i];
  const auto& y1 = y[i + 1];
  return InterpLinear(x0, x1, y0, y1, xq);
}


template <class T, int Dim>
Curve<T, Dim> InterpLinear(const std::vector<T>& x, const Curve<T, Dim>& y, const std::vector<T>& xq) {
  Curve<T, Dim> res;
  res.reverse(xq.size());

  for (auto xi : x) {
    res.emplace_back(InterpLinear(x, y, xi));
  }

  return res;
}


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
std::tuple<Curve<T, Dim>, std::vector<T>, std::vector<T>>  // (interp_pts, interp_ss, refined_s0)
InterpCurve(const Curve<T, Dim>& pts, double ds) {
  assert(pts.size() > 1);

  int pt_num = static_cast<int>(pts.size());
  std::vector<T> s(pt_num, 0);
  std::vector<T> ss;

  for (int i = 1; i < pt_num; i++) {
    auto len = (pts[i] - pts[i - 1]).norm();
    s[i] = s[i - 1] + len;
  }

  Curve<T, Dim> res;
  double diff = ds;
  while (diff > ds * 0.1) {
    auto last_len = s.back();
    auto spline = InterpSpline(s, pts);

    ss.clear();
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
  return { res, ss, s };
}

}  // namespace halo_pm

#endif  // GEO_SPLINE_H_
