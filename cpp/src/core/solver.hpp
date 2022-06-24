#ifndef ODE_SOLVER_H_
#define ODE_SOLVER_H_

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <tuple>

#include "core/math.hpp"
#include "core/types.hpp"
#include "geo/geo.hpp"
#include "util/log.hpp"

namespace halo_pm {

template <class T, int OutputDim, int InputDim>
using Func = std::function<Vec<T, OutputDim>(const Vec<T, InputDim>&)>;

template <class T, int OutputDim, int InputDim>
using FuncAndDiff = std::function<std::tuple<Vec<T, OutputDim>, Mat<T, OutputDim, InputDim>>(const Vec<T, InputDim>&)>;


struct SolverOption {
  static constexpr double kDefaultAbsEps = 5e-6;
  static constexpr double kDefaultRelEps = 1e-7;
  static constexpr double kDefaultStepH = 0.05;
  static constexpr size_t kDefaultMaxEval = 15;
  static constexpr size_t kDefaultMaxPts = 100;

  double abs_eps_;
  double rel_eps_;
  double h_;
  size_t max_eval_;
  size_t max_pts_;

  SolverOption()
      : abs_eps_(kDefaultAbsEps), rel_eps_(kDefaultRelEps), h_(kDefaultStepH), max_eval_(kDefaultMaxEval),
        max_pts_(kDefaultMaxPts) {}
};


// =============== FindSolution ===============
struct SolutionStatus {
  bool solved_;
  size_t func_eval_cnt_;

  SolutionStatus() : solved_(false), func_eval_cnt_(0){};
};

template <class T, int OutputDim, int InputDim>
std::tuple<Vec<T, InputDim>, SolutionStatus>                                // Result & status
FindSolution(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,           // The function
             const Vec<T, InputDim>& x_start, const Vec<T, OutputDim>& yq,  // Start point & target value
             SolverOption option = SolverOption{}) {                        // Options
  LOG_DEBUG("yq: %s", ObjLogFormatter<Vec<T, OutputDim>>{ yq }.Format());

  SolutionStatus status{};

  Vec<T, InputDim> x0 = x_start;
  auto [y, jac] = func_jac(x0);
  status.func_eval_cnt_++;
  LOG_DEBUG("x0: %s, y0: %s", ObjLogFormatter<Vec<T, InputDim>>{ x_start }.Format(),
            ObjLogFormatter<Vec<T, OutputDim>>{ y }.Format());

  if (y.hasNaN()) {
    return std::make_tuple(Vec<T, InputDim>{}, status);
  }

  constexpr double kHMin = 0.1;

  Vec<T, OutputDim> dy = yq - y;
  auto x = x0;
  Vec<T, InputDim> dx = x;
  while ((dy.norm() > option.abs_eps_ && dx.norm() / x0.norm() > option.rel_eps_) &&
         status.func_eval_cnt_ <= option.max_eval_) {
    LOG_DEBUG("dy: %s", ObjLogFormatter<Vec<T, OutputDim>>{ dy }.Format());

    // Find deepest gradient
    Eigen::JacobiSVD svd(jac, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto rank = svd.rank();
    LOG_DEBUG("jac rank: %d", rank);
    const auto& u = svd.matrixU();
    const auto& v = svd.matrixV();
    auto s = svd.singularValues();
    for (int i = 0; i < rank; i++) {
      s(i) = 1.0f / s(i);
    }
    dx = v.leftCols(rank) * s.head(rank).asDiagonal() * u.leftCols(rank).transpose() * dy;
    LOG_DEBUG("jac: \n%s", ObjLogFormatter<Mat<T, OutputDim, InputDim>>{ jac }.Format());
    LOG_DEBUG("dx: %s, |dx|: %.4e, |dx|/|x|: %.4e", ObjLogFormatter<Vec<T, InputDim>>{ dx }.Format(), dx.norm(),
              dx.norm() / x.norm());

    // Linear search
    constexpr float kBeta = 0.6f;
    float h = 1.0f;
    for (;; h *= kBeta) {
      x = x0 + dx * h;
      std::tie(y, std::ignore) = func_jac(x);
      status.func_eval_cnt_++;
      LOG_DEBUG("h: %.6f, x: %s, y: %s", h, ObjLogFormatter<Vec<T, InputDim>>{ x }.Format(),
                ObjLogFormatter<Vec<T, OutputDim>>{ y }.Format());

      if (!y.hasNaN() &&
          ((yq - y).norm() <= option.abs_eps_ || dx.norm() / x0.norm() <= option.rel_eps_ ||
           status.func_eval_cnt_ > option.max_eval_ || h < kHMin || (yq - y).norm() <= dy.norm() * h * 0.8)) {
        x0 = x;
        dy = yq - y;
        std::tie(std::ignore, jac) = func_jac(x);
        status.func_eval_cnt_++;
        break;
      }
    }

    if (h < kHMin) {
      break;
    }
  }
  status.solved_ = dy.norm() <= option.abs_eps_;
  return std::make_tuple(x0, status);
}


// =============== FindContour ===============
struct ContourStatus {
  bool closed_;
  size_t func_eval_cnt_;

  ContourStatus() : closed_(false), func_eval_cnt_(0) {}
};


template <class T, int Dim>
std::tuple<Vec<T, Dim>, Vec<T, Dim>>  // (new_x, dx)
rk4(const Func<T, Dim, Dim>& dx_func, const Vec<T, Dim>& x0, int direction, double h) {
  auto k1 = dx_func(x0);
  auto k2 = dx_func(x0 + direction * h / 2.0 * k1);
  auto k3 = dx_func(x0 + direction * h / 2.0 * k2);
  auto k4 = dx_func(x0 + direction * h * k3);

  auto dx = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  auto new_x = x0 + dx * direction * h;
  return std::make_tuple(new_x, dx);
}


template <class T, int InputDim, int OutputDim>
std::tuple<Vec<T, OutputDim>, Vec<T, InputDim>>  // (f, dx)
WrapToFdx(const FuncAndDiff<T, OutputDim, InputDim>& fdy, const Vec<T, InputDim>& x, const Vec<T, InputDim>& ref_dx) {
  auto [y, j] = fdy(x);
  Vec<T, InputDim> dx;
  dx.fill(std::numeric_limits<T>::quiet_NaN());

  if (y.hasNaN()) {
    return std::make_tuple(y, dx);
  }

  Eigen::JacobiSVD svd(j, Eigen::ComputeFullV);
  const auto& v = svd.matrixV();
  dx = v.rightCols(1);

  if (ref_dx.squaredNorm() > 1e-6 && dx.dot(ref_dx) < 0) {
    dx = -dx;
  }
  return std::make_tuple(y, dx);
}


template <class T, int OutputDim, int InputDim>
std::tuple<Curve<T, InputDim>, ContourStatus>                         // Result & status
SearchDirection(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,  // Function
                const Vec<T, InputDim>& x_start, int direction,       // Start point & direction
                SolverOption option = SolverOption{}) {               // Option
  LOG_DEBUG("SearchDirection!");
  Curve<T, InputDim> res;
  ContourStatus status;
  res.emplace_back(x_start);

  auto h = option.h_;
  auto h_min = h / 16.0;
  auto h_max = h * 16.0;
  int h_extension_cnt = 0;

  Vec<T, OutputDim> y0;
  Vec<T, InputDim> ref_dx;
  std::tie(y0, ref_dx) = WrapToFdx(func_jac, x_start, Vec<T, InputDim>{});
  status.func_eval_cnt_++;

  while (res.size() < option.max_pts_) {
    LOG_DEBUG("res.size: %zu", res.size());
    const auto& x0 = res.back();
    Vec<T, InputDim> x2;
    LOG_DEBUG("x0: %s", ObjLogFormatter<Vec<T, InputDim>>{ x0 }.Format());

    bool x2_solved = false;
    bool retry = false;
    constexpr T kBendingTh = 20.0 * kDegree2Rad;

    // Start x2-loop to find new point
    while (!x2_solved && h > h_min) {
      // Call RK4 for first try.
      Func<T, InputDim, InputDim> dx_func = [=, &func_jac, &ref_dx](const Vec<T, InputDim>& x) -> Vec<T, InputDim> {
        auto [y, dx] = WrapToFdx(func_jac, x, ref_dx);
        return dx;
      };  // func: x -> dx
      auto [x1, ref_dx1] = rk4(dx_func, x0, direction, h);
      status.func_eval_cnt_ += 4;
      LOG_DEBUG("x1: %s, h: %.6f, dx: %s", ObjLogFormatter<Vec<T, InputDim>>{ x1 }.Format(), h,
                ObjLogFormatter<Vec<T, InputDim>>{ ref_dx }.Format());

      // If NaN, shrink h
      if (x1.hasNaN()) {
        LOG_DEBUG("x1 has nan. shrink h");
        h *= 0.5;
        h_extension_cnt = 0;
        continue;
      }

      // Calculate error for first try.
      Vec<T, OutputDim> y1;
      std::tie(y1, std::ignore) = func_jac(x1);
      status.func_eval_cnt_++;
      auto x1_error = (y1 - y0).norm();
      LOG_DEBUG("y1: %s, x1_error: %.6f", ObjLogFormatter<Vec<T, OutputDim>>{ y1 }.Format(), x1_error);

      // Check bending angle
      T bending = 0;
      if (res.size() > 1) {
        auto it = res.crbegin();
        bending = BendingAngle(x1, *it, *(it + 1));
      }
      if (bending > kBendingTh) {
        h *= 0.5;
        h_extension_cnt = 0;
        retry = true;
        res.pop_back();
        LOG_DEBUG("bending %.6f too large. shrink h to %.6f and pop one point.", bending, h);
        break;
      }
      if (x1_error > std::max(option.abs_eps_ * 5e3, h * 0.1)) {
        h *= 0.5;
        h_extension_cnt = 0;
        LOG_DEBUG("x1_error too large. shrink h to %.6f", h);
        continue;
      }

      // Check error. Refine if neccessary
      x2 = x1;
      x2_solved = true;
      if (x1_error > option.abs_eps_) {
        LOG_DEBUG("x1_error: %.6f too large. refine x2.", x1_error, x2_solved);
        SolutionStatus x2_status;
        std::tie(x2, x2_status) = FindSolution(func_jac, x1, y0, option);
        status.func_eval_cnt_ += x2_status.func_eval_cnt_;
        x2_solved = x2_status.solved_;
        LOG_DEBUG("x2: %s solved: %d", ObjLogFormatter<Vec<T, InputDim>>{ x2 }.Format(), x2_solved);
      }
      ref_dx = ref_dx1;

      // Adaptive scheme
      T dx1 = (x1 - x2).norm();
      if (!x2_solved || dx1 > h * 0.05) {
        h *= 0.5;
        h_extension_cnt = 0;
        LOG_DEBUG("can't solve x2 or x2 changes too much. shrink h to %.6f", h);
      } else if (h < h_max * 0.5 && (dx1 < h * 0.01 || x1_error < std::max(option.abs_eps_ * 100, h * 0.002))) {
        h_extension_cnt++;
        if (h_extension_cnt > 2) {
          h *= 2;
          LOG_DEBUG("extent h to %.6f", h);
        }
      }
    }  // x2-loop end
    if (retry) {
      continue;
    }
    if (!x2_solved) {
      break;
    }

    res.emplace_back(x2);
    if (res.size() > 2) {
      std::tie(status.closed_, std::ignore) = CheckLoopAndReduce(res, option.h_ * 0.05, option.h_ * 0.5);
    }
    if (status.closed_) {
      res.pop_back();
      break;
    }
  }

  if (status.closed_) {
    res.emplace_back(res.front());
  }
  if (res.size() < 2) {
    res.clear();
  }
  LOG_DEBUG("search direction: %d complete! closed: %d", direction, status.closed_);
  return std::make_tuple(res, status);
}

template <class T, int OutputDim, int InputDim>
std::tuple<Curve<T, InputDim>, ContourStatus>                     // Result & status
FindContour(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,  // Function
            const Vec<T, InputDim>& x_start,                      // Start point
            SolverOption option = SolverOption{}) {               // Option

  // 1. Search forward direction
  auto [contour_f, status_f] = SearchDirection(func_jac, x_start, 1, option);
  LOG_DEBUG("status_f.closed: %d", status_f.closed_);

  if (status_f.closed_) {
    return std::make_tuple(contour_f, status_f);
  }

  // 2. Search backward direction
  auto [contour_b, status_b] = SearchDirection(func_jac, x_start, -1, option);
  ContourStatus status;
  status.func_eval_cnt_ = status_f.func_eval_cnt_ + status_b.func_eval_cnt_;
  LOG_DEBUG("status_b.closed: %d", status_b.closed_);

  Curve<T, InputDim> contour;
  if (contour_b.size() < 2) {
    contour = contour_f;
  } else if (status_b.closed_) {
    contour = contour_b;
  } else {
    contour = contour_b;
    std::reverse(contour.begin(), contour.end());
    contour.pop_back();  // remove start (contour_f has the same start)
    contour.insert(contour.end(), contour_f.begin(), contour_f.end());
  }
  // If empty
  if (contour.size() < 2) {
    status.closed_ = false;
    return std::make_tuple(Curve<T, InputDim>{}, status);
  }

  // Check close loop
  else {
    std::tie(status.closed_, contour) = CheckLoopAndReduce(contour, option.h_ * 0.1, -1);
    return std::make_tuple(contour, status);
  }
}

}  // namespace halo_pm

#endif  // ODE_SOLVER_H_
