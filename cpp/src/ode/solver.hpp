#ifndef ODE_SOLVER_H_
#define ODE_SOLVER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>

#include "core/geo.hpp"
#include "core/math.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

namespace halo_pm {

template <class T, int OutputDim, int InputDim>
using Func = std::function<Vec<T, OutputDim>(const Vec<T, InputDim>&)>;

template <class T, int OutputDim, int InputDim>
using FuncAndDiff = std::function<std::tuple<Vec<T, OutputDim>, Mat<T, OutputDim, InputDim>>(const Vec<T, InputDim>&)>;


struct SolverOption {
  static constexpr double kDefaultEps = 1e-6;
  static constexpr double kDefaultContourH = 0.05;
  static constexpr size_t kDefaultMaxEval = 15;
  static constexpr size_t kDefaultMaxPts = 100;

  double eps_;
  double h_;
  size_t max_eval_;
  size_t max_pts_;

  SolverOption() : eps_(kDefaultEps), h_(kDefaultContourH), max_eval_(kDefaultMaxEval), max_pts_(kDefaultMaxPts) {}
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
  while (dy.norm() > option.eps_ && status.func_eval_cnt_ <= option.max_eval_) {
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
    Vec<T, InputDim> dx = v.leftCols(rank) * s.head(rank).asDiagonal() * u.leftCols(rank).transpose() * dy;
    LOG_DEBUG("jac: \n%s", ObjLogFormatter<Mat<T, OutputDim, InputDim>>{ jac }.Format());
    LOG_DEBUG("dx: %s", ObjLogFormatter<Vec<T, InputDim>>{ dx }.Format());

    // Linear search
    constexpr float kBeta = 0.6f;
    float h = 1.0f;
    for (;; h *= kBeta) {
      x = x0 + dx * h;
      std::tie(y, std::ignore) = func_jac(x);
      status.func_eval_cnt_++;
      LOG_DEBUG("h: %.6f, x: %s, y: %s", h, ObjLogFormatter<Vec<T, InputDim>>{ x }.Format(),
                ObjLogFormatter<Vec<T, OutputDim>>{ y }.Format());

      if (!y.hasNaN() && ((yq - y).norm() <= option.eps_ || status.func_eval_cnt_ > option.max_eval_ || h < kHMin ||
                          (yq - y).norm() <= dy.norm() * h * 0.8)) {
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
  status.solved_ = dy.norm() < option.eps_;
  return std::make_tuple(x0, status);
}


// =============== FindContour ===============
struct ContourStatus {
  bool closed_;
  size_t func_eval_cnt_;

  ContourStatus() : closed_(false), func_eval_cnt_(0) {}
};


template <class T, int OutputDim, int InputDim>
std::tuple<Curve<T, InputDim>, ContourStatus>                         // Result & status
SearchDirection(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,  // Function
                const Vec<T, InputDim>& x_start, int direction,       // Start point & direction
                SolverOption option = SolverOption{}) {               // Option
}

template <class T, int OutputDim, int InputDim>
std::tuple<Curve<T, InputDim>, ContourStatus>                     // Result & status
FindContour(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,  // Function
            const Vec<T, InputDim>& x_start,                      // Start point
            SolverOption option = SolverOption{}) {               // Option

  // 1. Search forward direction
  auto [contour_f, status_f] = SearchDirection(func_jac, x_start, 1, option);

  if (status_f.closed_) {
    return std::make_tuple(contour_f, status_f);
  }

  // 2. Search backward direction
  auto [contour_b, status_b] = SearchDirection(func_jac, x_start, -1, option);
  ContourStatus status;
  status.func_eval_cnt_ = status_f.func_eval_cnt_ + status_b.func_eval_cnt_;

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
    std::make_tuple(contour, status);
  }
}

}  // namespace halo_pm

#endif  // ODE_SOLVER_H_
