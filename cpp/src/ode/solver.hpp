#ifndef ODE_SOLVER_H_
#define ODE_SOLVER_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>

#include "core/math.hpp"
#include "core/types.hpp"

namespace halo_pm {

template <class T, size_t OutputDim, size_t InputDim>
using Func = std::function<Vec<float, OutputDim>(const Vec<float, InputDim>&)>;

template <class T, size_t OutputDim, size_t InputDim>
using FuncAndDiff = std::function<std::tuple<Vec<T, OutputDim>, Mat<T, OutputDim, InputDim>>(const Vec<T, InputDim>&)>;


struct SolverOption {
  static constexpr double kDefaultEps = 1e-8;
  static constexpr size_t kDefaultMaxEval = 15;
  static constexpr size_t kDefaultMaxPts = 100;

  double eps_;
  size_t max_eval_;
  size_t max_pts_;

  SolverOption() : eps_(kDefaultEps), max_eval_(kDefaultMaxEval), max_pts_(kDefaultMaxPts) {}
};


struct SolutionStatus {
  bool solved_;
  size_t func_eval_cnt_;

  SolutionStatus() : solved_(false), func_eval_cnt_(0){};
};

template <class T, size_t OutputDim, size_t InputDim>
std::tuple<Vec<T, OutputDim>, SolutionStatus>                          // Result & status
FindSolution(const FuncAndDiff<T, OutputDim, InputDim>& func_jac,      // The function
             const Vec<T, InputDim>& x0, const Vec<T, OutputDim>& yq,  // Start point & target value
             SolverOption option = SolverOption{}) {                   // Options
  SolutionStatus status{};

  auto [y0, jac] = func_jac(x0);
  status.func_eval_cnt_++;
  if (std::any_of(y0.data_, y0.data_ + OutputDim, std::isnan)) {
    return std::make_tuple(Vec<T, OutputDim>{}, status);
  }

  constexpr double kHMin = 0.1;

  auto dy = yq - y0;
  auto x = x0;
  double h = 1.0;
  while (dy.norm() > option.eps_ && status.func_eval_cnt_ < option.max_eval_ && h > kHMin) {
    // Find deepest gradient
    using JacMat = Eigen::Matrix<float, OutputDim, InputDim>;
    Eigen::Map<JacMat> jac_map(jac.data_);
    Eigen::JacobiSVD<JacMat, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(jac_map);
  }
}

}  // namespace halo_pm

#endif  // ODE_SOLVER_H_
