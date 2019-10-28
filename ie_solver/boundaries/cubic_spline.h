// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_
#define IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {

class CubicSpline : public CubicBoundary {
 public:
  void initialize(int N, BoundaryCondition bc);

  void get_spline_points(std::vector<double>* x0_spline_points,
                         std::vector<double>* x1_spline_points);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_
