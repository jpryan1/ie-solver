// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {

class Ex1Boundary : public CubicBoundary {
 public:
  void initialize(int N, BoundaryCondition bc);

  void get_spline_points(std::vector<double>* outer_x0_spline_points,
                         std::vector<double>* outer_x1_spline_points);
  void get_star_spline_points(double x, double y,
                              std::vector<double>* star_x0_spline_points,
                              std::vector<double>* star_x1_spline_points);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_
