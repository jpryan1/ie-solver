// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_EX2BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_EX2BOUNDARY_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {

class Ex2Boundary : public CubicBoundary {
 public:
  void initialize(int N, BoundaryCondition bc);

  void get_spline_points(std::vector<double>* x0_spline_points,
                         std::vector<double>* x1_spline_points);

  void get_star_spline_points(double x, double y,
                              std::vector<double>* x0_points,
                              std::vector<double>* x1_points);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_EX2BOUNDARY_H_
