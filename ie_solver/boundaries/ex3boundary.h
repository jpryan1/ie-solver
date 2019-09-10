// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_EX3BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_EX3BOUNDARY_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {

class Ex3Boundary : public CubicBoundary {
 public:
  void initialize(int N, BoundaryCondition bc);

  void get_spline_points(std::vector<double>* outer_x0_spline_points,
                         std::vector<double>* outer_x1_spline_points);
                              
  void get_fin_spline_points(std::vector<double>* x0_points,
                                    std::vector<double>* x1_points);
                              
  double fin_theta = -M_PI/4.;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_EX3BOUNDARY_H_
