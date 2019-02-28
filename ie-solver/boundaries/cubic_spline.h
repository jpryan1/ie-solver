// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_
#define IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_

#include <vector>
#include "ie-solver/boundaries/boundary.h"

namespace ie_solver {

class CubicSpline : public Boundary {
 public:
  void initialize(int N, BoundaryCondition bc);
  bool is_in_domain(const Vec2& a);

  void get_spline_points();
  void get_cubics();
  void interpolate();


  void find_real_roots_of_cubic(const std::vector<double>& y_cubic,
                                std::vector<double>* t_vals);
  int num_right_intersections(double x, double y, int index);

  int num_spline_points = 10;
  std::vector<double> x0_spline_points, x1_spline_points;

  // v0 + v1 x + v2 x^2 + v3 x^3
  std::vector<std::vector<double>> x0_cubics, x1_cubics;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_CUBIC_SPLINE_H_
