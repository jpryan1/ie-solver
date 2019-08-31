// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {

class Ex1Boundary : public Boundary {
 public:

  enum Ex1BC {
    TANGENTIAL,
    REVERSE_TANGENTIAL,
    NORMAL,
    REVERSE_NORMAL
  };

  void initialize(int N, BoundaryCondition bc);
  bool is_in_domain(const Vec2& a);

  void get_spline_points(std::vector<double>* outer_x0_spline_points,
                         std::vector<double>* outer_x1_spline_points);
  void get_star_spline_points(double x, double y,
                              std::vector<double>* star_x0_spline_points,
                              std::vector<double>* star_x1_spline_points);
  void get_cubics(const std::vector<double>& x0_points,
                  const std::vector<double>& x1_points,
                  std::vector<std::vector<double>>* x0_cubics,
                  std::vector<std::vector<double>>* x1_cubics);

  void interpolate(int bc_index, bool is_interior, int nodes_per_spline, Ex1BC bc,
                   const std::vector<std::vector<double>>& x0_cubics_,
                   const std::vector<std::vector<double>>& x1_cubics_);


  void find_real_roots_of_cubic(const std::vector<double>& y_cubic,
                                std::vector<double>* t_vals);
  int num_right_intersections(double x, double y, int index);

  // v0 + v1 x + v2 x^2 + v3 x^3
  std::vector<std::vector<double>> all_cubics_x0, all_cubics_x1;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_EX1BOUNDARY_H_
