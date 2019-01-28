// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_
#define IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_

#include "ie-solver/boundaries/boundary.h"

namespace ie_solver {

class RoundedSquare : public Boundary {
 public:
  void initialize(int N, BoundaryCondition bc);
  bool is_in_domain(const Vec2& a);

  void draw_line(int bc_index, int num_points, double start_x, double start_y,
                 double end_x, double end_y, bool normal_is_left);
  void draw_quarter_circle(int bc_index, int num_points, double start_x,
                           double start_y, double end_x, double end_y,
                           bool convex);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_
