// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_
#define IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_

#include "ie-solver/boundaries/boundary.h"

namespace ie_solver {

class RoundedSquare : public Boundary {
 public:
  enum BoundaryCondition {
    SINGLE_ELECTRON
  };

  void initialize(int N, int bc_enum);
  bool is_in_domain(const Vec2& a);

  void draw_line(int bc_index, int num_points, double start_x, double start_y,
                 double end_x, double end_y, bool normal_is_left, int bc_enum);
  void draw_quarter_circle(int bc_index, int num_points, double start_x,
                           double start_y, double end_x, double end_y,
                           bool convex, int bc_enum);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_ROUNDED_SQUARE_H_
