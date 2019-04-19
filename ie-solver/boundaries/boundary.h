// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_BOUNDARY_H_

#include <vector>
#include "ie-solver/vec2.h"
#include "ie-solver/ie_mat.h"

namespace ie_solver {

struct Hole {
  Vec2 center;
  double radius;
};

class Boundary {
 public:
  int perturbation_size = 0;
  std::vector<double> points, normals, curvatures, weights;
  std::vector<Hole> holes;
  ie_Mat boundary_values;
  enum BoundaryShape {
    CIRCLE,
    ROUNDED_SQUARE,
    ROUNDED_SQUARE_WITH_BUMP,
    SQUIGGLY,
    ELLIPSES,
    CUBIC_SPLINE
  };
  enum BoundaryCondition {
    SINGLE_ELECTRON,
    ALL_ONES,
    BUMP_FUNCTION,
    STOKES
  };
  BoundaryShape boundary_shape;
  BoundaryCondition boundary_condition;
  virtual void initialize(int n, BoundaryCondition bc) = 0;
  virtual bool is_in_domain(const Vec2& a) = 0;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_BOUNDARY_H_
