// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_ANNULUS_H_
#define IE_SOLVER_BOUNDARIES_ANNULUS_H_

#include "ie-solver/boundaries/boundary.h"

namespace ie_solver {

class Annulus : public Boundary {
 public:
  void initialize(int N, BoundaryCondition bc);
  bool is_in_domain(const Vec2& a);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_ANNULUS_H_
