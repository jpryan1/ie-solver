// Copyright 2019 John Paul Ryan
#ifndef _CIRCLE_H_
#define _CIRCLE_H_

#include "ie-solver/boundaries/boundary.h"

namespace ie_solver {

class Circle : public Boundary {
 public:
  enum BoundaryCondition {
    SINGLE_ELECTRON
  };

  void initialize(int N, int bc_enum);
  bool is_in_domain(const Vec2& a);
};

}  // namespace ie_solver

#endif
