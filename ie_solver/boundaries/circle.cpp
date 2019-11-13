// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie_solver/boundaries/circle.h"
#include "ie_solver/log.h"

namespace ie_solver {

void Circle::initialize(int N, BoundaryCondition bc) {
  boundary_shape = CIRCLE;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  for (int i = 0; i < N; i++) {
    double ang = i * 2.0 * M_PI / N;
    double x = 0.5 + 0.25 * cos(ang);
    double y = 0.5 + 0.25 * sin(ang);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(cos(ang));
    normals.push_back(sin(ang));
    curvatures.push_back(4);  // 1/r, r=0.25
    weights.push_back(M_PI / (N * 2));
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(N, 1);
    apply_boundary_condition(0, N, SINGLE_ELECTRON);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, N, bc);
  }
}

bool Circle::is_in_domain(const Vec2& a) {
  double x = a.a[0] - 0.5;
  double y = a.a[1] - 0.5;
  double eps = 1e-2;

  double dist = sqrt(pow(x, 2) + pow(y, 2));
  if (dist + eps > 0.25) return false;
  return true;
}

}  // namespace ie_solver
