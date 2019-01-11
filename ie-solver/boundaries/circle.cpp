// Copyright 2019 John Paul Ryan
#include <cmath>
#include "ie-solver/boundaries/circle.h"
#include "ie-solver/log.h"

namespace ie_solver {

void Circle::initialize(int N, int bc_enum) {
  boundary_condition = ie_Mat(N, 1);
  if (bc_enum != BoundaryCondition::SINGLE_ELECTRON) {
    LOG::ERROR("Circle boundary can only do single electron bc currently.");
  }

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

    double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
    boundary_condition.set(i, 0, potential);
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
