// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie_solver/boundaries/annulus.h"
#include "ie_solver/log.h"

namespace ie_solver {

void Annulus::initialize(int N, BoundaryCondition bc) {

  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  if (holes.size() == 0) {
    Hole hole;
    hole.center = Vec2(0.6, 0.5);
    hole.radius = 0.025;
    holes.push_back(hole);
    hole.center = Vec2(0.4, 0.5);
    holes.push_back(hole);
  }

  int hole_nodes = N / 5;
  int num_points = N + hole_nodes * holes.size();
  for (int i = 0; i < N; i++) {
    double ang = i * 2.0 * M_PI / N;
    double x = 0.5 + 0.25 * cos(ang);
    double y = 0.5 + 0.25 * sin(ang);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(cos(ang));
    normals.push_back(sin(ang));
    curvatures.push_back(4);
    weights.push_back(0.5 * M_PI / N);
  }

  for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
    Hole hole = holes[hole_idx];
    int start_idx = N + hole_nodes * hole_idx;
    int end_idx = N + hole_nodes * (hole_idx + 1);
    for (int i = start_idx; i < end_idx; i++) {
      double ang = (i - start_idx) * 2.0 * M_PI / (end_idx - start_idx);
      double x = hole.center.a[0] + hole.radius * cos(ang);
      double y = hole.center.a[1] + hole.radius * sin(ang);
      points.push_back(x);
      points.push_back(y);
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / hole.radius);  // 1/r
      weights.push_back((2 * hole.radius * M_PI) / (end_idx - start_idx));
    }
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(num_points, 1);
    apply_boundary_condition(0, N, ALL_ONES);
    apply_boundary_condition(N, num_points, ALL_ZEROS);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, num_points, bc);
  }
}


bool Annulus::is_in_domain(const Vec2& a) {
  double x = a.a[0] - 0.5;
  double y = a.a[1] - 0.5;
  double eps = 1e-2;

  double dist = sqrt(pow(x, 2) + pow(y, 2));
  if (dist + eps > 0.25) return false;
  for (Hole hole : holes) {
    Vec2 r = a - hole.center;
    if (r.norm() - eps < hole.radius) return false;
  }
  return true;
}

}  // namespace ie_solver
