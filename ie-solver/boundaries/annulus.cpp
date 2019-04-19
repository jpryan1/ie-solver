// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/boundaries/annulus.h"
#include "ie-solver/log.h"

namespace ie_solver {

void Annulus::initialize(int N, BoundaryCondition bc) {
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  if (holes.size() == 0) {
    Hole hole;
    hole.center = Vec2(0.35, 0.35);
    hole.radius = 0.05;
    holes.push_back(hole);
    hole.center = Vec2(0.65, 0.65);
    hole.radius = 0.05;
    holes.push_back(hole);
  }
  int hole_dofs = 300;
  int num_points = N + hole_dofs * holes.size();

  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * num_points, 1);
  } else {
    boundary_values = ie_Mat(num_points, 1);
  }
  boundary_condition = bc;
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
    curvatures.push_back(4);
    weights.push_back(0.5 * M_PI / N);

    double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
    switch (boundary_condition) {
      case BoundaryCondition::SINGLE_ELECTRON:
        boundary_values.set(i, 0, potential);
        break;
      case BoundaryCondition::ALL_ONES:
        boundary_values.set(i, 0, 1.0);
        break;
      case BoundaryCondition::BUMP_FUNCTION: {
        double N = boundary_values.height();
        double x_val = -1 * ((N - 1.0 - i) / (N - 1.0))
                       + (i / (N - 1.0));
        potential = exp(-1.0 / (1.0 - pow(x_val, 2)));
        boundary_values.set(i, 0, potential);
        break;
      }
      case BoundaryCondition::STOKES:
        boundary_values.set(2 * i, 0, -normals[2 * i + 1]);
        boundary_values.set(2 * i + 1, 0, normals[2 * i]);
        break;
    }
  }
  for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
    Hole hole = holes[hole_idx];
    int start_idx = N + hole_dofs * hole_idx;
    int end_idx = N + hole_dofs * (hole_idx + 1);
    for (int i = start_idx; i < end_idx; i++) {
      double ang = (i - start_idx) * 2.0 * M_PI / (end_idx - start_idx);
      double x = hole.center.a[0] + hole.radius * cos(ang);
      double y = hole.center.a[1] + hole.radius * sin(ang);
      points.push_back(x);
      points.push_back(y);
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / hole.radius); // 1/r, r=0.05
      weights.push_back((2 * hole.radius * M_PI) / (end_idx - start_idx));

      double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
      switch (boundary_condition) {
        case BoundaryCondition::SINGLE_ELECTRON:
          boundary_values.set(i, 0, potential);
          break;
        case BoundaryCondition::ALL_ONES:
          boundary_values.set(i, 0, 1.0);
          break;
        case BoundaryCondition::BUMP_FUNCTION: {
          double N = boundary_values.height();
          double x_val = -1 * ((N - 1.0 - i) / (N - 1.0))
                         + (i  / (N - 1.0));
          potential = exp(-1.0 / (1.0 - pow(x_val, 2)));
          boundary_values.set(i , 0, potential);
          break;
        }
        case BoundaryCondition::STOKES:
          boundary_values.set(2 * i , 0, -1.5 * normals[2 * i  + 1]);
          boundary_values.set(2 * i  + 1, 0, 1.5 * normals[2 * i ]);

          break;
      }
    }
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
