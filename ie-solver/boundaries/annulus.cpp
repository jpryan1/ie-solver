// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/boundaries/annulus.h"
#include "ie-solver/log.h"

namespace ie_solver {

void Annulus::initialize(int N, BoundaryCondition bc) {
  //            SHAPE PROTOCOL
  //
  // Except in trivial cases, we need to be careful if we want to evenly
  // distribute discretization points on the boundary, while compartmentalizing
  // drawing functions. For that purpose, we establish a SCALE_UNIT and require
  // that all draw functions take some multiple of this unit. For example, if
  // the ratio of shape A's points to shape B's points should be x/y, then shape
  // A can have x SCALE_UNIT's and shape B can have y.
  //
  //      Annulus
  // Since we allow the holder of this boundary to change it, we will not
  // enforce rigid units here, but instead opt to let N define the num of nodes
  // on the outer boundary, then add nodes to inner boundaries as necessary

  boundary_condition = bc;
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

  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * num_points, 1);
  } else {
    boundary_values = ie_Mat(num_points, 1);
  }

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
    double potential;
    switch (boundary_condition) {
      case BoundaryCondition::SINGLE_ELECTRON:
        potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
        boundary_values.set(i, 0, potential);
        break;
      case BoundaryCondition::ALL_ONES:
        boundary_values.set(i, 0, 1.0);
        break;
      case BoundaryCondition::BUMP_FUNCTION: {
        double N = boundary_values.height();
        // This x_val is -1 at i=0 and +1 at i=N-1
        double x_val = ((2.0 * i + 1.0 - N) / (N - 1.0));
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

      double potential;
      switch (boundary_condition) {
        case BoundaryCondition::SINGLE_ELECTRON:
          potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
          boundary_values.set(i, 0, potential);
          break;
        case BoundaryCondition::ALL_ONES:
          boundary_values.set(i, 0, 0.0);
          break;
        case BoundaryCondition::BUMP_FUNCTION: {
          double N = boundary_values.height();
          // This x_val is -1 at i=0 and +1 at i=N-1
          double x_val = ((2.0 * i + 1.0 - N) / (N - 1.0));
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
