// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/boundaries/annulus.h"
#include "ie-solver/log.h"

namespace ie_solver {

void Annulus::initialize(int N, BoundaryCondition bc) {
  N = 6 * (N / 6);


  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * N, 1);
  } else {
    boundary_values = ie_Mat(N, 1);
  }
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();


  for (int i = 0; i < 5 * (N / 6); i++) {
    double ang = i * 2.0 * M_PI / (5 * (N / 6));
    double x = 0.5 + 0.25 * cos(ang);
    double y = 0.5 + 0.25 * sin(ang);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(cos(ang));
    normals.push_back(sin(ang));
    curvatures.push_back(4);
    weights.push_back(0.5 * M_PI / (5 * (N / 6)));

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


  for (int i = 5 * (N / 6); i < N; i++) {
    double ang = (i - (5 * (N / 6))) * 2.0 * M_PI / (N / 6);
    double x = 0.5 + 0.05 * cos(ang);
    double y = 0.5 + 0.05 * sin(ang);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(-cos(ang));
    normals.push_back(-sin(ang));
    curvatures.push_back(-20);  // 1/r, r=0.05
    weights.push_back((0.1 * M_PI) / (N / 6));

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
        boundary_values.set(2 * i , 0, -1.5*normals[2 * i  + 1]);
        boundary_values.set(2 * i  + 1, 0, 1.5*normals[2 * i ]);
        break;
    }
  }
}


bool Annulus::is_in_domain(const Vec2& a) {
  double x = a.a[0] - 0.5;
  double y = a.a[1] - 0.5;
  double eps = 1e-2;

  double dist = sqrt(pow(x, 2) + pow(y, 2));
  if (dist + eps > 0.25) return false;
  if (dist - eps < 0.05) return false;
  return true;
}

}  // namespace ie_solver
