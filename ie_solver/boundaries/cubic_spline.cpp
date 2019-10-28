// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/log.h"

#define NODES_PER_SPLINE 500
#define NUM_SPLINE_POINTS 10

namespace ie_solver {


void CubicSpline::get_spline_points(std::vector<double>* x0_spline_points,
                                    std::vector<double>* x1_spline_points) {
  for (int i = 0; i < 10; i++) {
    double ang = 2 * M_PI * (i / 10.0);

    double x = 0.5 + 0.2 * cos(ang);
    double y = 0.5 + 0.2 * sin(ang);
    if (i == 0) x += 0.2;
    if (i == 3) y += 0.1;
    if (i == 5) x -= 0.2;
    if (i == 8) y -= 0.05;

    x0_spline_points->push_back(x);
    x1_spline_points->push_back(y);
  }
}


void CubicSpline::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  // Currently, just get a working interpolation going
  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * NODES_PER_SPLINE * NUM_SPLINE_POINTS, 1);
  } else {
    boundary_values = ie_Mat(NODES_PER_SPLINE * NUM_SPLINE_POINTS, 1);
  }
  std::vector<double> x0_spline_points, x1_spline_points;
  get_spline_points(&x0_spline_points, &x1_spline_points);
  std::vector<std::vector<double>> x0_cubics, x1_cubics;

  get_cubics(x0_spline_points, x1_spline_points, &x0_cubics, &x1_cubics);
  interpolate(0, false,  NODES_PER_SPLINE, boundary_condition, x0_cubics,
              x1_cubics);
  num_outer_nodes = NODES_PER_SPLINE * NUM_SPLINE_POINTS;
}

}  // namespace ie_solver
