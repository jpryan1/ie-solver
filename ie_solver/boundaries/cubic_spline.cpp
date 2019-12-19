// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/log.h"

#define NUM_SPLINE_POINTS 20

namespace ie_solver {


void CubicSpline::get_spline_points(std::vector<double>* x0_spline_points,
                                    std::vector<double>* x1_spline_points) {
  for (int i = 0; i < NUM_SPLINE_POINTS; i++) {
    double ang = 2 * M_PI * (i / (NUM_SPLINE_POINTS + 0.));

    double x =  0.375 * cos(ang) * (sin(5 * ang) + 4);
    double y =  0.375 * sin(ang) * (sin(5 * ang) + 4);

    x0_spline_points->push_back(0.5 + x);
    x1_spline_points->push_back(0.5 + y);
  }
}


void CubicSpline::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  int NODES_PER_SPLINE = N / NUM_SPLINE_POINTS;
  // Currently, just get a working interpolation going

  std::vector<double> x0_spline_points, x1_spline_points;
  get_spline_points(&x0_spline_points, &x1_spline_points);
  std::vector<std::vector<double>> x0_cubics, x1_cubics;

  get_cubics(x0_spline_points, x1_spline_points, &x0_cubics, &x1_cubics);
  interpolate(false,  NODES_PER_SPLINE, x0_cubics,
              x1_cubics);
  num_outer_nodes = NODES_PER_SPLINE * NUM_SPLINE_POINTS;

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(weights.size(), 1);
    apply_boundary_condition(0, N, SINGLE_ELECTRON);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}

}  // namespace ie_solver
