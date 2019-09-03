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


void CubicSpline::get_spline_points( std::vector<double>* x0_spline_points,
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
  get_cubics(x0_spline_points, x1_spline_points, &all_cubics_x0, &all_cubics_x1);
  interpolate(0, false,  NODES_PER_SPLINE, boundary_condition, all_cubics_x0, all_cubics_x1);
  //            SHAPE PROTOCOL
  //
  // Except in trivial cases, we need to be careful if we want to evenly
  // distribute discretization points on the boundary, while compartmentalizing
  // drawing functions. For that purpose, we establish a SCALE_UNIT and require
  // that all draw functions take some multiple of this unit. For example, if
  // the ratio of shape A's points to shape B's points should be x/y, then shape
  // A can have x SCALE_UNIT's and shape B can have y.
  //
  //    Rounded Square
  //
  //    Shape       Num of SCALE_UNIT's     Num of Shape in Boundary
  //    line        2                       24
  //    corner      3                       4
  // __________________________________________
  //
  //    Num of SCALE_UNIT's = 2*24 + 3*4 = 60
}



bool CubicSpline::is_in_domain(const Vec2& a) {
  const double* v = a.a;
  int intersections = 0;
  for (int i = 0; i < NUM_SPLINE_POINTS; i++) {
    int right_intersections = num_right_intersections(v[0], v[1], i);
    if (right_intersections == -1) return false;
    intersections += right_intersections;
  }
  if (intersections % 2 == 0) {
    return false;
  }

  return true;
}

}  // namespace ie_solver
