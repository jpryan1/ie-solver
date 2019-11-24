// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex2boundary.h"
#include "ie_solver/log.h"

#define OUTER_NUM_SPLINE_POINTS 28
// #define FIN_SPLINE_POINTS 4
// #define FIN_RAD 0.075
namespace ie_solver {


// USE ROUNDED SQUARE FOR EXTERIOR BOUNDARY

// move boundaries outwards to make laminar flow
void Ex2Boundary::get_spline_points(std::vector<double>* x0_points,
                                    std::vector<double>* x1_points) {
  for (int i = 0; i < 12; i++) {
    x0_points->push_back(-1 + 3 * (i / 12.0));
    x1_points->push_back(0.25);
  }
  for (int i = 0; i < 2; i++) {
    x0_points->push_back(2);
    x1_points->push_back(0.25 + 0.5 * i / 2.0);
  }

  for (int i = 0; i < 12; i++) {
    x0_points->push_back(2 - 3 * (i / 12.0));
    x1_points->push_back(0.75);
  }
  for (int i = 0; i < 2; i++) {
    x0_points->push_back(-1);
    x1_points->push_back(0.75 - 0.5 * i / 2.0);
  }
}

// void Ex2Boundary::get_fin_spline_points(std::vector<double>* x0_points,
//                                         std::vector<double>* x1_points) {
//   x0_points->push_back(0.025);
//   x1_points->push_back(0.);

//   x0_points->push_back(0.);
//   x1_points->push_back(0.075);

//   x0_points->push_back(-0.025);
//   x1_points->push_back(0.);

//   x0_points->push_back(0.);
//   x1_points->push_back(-0.075);

//   // Rotate by perturbation_parameters[0]

//   for (int i = x0_points->size() - 4; i < x0_points->size(); i++) {
//     double temp = cos(perturbation_parameters[0]) * (*x0_points)[i]
//                   - sin(perturbation_parameters[0]) * (*x1_points)[i];
//     (*x1_points)[i] = 0.5 + sin(perturbation_parameters[0]) * (*x0_points)[i]
//                       + cos(perturbation_parameters[0]) * (*x1_points)[i];
//     (*x0_points)[i] = 0.5 + temp;
//   }
// }


void Ex2Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::EX2;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  holes.clear();

  if (perturbation_parameters.size() == 0) {
    for (int i = 0; i <= 10; i++) perturbation_parameters.push_back(0.5);
  }

  int OUTER_NODES_PER_SPLINE = (3 * N / 4) / OUTER_NUM_SPLINE_POINTS;
  int NUM_CIRCLE_POINTS = (N / 4) / 8;
  // int FIN_NODES_PER_SPLINE = ((N / 4) / 4) / FIN_SPLINE_POINTS;

  Hole circle;
  // fin.center = Vec2(0.5, 0.5);
  // fin.radius = FIN_RAD;
  // fin.num_nodes = FIN_NODES_PER_SPLINE * FIN_SPLINE_POINTS;
  // holes.push_back(fin);

  circle.radius = 0.025;
  circle.num_nodes =  NUM_CIRCLE_POINTS;
  for (int i = 0; i <= 10; i++) {
    double x = i / 10.;
    circle.center = Vec2(x, perturbation_parameters[i]);
    holes.push_back(circle);
  }


  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  interpolate(false, OUTER_NODES_PER_SPLINE,
              outer_x0_cubics, outer_x1_cubics);
  num_outer_nodes = OUTER_NODES_PER_SPLINE * OUTER_NUM_SPLINE_POINTS;

  // std::vector<double> fin_x0_spline_points, fin_x1_spline_points;
  // get_fin_spline_points(&fin_x0_spline_points, &fin_x1_spline_points);

  // std::vector<std::vector<double>> fin_x0_cubics, fin_x1_cubics;
  // get_cubics(fin_x0_spline_points, fin_x1_spline_points,
  //            &fin_x0_cubics, &fin_x1_cubics);

  // interpolate(true, FIN_NODES_PER_SPLINE,
  //             fin_x0_cubics, fin_x1_cubics);
  for (int i = 0; i < holes.size(); i++) {
    Hole circle = holes[i];
    for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
      double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
      points.push_back(circle.center.a[0] + circle.radius * cos(ang));
      points.push_back(circle.center.a[1] + circle.radius * sin(ang));
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / circle.radius);
      weights.push_back(2.0 * M_PI * circle.radius / NUM_CIRCLE_POINTS);
    }
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(weights.size() * 2, 1);
    apply_boundary_condition(0, weights.size(), LEFT_TO_RIGHT_FLOW);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}


}  // namespace ie_solver
