// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex4boundary.h"
#include "ie_solver/log.h"

#define FIN_RAD 0.1
namespace ie_solver {


// USE ROUNDED SQUARE FOR EXTERIOR BOUNDARY

// move boundaries outwards to make laminar flow
void Ex4Boundary::get_spline_points(std::vector<double>* x0_points,
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

void Ex4Boundary::get_fin_spline_points(std::vector<double>* x0_points,
                                        std::vector<double>* x1_points) {
  x0_points->push_back(FIN_RAD / 3);
  x1_points->push_back(0.);

  x0_points->push_back(0.);
  x1_points->push_back(FIN_RAD);

  x0_points->push_back(-FIN_RAD / 3);
  x1_points->push_back(0.);

  x0_points->push_back(0.);
  x1_points->push_back(-FIN_RAD);

  // Rotate by perturbation_parameters[0]

  for (int i = x0_points->size() - 4; i < x0_points->size(); i++) {
    double temp = cos(perturbation_parameters[0]) * (*x0_points)[i]
                  - sin(perturbation_parameters[0]) * (*x1_points)[i];
    (*x1_points)[i] = 0.5 + sin(perturbation_parameters[0]) * (*x0_points)[i]
                      + cos(perturbation_parameters[0]) * (*x1_points)[i];
    (*x0_points)[i] = 0.5 + temp;
  }
}


void Ex4Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::EX4;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  holes.clear();

  if (perturbation_parameters.size() == 0) {
    perturbation_parameters.push_back(0);
  }


  int OUTER_NUM_SPLINE_POINTS = 28;
  int FIN_SPLINE_POINTS = 4;

  int OUTER_NODES_PER_SPLINE = (3 * N / 4) / OUTER_NUM_SPLINE_POINTS;
  int FIN_NODES_PER_SPLINE = (N / 4) / FIN_SPLINE_POINTS;

  num_outer_nodes = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;

  Hole fin;
  fin.center = Vec2(0.5, 0.5);
  fin.radius = FIN_RAD;
  fin.num_nodes = FIN_SPLINE_POINTS * FIN_NODES_PER_SPLINE;
  holes.push_back(fin);

  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  interpolate(false, OUTER_NODES_PER_SPLINE,
              outer_x0_cubics, outer_x1_cubics);

  std::vector<double> fin_x0_spline_points, fin_x1_spline_points;
  get_fin_spline_points(&fin_x0_spline_points, &fin_x1_spline_points);

  std::vector<std::vector<double>> fin_x0_cubics, fin_x1_cubics;
  get_cubics(fin_x0_spline_points, fin_x1_spline_points,
             &fin_x0_cubics, &fin_x1_cubics);

  interpolate(true, FIN_NODES_PER_SPLINE,
              fin_x0_cubics, fin_x1_cubics);


  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(weights.size() * 2, 1);
    apply_boundary_condition(0, weights.size(), LEFT_TO_RIGHT_FLOW);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}


}  // namespace ie_solver
