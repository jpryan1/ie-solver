// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex3boundary.h"
#include "ie_solver/log.h"

#define OUTER_NUM_SPLINE_POINTS 28
#define FIN_SPLINE_POINTS 4
#define FIN_RAD 0.075
namespace ie_solver {


// USE ROUNDED SQUARE FOR EXTERIOR BOUNDARY

// move boundaries outwards to make laminar flow
void Ex3Boundary::get_spline_points(std::vector<double>* x0_points,
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

void Ex3Boundary::get_fin_spline_points(std::vector<double>* x0_points,
                                        std::vector<double>* x1_points) {
  x0_points->push_back(0.025);
  x1_points->push_back(0.);

  x0_points->push_back(0.);
  x1_points->push_back(0.075);

  x0_points->push_back(-0.025);
  x1_points->push_back(0.);

  x0_points->push_back(0.);
  x1_points->push_back(-0.075);

  // Rotate by perturbation_parameter

  for (int i = x0_points->size() - 4; i < x0_points->size(); i++) {
    double temp = cos(perturbation_parameter) * (*x0_points)[i] - sin(
                    perturbation_parameter) * (*x1_points)[i];
    (*x1_points)[i] = 0.5 + sin(perturbation_parameter) * (*x0_points)[i] + cos(
                        perturbation_parameter) * (*x1_points)[i];
    (*x0_points)[i] = 0.5 + temp;
  }
}


void Ex3Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  all_cubics_x0.clear();
  all_cubics_x1.clear();

  //
  int OUTER_NODES_PER_SPLINE = (3 * N / 4) / OUTER_NUM_SPLINE_POINTS;
  int NUM_CIRCLE_POINTS = (N / 4) / 8;
  int FIN_NODES_PER_SPLINE = ((N / 4) / 4) / FIN_SPLINE_POINTS;
  if (holes.size() == 0) {
    Hole circle1, circle2, circle3, circle4, circle5, circle6, fin;
    circle1.center = Vec2(0.2, 0.4);
    circle1.radius = 0.025;
    holes.push_back(circle1);
    circle2.center = Vec2(0.3, 0.4);
    circle2.radius = 0.025;
    holes.push_back(circle2);

    circle3.center = Vec2(0.4, 0.4);
    circle3.radius = 0.025;
    holes.push_back(circle3);
    circle4.center = Vec2(0.6, 0.6);
    circle4.radius = 0.025;
    holes.push_back(circle4);

    circle5.center = Vec2(0.7, 0.6);
    circle5.radius = 0.025;
    holes.push_back(circle5);
    circle6.center = Vec2(0.8, 0.6);
    circle6.radius = 0.025;
    holes.push_back(circle6);
    fin.center = Vec2(0.5, 0.5);
    fin.radius = FIN_RAD;
    holes.push_back(fin);
  }

  int total_num = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE +
                  + 6 * NUM_CIRCLE_POINTS + FIN_SPLINE_POINTS * FIN_NODES_PER_SPLINE;
  boundary_values = ie_Mat(2 * total_num, 1);

  int bc_index = 0;

  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  for (int i = 0; i < outer_x0_cubics.size(); i++) {
    all_cubics_x0.push_back(outer_x0_cubics[i]);
    all_cubics_x1.push_back(outer_x1_cubics[i]);
  }

  interpolate(bc_index, false, OUTER_NODES_PER_SPLINE, LEFT_TO_RIGHT_FLOW,
              outer_x0_cubics, outer_x1_cubics);
  bc_index +=  OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;
  num_outer_nodes = OUTER_NODES_PER_SPLINE * OUTER_NUM_SPLINE_POINTS;

  std::vector<double> fin_x0_spline_points, fin_x1_spline_points;
  get_fin_spline_points(&fin_x0_spline_points, &fin_x1_spline_points);

  std::vector<std::vector<double>> fin_x0_cubics, fin_x1_cubics;
  get_cubics(fin_x0_spline_points, fin_x1_spline_points,
             &fin_x0_cubics, &fin_x1_cubics);

  for (int i = 0; i < fin_x0_cubics.size(); i++) {
    all_cubics_x0.push_back(fin_x0_cubics[i]);
    all_cubics_x1.push_back(fin_x1_cubics[i]);
  }

  interpolate(bc_index, true, FIN_NODES_PER_SPLINE, NO_SLIP,
              fin_x0_cubics, fin_x1_cubics);
  bc_index +=  fin_x0_cubics.size() * FIN_NODES_PER_SPLINE;

  for (int i = 0; i < holes.size() - 1; i++) {
    Hole circle = holes[i];
    for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
      double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
      points.push_back(circle.center.a[0] + circle.radius * cos(ang));
      points.push_back(circle.center.a[1] + circle.radius * sin(ang));
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / circle.radius);
      weights.push_back(2.0 * M_PI * circle.radius / NUM_CIRCLE_POINTS);

      boundary_values.set(2 * bc_index, 0, 0);
      boundary_values.set(2 * bc_index + 1, 0, 0);

      bc_index++;
    }
  }
}


}  // namespace ie_solver
