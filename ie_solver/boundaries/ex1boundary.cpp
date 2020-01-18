// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/log.h"

namespace ie_solver {


void Ex1Boundary::get_spline_points(std::vector<double>* x0_points,
                                    std::vector<double>* x1_points) {
  for (int i = 0; i < 3; i++) {
    x0_points->push_back(i / 3.0);
    x1_points->push_back(0);

    x0_points->push_back((i / 3.0) + (1.0 / 6.0));
    x1_points->push_back(0.05);
  }

  for (int i = 0; i < 3; i++) {
    x0_points->push_back(1.0);
    x1_points->push_back(i / 3.0);

    x0_points->push_back(0.95);
    x1_points->push_back((i / 3.0) + (1.0 / 6.0));
  }

  for (int i = 0; i < 3; i++) {
    x0_points->push_back(1.0 - (i / 3.0));
    x1_points->push_back(1.0);

    x0_points->push_back(1.0 - (1.0 / 6.0) - (i / 3.0));
    x1_points->push_back(0.95);
  }

  for (int i = 0; i < 3; i++) {
    x0_points->push_back(0);
    x1_points->push_back(1.0 - (i / 3.0));

    x0_points->push_back(0.05);
    x1_points->push_back(1.0 - (1.0 / 6.0) - (i / 3.0));
  }
}


void Ex1Boundary::get_star_spline_points(double x, double y,
    std::vector<double>* x0_points, std::vector<double>* x1_points) {

  for (int i = 0; i < 20; i++) {
    double ang = 2 * M_PI * (i / (20.));

    double xc =  0.015 * cos(ang) * (sin(5 * ang) + 4);
    double yc =  0.015 * sin(ang) * (sin(5 * ang) + 4);

    x0_points->push_back(x + xc);
    x1_points->push_back(y + yc);
  }
}


void Ex1Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::EX1;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  holes.clear();

  if (perturbation_parameters.size() == 0) {
    perturbation_parameters.push_back(M_PI / 2.);
  }

  int OUTER_NUM_SPLINE_POINTS = 24;
  int STAR_NUM_SPLINE_POINTS = 20;

  int OUTER_NODES_PER_SPLINE = (2 * N / 3) / OUTER_NUM_SPLINE_POINTS;
  int STAR_NODES_PER_SPLINE = (N / 12) / STAR_NUM_SPLINE_POINTS;

  int NUM_CIRCLE_POINTS = (N / 12);

  num_outer_nodes = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;

  Hole star1, star2, circle1, circle2;

  star1.center = Vec2(0.5 + 0.3 * cos(M_PI + perturbation_parameters[0]),
                      0.5 + 0.3 * sin(M_PI + perturbation_parameters[0]));
  star1.radius = 0.05;
  star1.num_nodes =  STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;
  holes.push_back(star1);
  star2.center = Vec2(0.4, 0.5);
  star2.radius = 0.05;
  star2.num_nodes =  STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;
  holes.push_back(star2);
  circle1.center = Vec2(0.6, 0.5);
  circle1.radius = 0.05;
  circle1.num_nodes =  NUM_CIRCLE_POINTS;
  holes.push_back(circle1);
  circle2.center = Vec2(0.5 + 0.3 * cos(perturbation_parameters[0]),
                        0.5 + 0.3 * sin(perturbation_parameters[0]));
  circle2.radius = 0.05;
  circle2.num_nodes =  NUM_CIRCLE_POINTS;
  holes.push_back(circle2);

  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  interpolate(false, OUTER_NODES_PER_SPLINE, outer_x0_cubics, outer_x1_cubics);
  std::vector<double> star_x0_points, star_x1_points;
  get_star_spline_points(star1.center.a[0], star1.center.a[1], &star_x0_points,
                         &star_x1_points);

  std::vector<std::vector<double>> star_x0_cubics, star_x1_cubics;
  get_cubics(star_x0_points, star_x1_points,
             &star_x0_cubics, &star_x1_cubics);

  interpolate(true, STAR_NODES_PER_SPLINE, star_x0_cubics, star_x1_cubics);

  std::vector<double> star2_x0_points, star2_x1_points;
  get_star_spline_points(star2.center.a[0], star2.center.a[1], &star2_x0_points,
                         &star2_x1_points);

  std::vector<std::vector<double>> star2_x0_cubics, star2_x1_cubics;
  get_cubics(star2_x0_points, star2_x1_points,
             &star2_x0_cubics, &star2_x1_cubics);

  interpolate(true, STAR_NODES_PER_SPLINE, star2_x0_cubics, star2_x1_cubics);

  for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
    double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
    points.push_back(circle1.center.a[0] + circle1.radius * cos(ang));
    points.push_back(circle1.center.a[1] + circle1.radius * sin(ang));
    normals.push_back(-cos(ang));
    normals.push_back(-sin(ang));
    curvatures.push_back(-1.0 / circle1.radius);
    weights.push_back(2.0 * M_PI * circle1.radius / NUM_CIRCLE_POINTS);
  }

  for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
    double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
    points.push_back(circle2.center.a[0] + circle2.radius * cos(ang));
    points.push_back(circle2.center.a[1] + circle2.radius * sin(ang));
    normals.push_back(-cos(ang));
    normals.push_back(-sin(ang));
    curvatures.push_back(-1.0 / circle2.radius);
    weights.push_back(2.0 * M_PI * circle2.radius / NUM_CIRCLE_POINTS);
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(weights.size() * 2, 1);
    int b1 = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;
    int b2 = b1 + STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;
    int b3 = b2 + STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;
    int b4 = b3 + NUM_CIRCLE_POINTS;
    int b5 = b4 + NUM_CIRCLE_POINTS;
    apply_boundary_condition(0, b1, TANGENT_VEC);
    apply_boundary_condition(b1, b2, REVERSE_NORMAL_VEC);
    apply_boundary_condition(b2, b3, NORMAL_VEC);
    apply_boundary_condition(b3, b4, REVERSE_NORMAL_VEC);
    apply_boundary_condition(b4, b5, NORMAL_VEC);

  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}


}  // namespace ie_solver
