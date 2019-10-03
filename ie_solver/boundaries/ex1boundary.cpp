// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/log.h"

#define OUTER_NUM_SPLINE_POINTS 24
#define STAR_NUM_SPLINE_POINTS 8

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
  double longer = 0.05;
  double shorter = 0.02;

  x0_points->push_back(x);
  x1_points->push_back(y - longer);

  x0_points->push_back(x + shorter);
  x1_points->push_back(y - shorter);

  x0_points->push_back(x + longer);
  x1_points->push_back(y);

  x0_points->push_back(x + shorter);
  x1_points->push_back(y + shorter);

  x0_points->push_back(x);
  x1_points->push_back(y + longer);

  x0_points->push_back(x - shorter);
  x1_points->push_back(y + shorter);

  x0_points->push_back(x - longer);
  x1_points->push_back(y);

  x0_points->push_back(x - shorter);
  x1_points->push_back(y - shorter);
}


void Ex1Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::EX1;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  int STAR_NODES_PER_SPLINE = (N / 12) / STAR_NUM_SPLINE_POINTS;
  int NUM_CIRCLE_POINTS = (N / 12);
  int OUTER_NODES_PER_SPLINE = (2 * N / 3) / OUTER_NUM_SPLINE_POINTS;
  Hole star1, star2, circle1, circle2;

  if (holes.size() == 0) {
    star1.center = Vec2(0.2, 0.5);
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
    circle2.center = Vec2(0.8, 0.5);
    circle2.radius = 0.05;
    circle2.num_nodes =  NUM_CIRCLE_POINTS;
    holes.push_back(circle2);
  } else {
    star1 = holes[0];
    star2 = holes[1];
    circle1 = holes[2];
    circle2 = holes[3];
  }

  int total_num = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE +
                  2 * STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE
                  + 2 * NUM_CIRCLE_POINTS;
  boundary_values = ie_Mat(2 * total_num, 1);

  int bc_index = 0;

  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  interpolate(bc_index, false, OUTER_NODES_PER_SPLINE , TANGENT_VEC,
              outer_x0_cubics, outer_x1_cubics);
  bc_index +=  OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;
  num_outer_nodes = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;
  std::vector<double> star_x0_points, star_x1_points;
  get_star_spline_points(star1.center.a[0], star1.center.a[1], &star_x0_points,
                         &star_x1_points);

  std::vector<std::vector<double>> star_x0_cubics, star_x1_cubics;
  get_cubics(star_x0_points, star_x1_points,
             &star_x0_cubics, &star_x1_cubics);

  interpolate(bc_index, true, STAR_NODES_PER_SPLINE, REVERSE_NORMAL_VEC,
              star_x0_cubics, star_x1_cubics);
  bc_index += STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;

  std::vector<double> star2_x0_points, star2_x1_points;
  get_star_spline_points(star2.center.a[0], star2.center.a[1], &star2_x0_points,
                         &star2_x1_points);

  std::vector<std::vector<double>> star2_x0_cubics, star2_x1_cubics;
  get_cubics(star2_x0_points, star2_x1_points,
             &star2_x0_cubics, &star2_x1_cubics);

  interpolate(bc_index, true, STAR_NODES_PER_SPLINE, NORMAL_VEC,
              star2_x0_cubics, star2_x1_cubics);
  bc_index += STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;

  for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
    double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
    points.push_back(circle1.center.a[0] + circle1.radius * cos(ang));
    points.push_back(circle1.center.a[1] + circle1.radius * sin(ang));
    normals.push_back(-cos(ang));
    normals.push_back(-sin(ang));
    curvatures.push_back(-1.0 / circle1.radius);
    weights.push_back(2.0 * M_PI * circle1.radius / NUM_CIRCLE_POINTS);
    boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index ]);
    boundary_values.set(2 * bc_index + 1, 0, -normals[2 * bc_index + 1]);
    bc_index++;
  }

  for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
    double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
    points.push_back(circle2.center.a[0] + circle2.radius * cos(ang));
    points.push_back(circle2.center.a[1] + circle2.radius * sin(ang));
    normals.push_back(-cos(ang));
    normals.push_back(-sin(ang));
    curvatures.push_back(-1.0 / circle2.radius);
    weights.push_back(2.0 * M_PI * circle2.radius / NUM_CIRCLE_POINTS);
    boundary_values.set(2 * bc_index, 0, normals[2 * bc_index ]);
    boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index + 1]);
    bc_index++;
  }

}


}  // namespace ie_solver
