// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex3boundary.h"
#include "ie_solver/log.h"

#define OUTER_NUM_SPLINE_POINTS 12
#define OUTER_NODES_PER_SPLINE 100
#define FIN_NODES_PER_SPLINE 32
#define FIN_NUM_SPLINE_POINTS 8
#define NUM_CIRCLE_POINTS 128
#define FIN_RAD ((0.2*sqrt(2)-0.075)/2.0)
namespace ie_solver {



// USE ROUNDED SQUARE FOR EXTERIOR BOUNDARY

// move boundaries outwards to make laminar flow
void Ex3Boundary::get_spline_points(std::vector<double>* x0_points,
                                    std::vector<double>* x1_points) {
  for (int i = 0; i < 4; i++) {
    x0_points->push_back(i / 4.0);
    x1_points->push_back(0.25);
  }
  for (int i = 0; i < 2; i++) {
    x0_points->push_back(1);
    x1_points->push_back(0.25 + 0.5 * i / 2.0);
  }

  for (int i = 0; i < 4; i++) {
    x0_points->push_back(1 - i / 4.0);
    x1_points->push_back(0.75);
  }
  for (int i = 0; i < 2; i++) {
    x0_points->push_back(0);
    x1_points->push_back(0.75 - 0.5 * i / 2.0);
  }
}


void Ex3Boundary::get_star_spline_points(double x, double y,
    std::vector<double>* x0_points, std::vector<double>* x1_points) {
  double longer = 0.025;
  double shorter = 0.01;

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


void Ex3Boundary::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  all_cubics_x0.clear();
  all_cubics_x1.clear();

  Hole circle1, circle2, circle3, circle4, circle5, circle6, circle7;//fin;
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


  circle6.center = Vec2(0.5, 0.5);
  circle6.radius = 0.075;
  holes.push_back(circle6);

  // fin.center = Vec2(0.5, 0.5);
  // fin.radius = FIN_RAD;
  // holes.push_back(fin);


  int total_num = OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE +
                  + 7 * NUM_CIRCLE_POINTS;
  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * total_num, 1);
  } else {
    boundary_values = ie_Mat(total_num, 1);
  }
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

      boundary_values.set(2 * bc_index, 0, 0);
      boundary_values.set(2 * bc_index + 1, 0, 0);

      bc_index++;
    }
  }
}


bool Ex3Boundary::is_in_domain(const Vec2& a) {
  const double* v = a.a;


// hard code rectangle inclusion
  int intersections = 0;
  for (int i = 0; i < all_cubics_x1.size(); i++) {
    int right_intersections = num_right_intersections(v[0], v[1], i);
    if (right_intersections == -1) return false;
    intersections += right_intersections;
  }
  if (intersections % 2 == 0) {
    return false;
  }

  for (Hole hole : holes) {
    if(hole.radius == FIN_RAD) continue;
    Vec2 r = a - hole.center;
    if (r.norm() < hole.radius + 1e-2) {
      return false;
    }
  }
  return true;
}

}  // namespace ie_solver
