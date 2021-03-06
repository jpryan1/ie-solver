// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie_solver/boundaries/rounded_square.h"
#include "ie_solver/log.h"


namespace ie_solver {


void RoundedSquare::draw_line(int num_points,
                              double start_x, double start_y,
                              double end_x, double end_y,
                              bool normal_is_left) {
  // NOTE: The first and last weights are not added - that is done in initialize
  // A point is placed on start_, not on end_

  double normal_x = -(end_y - start_y);
  double normal_y = end_x - start_x;
  double norm = sqrt(pow(normal_x, 2) + pow(normal_y, 2));
  normal_x /= norm;
  normal_y /= norm;
  if (!normal_is_left) {
    normal_x *= -1;
    normal_y *= -1;
  }
  double weight = norm / num_points;

  for (int i = 0; i < num_points; i++) {
    double x = start_x + (end_x - start_x) * ((i + 0.0) / num_points);
    double y = start_y + (end_y - start_y) * ((i + 0.0) / num_points);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(normal_x);
    normals.push_back(normal_y);
    curvatures.push_back(0);
    if (i != 0) {
      weights.push_back(weight);
    }
  }
}


void RoundedSquare::draw_quarter_circle(int num_points,
                                        double start_x, double start_y,
                                        double end_x, double end_y,
                                        bool convex) {
  // NOTE: The first and last weights are not added - that is done in initialize
  // A point is placed on start_, not on end_
  // From start_ to end_, a clockwise quartercircle is drawn. If convex, then
  // normals point out from clock, else in.
  // Assumed for now that the slope of the line through the points is +/- 1

  double c_x, c_y, ang;
  if (!convex) {
    double swap = end_x;
    end_x = start_x;
    start_x = swap;
    swap = end_y;
    end_y = start_y;
    start_y = swap;
  }
  if (start_x < end_x) {
    if (start_y < end_y) {
      c_x = end_x;
      c_y = start_y;
      ang = M_PI;
    } else {
      c_x = start_x;
      c_y = end_y;
      ang = M_PI / 2.0;
    }
  } else {
    if (start_y < end_y) {
      c_x = start_x;
      c_y = end_y;
      ang = 3 * M_PI / 2.0;
    } else {
      c_x = end_x;
      c_y = start_y;
      ang = 0;
    }
  }
  if (!convex) {
    ang -= M_PI / 2.0;
  }

  double rad = fabs(end_x - start_x);
  double curvature = 1.0 / rad;
  if (!convex) {
    curvature *= -1;
  }
  double weight = 2.0 * M_PI * rad / (4.0 * num_points);

  for (int i = 0; i < num_points; i++) {
    double current_ang;
    if (convex) {
      current_ang = ang - ((i * M_PI) / (2 * num_points));
    } else {
      current_ang = ang + ((i * M_PI) / (2 * num_points));
    }
    double x = c_x + rad * cos(current_ang);
    double y = c_y + rad * sin(current_ang);
    points.push_back(x);
    points.push_back(y);
    double normal_x = cos(current_ang);
    double normal_y = sin(current_ang);
    if (!convex) {
      normal_x *= -1;
      normal_y *= -1;
    }
    normals.push_back(normal_x);
    normals.push_back(normal_y);
    curvatures.push_back(curvature);
    if (i != 0) {
      weights.push_back(weight);
    }
  }
}


void RoundedSquare::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::ROUNDED_SQUARE;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
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

  N = 60 * (N / 60);
  int SCALE_UNIT = N / 60;
  int line_points = 2 * SCALE_UNIT;
  int corner_points = 3 * SCALE_UNIT;

  // Line width = 0.1
  // Corner radius = 0.1
  double line_weight = 0.1 / line_points;
  double corner_weight = (2.0 * M_PI * 0.1 / 4.0) / corner_points;
  double middie = (line_weight + corner_weight) / 2.0;


  weights.push_back(middie);
  draw_line(6 * line_points, 0.8, 0.1, 0.2, 0.1, true);

  weights.push_back(middie);
  draw_quarter_circle(corner_points, 0.2, 0.1, 0.1, 0.2, true);

  weights.push_back(middie);
  draw_line(6 * line_points, 0.1, 0.2, 0.1, 0.8, true);

  weights.push_back(middie);
  draw_quarter_circle(corner_points, 0.1, 0.8, 0.2, 0.9, true);

  weights.push_back(middie);
  draw_line(6 * line_points, 0.2, 0.9, 0.8, 0.9, true);

  weights.push_back(middie);
  draw_quarter_circle(corner_points, 0.8, 0.9, 0.9, 0.8, true);

  weights.push_back(middie);
  draw_line(6 * line_points, 0.9, 0.8, 0.9, 0.2, true);

  weights.push_back(middie);
  draw_quarter_circle(corner_points, 0.9, 0.2, 0.8, 0.1, true);

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ie_Mat(weights.size(), 1);
    apply_boundary_condition(0, weights.size(), SINGLE_ELECTRON);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}


bool RoundedSquare::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  double eps = 1e-2;

  if (fabs(v[0] - 0.5) > 0.4 - eps
      || fabs(v[1] - 0.5) > 0.4 - eps) {
    return false;
  }

  if (v[0] <= 0.9 + eps && v[0] + eps >= 0.1 && v[1] + eps >= 0.2
      && v[1] <= 0.8 + eps) return true;
  if (v[0] <= 0.8  + eps && v[0] + eps >= 0.2 && v[1] + eps >= 0.1
      && v[1] <= 0.9 + eps) return true;

  double min = 1;
  Vec2 bl(0.2, 0.2);
  Vec2 br(0.8, 0.2);
  Vec2 tl(0.2, 0.8);
  Vec2 tr(0.8, 0.8);

  min = fmin(min, (bl - a).norm());
  min = fmin(min, (tr - a).norm());
  min = fmin(min, (tl - a).norm());
  min = fmin(min, (br - a).norm());
  if (min + eps > 0.1) return false;
  return true;
}

}  // namespace ie_solver
