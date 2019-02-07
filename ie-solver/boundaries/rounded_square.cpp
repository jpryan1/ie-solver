// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/boundaries/rounded_square.h"
#include "ie-solver/log.h"


namespace ie_solver {


void RoundedSquare::draw_line(int bc_index, int num_points,
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
    double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
    switch (boundary_condition) {
      case BoundaryCondition::SINGLE_ELECTRON:
        boundary_values.set(bc_index++, 0, potential);
        break;
      case BoundaryCondition::ALL_ONES:
        boundary_values.set(bc_index++, 0, 1.0);
        break;
    }
  }
}


void RoundedSquare::draw_quarter_circle(int bc_index, int num_points,
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
    double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
    switch (boundary_condition) {
      case BoundaryCondition::SINGLE_ELECTRON:
        boundary_values.set(bc_index++, 0, potential);
        break;
      case BoundaryCondition::ALL_ONES:
        boundary_values.set(bc_index++, 0, 1.0);
        break;
    }
  }
}


void RoundedSquare::initialize(int N, BoundaryCondition bc) {
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  boundary_values = ie_Mat(N, 1);

  // This square will have side length 0.9 and will have BL corner 0.05, 0.05

  // Sides will go from 0.1 to 0.9, rounded corners will have the rest

  // radius of rounded corner is 0.05

  int NUM_SIDE_POINTS = 8 * (N / 36);
  int NUM_CORN_POINTS = N / 36;

  double middie = (0.8 / NUM_SIDE_POINTS) * 0.5 + (2 * M_PI * 0.05 /
                  (4 * NUM_CORN_POINTS)) * 0.5;
  int bc_idx = 0;

  // bottom side
  weights.push_back(middie);
  draw_line(bc_idx, NUM_SIDE_POINTS, 0.9, 0.05, 0.1, 0.05, true);
  bc_idx += NUM_SIDE_POINTS;

  // bottom left corner
  weights.push_back(middie);
  draw_quarter_circle(bc_idx, NUM_CORN_POINTS, 0.1, 0.05, 0.05, 0.1, true);
  bc_idx += NUM_CORN_POINTS;

  // left side
  weights.push_back(middie);
  draw_line(bc_idx, NUM_SIDE_POINTS, 0.05, 0.1, 0.05, 0.9, true);
  bc_idx += NUM_SIDE_POINTS;

  // top left corner
  weights.push_back(middie);
  draw_quarter_circle(bc_idx, NUM_CORN_POINTS, 0.05, 0.9, 0.1, 0.95, true);
  bc_idx += NUM_CORN_POINTS;

  // top side
  weights.push_back(middie);
  draw_line(bc_idx, NUM_SIDE_POINTS, 0.1, 0.95, 0.9, 0.95, true);
  bc_idx += NUM_SIDE_POINTS;

  // top right corner
  weights.push_back(middie);
  draw_quarter_circle(bc_idx, NUM_CORN_POINTS, 0.9, 0.95, 0.95, 0.9, true);
  bc_idx += NUM_CORN_POINTS;

  // right side
  weights.push_back(middie);
  draw_line(bc_idx, NUM_SIDE_POINTS, 0.95, 0.9, 0.95, 0.1, true);
  bc_idx += NUM_SIDE_POINTS;

  // bottom right corner
  weights.push_back(middie);
  draw_quarter_circle(bc_idx, NUM_CORN_POINTS, 0.95, 0.1, 0.9, 0.05, true);
  bc_idx += NUM_CORN_POINTS;
}


bool RoundedSquare::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  double eps = 1e-2;

  if (fabs(v[0] - 0.5) > 0.45 - eps
      || fabs(v[1] - 0.5) > 0.45 - eps) {
    return false;
  }

  if (v[0] <= 0.95 + eps && v[0] + eps >= 0.05 && v[1] + eps >= 0.1
      && v[1] <= 0.9 + eps) return true;
  if (v[0] <= 0.9  + eps && v[0] + eps >= 0.1 && v[1] + eps >= 0.05
      && v[1] <= 0.95 + eps) return true;

  double min = 1;
  Vec2 bl(0.1, 0.1);
  Vec2 br(0.9, 0.1);
  Vec2 tl(0.1, 0.9);
  Vec2 tr(0.9, 0.9);

  min = fmin(min, (bl - a).norm());
  min = fmin(min, (tr - a).norm());
  min = fmin(min, (tl - a).norm());
  min = fmin(min, (br - a).norm());
  if (min + eps > 0.05) return false;
  return true;
}

}  // namespace ie_solver
