// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/boundaries/rounded_square_with_bump.h"
#include "ie-solver/log.h"


namespace ie_solver {
void RoundedSquareWithBump::reinitialize(int N, int bc_enum, int bump) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  bump_size = bump;
  initialize(N, bc_enum);
}


void RoundedSquareWithBump::draw_line(int bc_index, int num_points,
                                      double start_x, double start_y,
                                      double end_x, double end_y,
                                      bool normal_is_left, int bc_enum) {
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
    boundary_condition.set(bc_index++, 0, potential);
  }
}


void RoundedSquareWithBump::draw_quarter_circle(int bc_index, int num_points,
    double start_x, double start_y,
    double end_x, double end_y, bool convex,
    int bc_enum) {
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
    boundary_condition.set(bc_index++, 0, potential);
  }
}


void RoundedSquareWithBump::initialize(int N, int bc_enum) {
  if (bc_enum != BoundaryCondition::SINGLE_ELECTRON) {
    LOG::ERROR("Circle boundary can only do single electron bc currently.");
  }

  // For now we ignore N and force the number of discretization points.
  int line_points = 128;
  int corner_points = 192;  // = 128*1.5
  double line_weight = 0.1 / line_points;
  double corner_weight = (M_PI / 20.0) / corner_points;
  double middie = (line_weight + corner_weight) / 2.0;
  // The idea is that line_weight and corner_weight should hopefully be very
  // similar.
  N = 22 * line_points + 8 * corner_points;

  boundary_condition = ie_Mat(N, 1);

  int bc_index = 0;


  // line ed

  weights.push_back(middie);
  draw_line(bc_index, 6 * line_points, 0.8, 0.1, 0.2, 0.1, true, 0);
  bc_index += 6 * line_points;

  // corner d

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.2, 0.1, 0.1, 0.2, true, 0);
  bc_index += corner_points;

  // line dc

  weights.push_back(middie);
  draw_line(bc_index, line_points, 0.1, 0.2, 0.1, 0.3, true, 0);
  bc_index += line_points;

  // corner c

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.1, 0.3, 0.2, 0.4, true, 0);
  bc_index += corner_points;

  // line co
  double o_x = 0.2 + ((0.0 + bump_size) / line_points) * 0.1;
  weights.push_back(middie);
  draw_line(bc_index, bump_size, 0.2, 0.4, o_x, 0.4, true, 0);
  bc_index += line_points;

  // corner o1
  // Careful here, not middie between two quarter circles
  weights.push_back(corner_weight);
  draw_quarter_circle(bc_index, corner_points, o_x, 0.4, o_x + 0.1, 0.5, false,
                      0);
  bc_index += corner_points;

  // corner o2

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, o_x + 0.1, 0.5, o_x, 0.6, false,
                      0);
  bc_index += corner_points;

  // line ob
  weights.push_back(middie);
  draw_line(bc_index, line_points, 0.3, 0.6, 0.2, 0.6, true, 0);
  bc_index += line_points;

  // corner b

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.2, 0.6, 0.1, 0.7, true, 0);
  bc_index += corner_points;

  // line ba

  weights.push_back(middie);
  draw_line(bc_index, line_points, 0.1, 0.7, 0.1, 0.8, true, 0);
  bc_index += line_points;

  // corner a

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.1, 0.8, 0.2, 0.9, true, 0);
  bc_index += corner_points;

  // line af

  weights.push_back(middie);
  draw_line(bc_index, 6 * line_points, 0.2, 0.9, 0.8, 0.9, true, 0);
  bc_index += 6 * line_points;

  // corner f

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.8, 0.9, 0.9, 0.8, true, 0);
  bc_index += corner_points;

  // line fe

  weights.push_back(middie);
  draw_line(bc_index, 6 * line_points, 0.9, 0.8, 0.9, 0.2, true, 0);
  bc_index += 6 * line_points;

  // corner e

  weights.push_back(middie);
  draw_quarter_circle(bc_index, corner_points, 0.9, 0.2, 0.8, 0.1, true, 0);
  bc_index += corner_points;

  assert(weights.size() == N);
}

// TODO(John) consider generalizing this somehow.
bool RoundedSquareWithBump::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  double eps = 1e-2;

  double o_dist = sqrt(pow(v[0] - 0.4, 2) + pow(v[1] - 0.5, 2));
  if (o_dist - eps < 0.1) {
    return false;
  }
  if (v[0] - eps < 0.4 && v[1] - eps < 0.6 && v[1] + eps > 0.4) {
    return false;
  }
  if (v[0] + eps < 0.8 && v[0] - eps > 0.2 && v[1] + eps < 0.8
      && v[1] - eps > 0.2) {
    return true;
  }
  if (v[0] + eps > 0.9 || v[0] - eps < 0.1 || v[1] + eps > 0.9
      || v[1] - eps < 0.1) {
    return false;
  }
  // Now make sure not in any of six scraggly corners

  double f_dist = sqrt(pow(v[0] - 0.8, 2) + pow(v[1] - 0.8, 2));
  if (v[0] + eps > 0.8 && v[1] + eps > 0.8 && f_dist + eps > 0.1) {
    return false;
  }

  double e_dist = sqrt(pow(v[0] - 0.8, 2) + pow(v[1] - 0.2, 2));
  if (v[0] + eps > 0.8 && v[1] - eps < 0.2 && e_dist + eps > 0.1) {
    return false;
  }

  double d_dist = sqrt(pow(v[0] - 0.2, 2) + pow(v[1] - 0.2, 2));
  if (v[0] - eps < 0.2 && v[1] - eps < 0.2 && d_dist + eps > 0.1) {
    return false;
  }


  double a_dist = sqrt(pow(v[0] - 0.2, 2) + pow(v[1] - 0.8, 2));
  if (v[0] - eps < 0.2 && v[1] + eps > 0.8 && a_dist + eps > 0.1) {
    return false;
  }

  double b_dist = sqrt(pow(v[0] - 0.2, 2) + pow(v[1] - 0.7, 2));
  if (v[0] - eps < 0.2 && v[1] - eps < 0.7 && v[1] + eps > 0.6
      && b_dist + eps > 0.1) {
    return false;
  }
  double c_dist = sqrt(pow(v[0] - 0.2, 2) + pow(v[1] - 0.3, 2));
  if (v[0] - eps < 0.2 && v[1] - eps < 0.4 && v[1] + eps > 0.3
      && c_dist + eps > 0.1) {
    return false;
  }

  return true;
}

}  // namespace ie_solver
