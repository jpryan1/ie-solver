// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/log.h"

#define OUTER_NUM_SPLINE_POINTS 24
#define OUTER_NODES_PER_SPLINE 100
#define STAR_NODES_PER_SPLINE 16
#define STAR_NUM_SPLINE_POINTS 8
#define NUM_CIRCLE_POINTS 128

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


void Ex1Boundary::get_cubics(const std::vector<double>& x0_points,
                             const std::vector<double>& x1_points,
                             std::vector<std::vector<double>>* x0_cubics,
                             std::vector<std::vector<double>>* x1_cubics) {
  int num_spline_points = x0_points.size();
  ie_Mat x0_d_vec(num_spline_points, 1);
  ie_Mat x0_y_vec(num_spline_points, 1);

  ie_Mat x1_d_vec(num_spline_points, 1);
  ie_Mat x1_y_vec(num_spline_points, 1);
  for (int i = 0; i < num_spline_points; i++) {
    double after = x0_points[(i + 1) % num_spline_points];
    double before = x0_points[(i + num_spline_points - 1) %
                              num_spline_points];
    x0_y_vec.set(i, 0,  3 * (after - before));

    after = x1_points[(i + 1) % num_spline_points];
    before = x1_points[(i + num_spline_points - 1) % num_spline_points];
    x1_y_vec.set(i, 0,  3 * (after - before));
  }

  ie_Mat A(num_spline_points, num_spline_points);

  A.set(0, num_spline_points - 1, 1.0);
  A.set(num_spline_points - 1, 0, 1.0);

  for (int i = 0; i < num_spline_points; i++) {
    A.set(i, i, 4.0);
    if (i != 0) {
      A.set(i, i - 1, 1.0);
    }
    if (i != num_spline_points - 1) {
      A.set(i, i + 1, 1.0);
    }
  }

  A.left_multiply_inverse(x0_y_vec, &x0_d_vec);
  A.left_multiply_inverse(x1_y_vec, &x1_d_vec);


  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> cubic;
    cubic.push_back(x0_points[i]);
    cubic.push_back(x0_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x0_points[i_plus_one] - x0_points[i])
                    - 2 * x0_d_vec.get(i, 0)
                    - x0_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x0_points[i] - x0_points[i_plus_one])
                    + x0_d_vec.get(i, 0)
                    + x0_d_vec.get(i_plus_one, 0));

    x0_cubics->push_back(cubic);
  }
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> cubic;
    cubic.push_back(x1_points[i]);
    cubic.push_back(x1_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x1_points[i_plus_one] - x1_points[i])
                    - 2 * x1_d_vec.get(i, 0)
                    - x1_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x1_points[i] - x1_points[i_plus_one])
                    + x1_d_vec.get(i, 0)
                    + x1_d_vec.get(i_plus_one, 0));
    x1_cubics->push_back(cubic);
  }
}


void Ex1Boundary::interpolate(int bc_index, bool is_interior,
                              int nodes_per_spline, Ex1BC bc,
                              const std::vector<std::vector<double>>& x0_cubics,
                              const std::vector<std::vector<double>>& x1_cubics) {
  // Must fill points, normals, curvatures, weights.
  // Points = duh
  // Normals = tangent of points rotated 90 deg clockwise
  // Curvatures = (x'y'' - x''y') / (x'^2 + y'^2)^1.5
  // assert positive curvature when testing pl0x
  // Weights = uggggh just estimate it
  int start = bc_index;
  int num_spline_points = x0_cubics.size();
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> x_cubic = x0_cubics[i];
    std::vector<double> y_cubic = x1_cubics[i];
    for (int j = 0; j < nodes_per_spline; j++) {
      double t = j / (nodes_per_spline + 0.0);
      double x = x_cubic[0] + t * x_cubic[1] + pow(t, 2) * x_cubic[2]
                 + pow(t, 3) * x_cubic[3];
      double y = y_cubic[0] + t * y_cubic[1] + pow(t, 2) * y_cubic[2]
                 + pow(t, 3) * y_cubic[3];

      points.push_back(x);
      points.push_back(y);

      double x_prime = x_cubic[1] + 2 * t * x_cubic[2]
                       + 3 * pow(t, 2) * x_cubic[3];
      double y_prime = y_cubic[1] + 2 * t * y_cubic[2]
                       + 3 * pow(t, 2) * y_cubic[3];

      double x_prime_prime = 2 * x_cubic[2] + 6 * t * x_cubic[3];
      double y_prime_prime = 2 * y_cubic[2] + 6 * t * y_cubic[3];

      double curvature = (x_prime * y_prime_prime - x_prime_prime * y_prime)
                         / pow(sqrt(pow(x_prime, 2) + pow(y_prime, 2)), 3);

      double norm = sqrt(pow(x_prime, 2) + pow(y_prime, 2));
      x_prime /= norm;
      y_prime /= norm;
      if (is_interior) {
        curvatures.push_back(-curvature);
        normals.push_back(-y_prime);
        normals.push_back(x_prime);
      } else {
        curvatures.push_back(curvature);
        normals.push_back(y_prime);
        normals.push_back(-x_prime);
      }
      switch (bc) {
        case TANGENTIAL:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index]);
          break;
        case REVERSE_TANGENTIAL:
          boundary_values.set(2 * bc_index, 0, normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, -normals[2 * bc_index]);
          break;
        case NORMAL:
          boundary_values.set(2 * bc_index, 0, normals[2 * bc_index ]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index + 1]);
          break;
        case REVERSE_NORMAL:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index ]);
          boundary_values.set(2 * bc_index + 1, 0, -normals[2 * bc_index + 1]);
          break;
      }

      bc_index++;
    }
  }
  int end = bc_index;

  double dist1 = sqrt(pow(points[2 * start]
                          - points[2 * (end - 1)], 2)
                      + pow(points[2 * start + 1]
                            - points[2 * (end - 1) + 1], 2));

  double dist2 = sqrt(pow(points[2 * start]
                          - points[2 * (start + 1)], 2)
                      + pow(points[2 * start + 1]
                            - points[2 * (start + 1) + 1], 2));
  weights.push_back((dist1 + dist2) / 2.0);
  for (unsigned int i = start + 1; i < end - 1; i++) {
    dist1 = sqrt(pow(points[2 * i] - points[2 * (i - 1)], 2)
                 + pow(points[2 * i + 1] - points[2 * (i - 1) + 1], 2));
    dist2 = sqrt(pow(points[2 * i] - points[2 * (i + 1)], 2)
                 + pow(points[2 * i + 1] - points[2 * (i + 1) + 1], 2));
    weights.push_back((dist1 + dist2) / 2.0);
  }

  dist1 = sqrt(pow(points[2 * (end - 1)]
                   - points[2 * ((end - 1) - 1)], 2)
               + pow(points[2 * (end - 1) + 1]
                     - points[2 * ((end - 1) - 1) + 1], 2));
  dist2 = sqrt(pow(points[2 * (end - 1)]
                   - points[2 * (start)], 2)
               + pow(points[2 * (end - 1) + 1]
                     - points[2 * (start) + 1], 2));
  weights.push_back((dist1 + dist2) / 2.0);
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
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  all_cubics_x0.clear();
  all_cubics_x1.clear();
  Hole star1, star2, circle1, circle2;
  if (holes.size() == 0) {
    star1.center = Vec2(0.2, 0.5);
    star1.radius = 0.05;
    holes.push_back(star1);
    star2.center = Vec2(0.4, 0.5);
    star2.radius = 0.05;
    holes.push_back(star2);
    circle1.center = Vec2(0.6, 0.5);
    circle1.radius = 0.05;
    holes.push_back(circle1);
    circle2.center = Vec2(0.8, 0.5);
    circle2.radius = 0.05;
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

  interpolate(bc_index, false, OUTER_NODES_PER_SPLINE , TANGENTIAL,
              outer_x0_cubics, outer_x1_cubics);
  bc_index +=  OUTER_NUM_SPLINE_POINTS * OUTER_NODES_PER_SPLINE;

  std::vector<double> star_x0_points, star_x1_points;
  get_star_spline_points(star1.center.a[0], star1.center.a[1], &star_x0_points,
                         &star_x1_points);

  std::vector<std::vector<double>> star_x0_cubics, star_x1_cubics;
  get_cubics(star_x0_points, star_x1_points,
             &star_x0_cubics, &star_x1_cubics);

  for (int i = 0; i < star_x0_cubics.size(); i++) {
    all_cubics_x0.push_back(star_x0_cubics[i]);
    all_cubics_x1.push_back(star_x1_cubics[i]);
  }

  interpolate(bc_index, true, STAR_NODES_PER_SPLINE, REVERSE_NORMAL,
              star_x0_cubics, star_x1_cubics);
  bc_index += STAR_NUM_SPLINE_POINTS * STAR_NODES_PER_SPLINE;

  std::vector<double> star2_x0_points, star2_x1_points;
  get_star_spline_points(star2.center.a[0], star2.center.a[1], &star2_x0_points,
                         &star2_x1_points);

  std::vector<std::vector<double>> star2_x0_cubics, star2_x1_cubics;
  get_cubics(star2_x0_points, star2_x1_points,
             &star2_x0_cubics, &star2_x1_cubics);

  for (int i = 0; i < star2_x0_cubics.size(); i++) {
    all_cubics_x0.push_back(star2_x0_cubics[i]);
    all_cubics_x1.push_back(star2_x1_cubics[i]);
  }

  interpolate(bc_index, true, STAR_NODES_PER_SPLINE, NORMAL,
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


void Ex1Boundary::find_real_roots_of_cubic(const std::vector<double>& y_cubic,
    std::vector<double>* t_vals) {
  ie_Mat companion(3, 3);
  companion.set(1, 0, 1.0);
  companion.set(2, 1, 1.0);
  companion.set(0, 2, -y_cubic[0] / y_cubic[3]);
  companion.set(1, 2, -y_cubic[1] / y_cubic[3]);
  companion.set(2, 2, -y_cubic[2] / y_cubic[3]);
  *t_vals = companion.real_eigenvalues();
}


int Ex1Boundary::num_right_intersections(double x, double y, int index) {
  // Find intersection on spline, if no t in [0,1] return 0
  // For each t in that range, find corresponding x val, check if
  // greater than x.
  std::vector<double> x_cubic = all_cubics_x0[index];
  std::vector<double> y_cubic = all_cubics_x1[index];
  y_cubic[0] -= y;
  std::vector<double> t_vals;
  find_real_roots_of_cubic(y_cubic, &t_vals);
  int intersections = 0;
  for (double t : t_vals) {
    if (t > 0 && t < 1) {
      double dif = x_cubic[0] + t * x_cubic[1] + pow(t, 2) * x_cubic[2]
                   + pow(t, 3) * x_cubic[3] - x;
      if (fabs(dif) < 0.01) return -1;
      if (dif > 0) {
        intersections++;
      }
    }
  }
  return intersections;
}


bool Ex1Boundary::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  int intersections = 0;
  for (int i = 0; i < OUTER_NUM_SPLINE_POINTS; i++) {
    int right_intersections = num_right_intersections(v[0], v[1], i);
    if (right_intersections == -1) return false;
    intersections += right_intersections;
  }
  if (intersections % 2 == 0) {
    return false;
  }
  for (Hole hole : holes) {
    Vec2 r = a - hole.center;
    if (r.norm() < hole.radius + 1e-2) {
      return false;
    }
  }


  return true;
}

}  // namespace ie_solver
