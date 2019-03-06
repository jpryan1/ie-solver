// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "ie-solver/boundaries/cubic_spline.h"
#include "ie-solver/log.h"


namespace ie_solver {


void CubicSpline::get_spline_points() {
  for (int i = 0; i < 10; i++) {
    double ang = 2 * M_PI * (i / 10.0);

    double x = 0.5 + 0.2 * cos(ang);
    double y = 0.5 + 0.2 * sin(ang);
    if (i == 0) x += 0.2;
    if (i == 3) y += 0.2;
    if (i == 5) x -= 0.2;
    if (i == 8) y -= 0.2;

    x0_spline_points.push_back(x);
    x1_spline_points.push_back(y);
  }
}


void CubicSpline::get_cubics() {
  ie_Mat x0_d_vec(num_spline_points, 1);
  ie_Mat x0_y_vec(num_spline_points, 1);

  ie_Mat x1_d_vec(num_spline_points, 1);
  ie_Mat x1_y_vec(num_spline_points, 1);

  for (int i = 0; i < num_spline_points; i++) {
    double after = x0_spline_points[(i + 1) % num_spline_points];
    double before = x0_spline_points[(i + num_spline_points - 1) %
                                     num_spline_points];
    x0_y_vec.set(i, 0,  3 * (after - before));

    after = x1_spline_points[(i + 1) % num_spline_points];
    before = x1_spline_points[(i + num_spline_points - 1) % num_spline_points];
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
    cubic.push_back(x0_spline_points[i]);
    cubic.push_back(x0_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x0_spline_points[i_plus_one] - x0_spline_points[i])
                    - 2 * x0_d_vec.get(i, 0)
                    - x0_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x0_spline_points[i] - x0_spline_points[i_plus_one])
                    + x0_d_vec.get(i, 0)
                    + x0_d_vec.get(i_plus_one, 0));

    x0_cubics.push_back(cubic);
  }

  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> cubic;
    cubic.push_back(x1_spline_points[i]);
    cubic.push_back(x1_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x1_spline_points[i_plus_one] - x1_spline_points[i])
                    - 2 * x1_d_vec.get(i, 0)
                    - x1_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x1_spline_points[i] - x1_spline_points[i_plus_one])
                    + x1_d_vec.get(i, 0)
                    + x1_d_vec.get(i_plus_one, 0));
    x1_cubics.push_back(cubic);
  }
}


void CubicSpline::interpolate() {
  // Must fill points, normals, curvatures, weights.
  // Points = duh
  // Normals = tangent of points rotated 90 deg clockwise
  // Curvatures = (x'y'' - x''y') / (x'^2 + y'^2)^1.5
  // assert positive curvature when testing pl0x
  // Weights = uggggh just estimate it
  int bc_index = 0;
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> x_cubic = x0_cubics[i];
    std::vector<double> y_cubic = x1_cubics[i];
    for (int j = 0; j < 100; j++) {
      double t = j / 100.0;
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
      curvatures.push_back(curvature);

      double norm = sqrt(pow(x_prime, 2) + pow(y_prime, 2));
      x_prime /= norm;
      y_prime /= norm;
      normals.push_back(y_prime);
      normals.push_back(-x_prime);

      double potential = log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI);
      switch (boundary_condition) {
        case BoundaryCondition::SINGLE_ELECTRON:
          boundary_values.set(bc_index++, 0, potential);
          break;
        case BoundaryCondition::ALL_ONES:
          boundary_values.set(bc_index++, 0, 1.0);
          break;
        case BoundaryCondition::BUMP_FUNCTION: {
          double N = boundary_values.height();
          double x_val = -1 * ((N - 1.0 - bc_index) / (N - 1.0))
                         + (bc_index / (N - 1.0));
          potential = exp(-1.0 / (1.0 - pow(x_val, 2)));
          boundary_values.set(bc_index++, 0, potential);
          break;
        }
        case BoundaryCondition::STOKES:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index]);
          bc_index++;
          break;
      }
    }
  }

  for (int i = 0; i < curvatures.size(); i++) {
    int before = (i + curvatures.size() - 1) % curvatures.size();
    int after = (i + 1) % curvatures.size();
    double dist1 = sqrt(pow(points[2 * i] - points[2 * before], 2)
                        + pow(points[2 * i + 1] - points[2 * before + 1], 2));
    double dist2 = sqrt(pow(points[2 * i] - points[2 * after], 2)
                        + pow(points[2 * i + 1] - points[2 * after + 1], 2));
    weights.push_back((dist1 + dist2) / 2.0);
  }
}


void CubicSpline::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::CUBIC_SPLINE;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  // Currently, just get a working interpolation going
  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * 1000, 1);
  } else {
    boundary_values = ie_Mat(1000, 1);
  }  get_spline_points();
  get_cubics();
  interpolate();
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
}

void CubicSpline::find_real_roots_of_cubic(const std::vector<double>& y_cubic,
    std::vector<double>* t_vals) {
  ie_Mat companion(3, 3);
  companion.set(1, 0, 1.0);
  companion.set(2, 1, 1.0);
  companion.set(0, 2, -y_cubic[0] / y_cubic[3]);
  companion.set(1, 2, -y_cubic[1] / y_cubic[3]);
  companion.set(2, 2, -y_cubic[2] / y_cubic[3]);
  *t_vals = companion.real_eigenvalues();

}


int CubicSpline::num_right_intersections(double x, double y, int index) {
  // Find intersection on spline, if no t in [0,1] return 0
  // For each t in that range, find corresponding x val, check if
  // greater than x.
  std::vector<double> x_cubic = x0_cubics[index];
  std::vector<double> y_cubic = x1_cubics[index];
  y_cubic[0] -= y;
  std::vector<double> t_vals;
  find_real_roots_of_cubic(y_cubic, &t_vals);
  int intersections = 0;
  for (double t : t_vals) {
    if (t > 0 && t < 1) {
      if (x_cubic[0] + t * x_cubic[1] + pow(t, 2)*x_cubic[2]
          + pow(t, 3)*x_cubic[3] > x) {
        intersections++;
      }
    }
  }
  return intersections;
}


bool CubicSpline::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  int intersections = 0;
  for (int i = 0; i < num_spline_points; i++) {
    intersections += num_right_intersections(v[0], v[1], i);
  }
  if (intersections % 2 == 0) {
    return false;
  }

  return true;
}

}  // namespace ie_solver
